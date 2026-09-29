use super::*;
/// **What one z-line says to one cell** — the horizontal line's answer, as the cell reads it.
///
/// The plane arrangement wrote a label on every circle a cylinder leaves on a ⊥ class, whether or
/// not a face was kept there (`emit_faces`' `disk_labels`/`arc_labels`, collected outside the keep
/// filter). So a cell's chamber is not *decided* here; it is **read** off the line at either end,
/// a disk or an arc ([`Chart::read_cell`] is the one reader).
pub(crate) enum End<'a> {
    /// The circle on this line is uncut (or cut, but every arc reads the same for this cell): one
    /// label for the whole disk, whatever the sector.
    Disk(Label),
    /// The circle is cut and the cell's sector is covered by **one or more** of its arcs: the
    /// rows the band road's `arc_at` joined on. More than one means the sector spans a rim node
    /// the chart has no vertical line for (see [`Chart::arc_around`]), and then every arc must
    /// answer the same or this end says nothing — the rule the [`End::Disk`] arm of a cut circle
    /// already follows.
    Exact(Vec<&'a ArcLabel>),
    /// The circle is cut and the sector still has no arcs to read: its two rulings are one and the
    /// same chart line (the sector is the whole circle less a ruling), or its two ends snap to a
    /// single rim node, or the rim's θ order cannot be formed, or a row of the run is missing from
    /// the split's arcs. A whole-circle cell whose arcs disagree lands here too. ☑ Counted
    /// (`end_other`); a sector that merely spans several arcs is not here — it is an
    /// [`End::Exact`] run.
    Other,
    /// The circle is cut, the cell's sector runs between two **adjacent** rim nodes, and no arc
    /// covers it: the piece of the circle no face of this class has as a rim (the three
    /// quarters a fillet's quarter leaves, the half a slot's end leaves; the split drops such a
    /// piece rather than keep a phantom edge). The lateral face is **not here**: this end decides
    /// existence, not a chamber.
    Uncovered,
    /// No circle of this cylinder on this line at all: the line is a ⊥ class outside every
    /// lateral face's span (`circle_on_class` leaves a circle on every class *within* a span).
    NoCircle,
}

impl End<'_> {
    /// **What this end says about the cell's chamber**, or `None` when it says nothing — which
    /// now includes a run of arcs that do not agree. The read is [`read_bits`], the one
    /// spelling, and it lives here so the two consumers below cannot drift into two.
    ///
    /// ☑ **Measured: the run case here is not exercised by the suite.** Making this answer `None`
    /// whenever the run holds two or more arcs leaves all 340 lib tests green — because the cells
    /// a run decides today are decided **absent** by the existence read below, and an absent cell
    /// is never asked for its chamber. This stays the general rule rather than a written-out
    /// refusal because it *is* the whole-circle arm's rule (see [`End::Disk`]'s site) with the
    /// sector's own arcs in place of the circle's; narrowing it would be the second spelling.
    fn chamber(&self, side: SolidSide, above: bool) -> Option<(bool, bool)> {
        let one = |l: &Label| read_bits(l, side, above);
        match self {
            End::Disk(l) => Some(one(l)),
            End::Exact(arcs) => {
                let mut it = arcs.iter().map(|a| one(&a.label));
                let first = it.next()?;
                it.all(|b| b == first).then_some(first)
            }
            End::Other | End::Uncovered | End::NoCircle => None,
        }
    }
}

/// Why [`Chart::arc_around`] could not hand back a run.
enum RunFail {
    /// Every piece of the run is one no contribution covers: the face is not over this sector
    /// — [`End::Uncovered`].
    AllMissing,
    /// The run cannot be named — [`End::Other`].
    Other,
}

/// One cell, read.
/// `ends`, `src2_disagree` and `exist_disagree` are the census's readers; production reads
/// `chamber`, `present` and `emit` (the emitter) and pays for the rest only as a copy.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct CellRead<'a> {
    pub(crate) ends: [End<'a>; 2],
    /// `(in_own, in_other)` on the cylinder's inside, agreed by every end that could speak.
    /// `None` when no end could, or when two ends disagreed (`src2_disagree`).
    pub(crate) chamber: Option<(bool, bool)>,
    pub(crate) src2_disagree: bool,
    /// Is this lateral face here at all — the existence question, answered by the trace
    /// where an end is cut ([`face_spans`]) and by the row's span where none is.
    pub(crate) present: bool,
    /// The trace says the face is here while no row's span covers the cell — the one direction
    /// of disagreement that is a defect (the other is a hole in a spanned face).
    pub(crate) exist_disagree: bool,
    /// Would the cell be a result face: present, with a chamber, and `keep` differing across the
    /// wall. `None` when the chamber is unknown.
    pub(crate) emit: Option<bool>,
}

impl Chart {
    /// **Whether the lateral ends at station `i` on both sides** — a station two secant walls
    /// state, where the result keeps no material **between** the two walls. The chamber there
    /// is the far side of the wall that bounds the `+θ` sector ([`StationTwin::bounds_plus`]),
    /// read off that wall's label; `keep` says whether it is in the result. Kept, the material
    /// runs across the line and so does the lateral; not kept, what lies on the two sides meets on
    /// the line alone. `false` for every station one wall states; `None` when the label or the
    /// side could not be read.
    pub(crate) fn slit_at(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        i: usize,
        side: SolidSide,
        kind: BoolKind,
    ) -> Option<bool> {
        let seg = self.theta.get(i)?;
        let Some(tw) = seg.twin else {
            return Some(false);
        };
        let (label, wall) = if tw.bounds_plus {
            (tw.label?, tw.wall)
        } else {
            (seg.label?, seg.wall)
        };
        let sd = ruling_side_of(jd, k, def, wall, seg.end[0])?;
        let plus_above = crate::arrangement::plus_theta_is_above(jd, wall, sd)?;
        let (own, other) = read_bits(&label, side, !plus_above);
        Some(!keep_for(kind, side, own, other))
    }

    /// The node of ruling `i` **on** the z-line `t`, from the carried `end ↔ z` pairing — or
    /// `None` when the ruling merely passes `t` without a node there.
    fn node_on(&self, i: usize, t: Rat) -> Option<combinatorics::NodeId> {
        let s = &self.theta[i];
        (0..2).find(|&e| s.z[e] == t).map(|e| s.end[e])
    }

    /// **The run of arcs of a cut rim that covers the sector `[x, y)`** when the sector's rulings
    /// are not themselves an adjacent node pair of that rim — a wall whose *face* stops short of this
    /// line leaves a ruling on the chart but no node on the circle (☑ the corner boss: two chords
    /// end at the plate's corner inside the footprint, so the rim has two nodes while the chart
    /// has four rulings, and three sectors lie inside the long arc).
    ///
    /// The θ order is asked of `arrangement::circular_order`, the one spelling the cells' own
    /// order comes from — never of a coordinate. A ruling's station on this line is its **name**
    /// ([`Chart::station_on`]): a ruling with a piece ending here is that piece's node, one
    /// without is the same canonical name a piece would have carried — and either way the rim
    /// is searched by name equality, so no point is handed to the order twice. The run is the
    /// rim nodes from the nearest one at or before `x` to the nearest at or after `y`, walked
    /// CCW — one arc when no rim node lies strictly inside the sector, and every arc it spans
    /// when some do. `None` when that walk yields no arc at all: the sector is the whole circle
    /// less one ruling, or the order cannot be formed, or a row of the run is not among `arcs`.
    /// ★ That last one is **silent here and named at the adjacent-pair site**
    /// (`RulingBoundNotYet`): a `None` falls back to the axial span, which is the conservative
    /// road, not a wrong answer — but the two sites state the same missing row differently, and
    /// that is worth one road one day.
    #[allow(clippy::too_many_arguments)]
    fn arc_around<'a>(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        c: usize,
        t: Rat,
        x: usize,
        y: usize,
        rim: &crate::draft::CutRim,
        arcs: &'a [ArcLabel],
        aliases: &crate::arrangement::Aliases,
    ) -> Result<Vec<&'a ArcLabel>, RunFail> {
        let m = rim.nodes.len();
        if x == y && m >= 2 {
            return Err(RunFail::Other);
        }
        let mut list: Vec<combinatorics::NodeId> = rim.nodes.clone();
        let mut index_of = |i: usize| -> Option<usize> {
            // A station the rim's own arc decomposition does not hold — the wall's face stops
            // short of the rim — joins the order as itself, under its canonical name, and the
            // search below asks which arc contains it. ★ Not under a node of *another* line
            // ([`Chart::station_on`]'s note): the order refuses that wherever the rim does hold
            // the station.
            let n = match self.node_on(i, t) {
                Some(n) => n,
                None => self.station_on(jd, k, def, c, i, aliases)?,
            };
            match rim.nodes.iter().position(|&r| r == n) {
                Some(p) => Some(p),
                None => {
                    list.push(n);
                    Some(list.len() - 1)
                }
            }
        };
        let Some(ix) = index_of(x) else {
            return Err(RunFail::Other);
        };
        let iy = if x == y {
            ix
        } else {
            match index_of(y) {
                Some(i) => i,
                None => return Err(RunFail::Other),
            }
        };
        let Ok((order, _)) = crate::arrangement::circular_order(jd, k, def, &list) else {
            return Err(RunFail::Other);
        };
        let n = order.len();
        let pos = |li: usize| order.iter().position(|&o| o == li);
        let (Some(px), Some(py)) = (pos(ix), pos(iy)) else {
            return Err(RunFail::Other);
        };
        let is_rim = |p: usize| order[p] < m;
        // The run's ends: the nearest rim node at or before `x` (clockwise), and at or after `y`.
        let (mut a, mut b) = (px, py);
        while !is_rim(a) {
            a = (a + n - 1) % n;
        }
        while !is_rim(b) {
            b = (b + 1) % n;
        }
        // ★★★★★ **A sector may span several arcs, and then the answer is all of them.** A rim node
        // strictly inside the sector is a θ the *rim* knows and the chart does not — the chart's
        // vertical lines come from `ruling_sweep`, which states a piece only where the lateral
        // face is, so a wall's ruling on the side the face does not reach is absent and the cell
        // is built across it. That is a defect of the **cell**, not of this lookup, and repairing
        // the chart is its own rung; what this can do meanwhile is refuse to pretend the sector
        // has one arc. It hands back the whole run — the rim nodes from `a` to `b`, walked CCW —
        // and the caller answers only when they **agree**, which is exactly what the whole-circle
        // arm of [`Chart::read_cell`] already does with a cut circle's arcs.
        let mut run: Vec<&'a ArcLabel> = Vec::new();
        let mut missing = 0usize;
        let mut cur = a;
        while cur != b {
            let nxt = {
                let mut q = (cur + 1) % n;
                while !is_rim(q) {
                    q = (q + 1) % n;
                }
                q
            };
            let (na, nb) = (list[order[cur]], list[order[nxt]]);
            match arcs.iter().find(|arc| arc.ends == [na, nb]) {
                Some(arc) => run.push(arc),
                // ★ A piece no contribution covers has no arc (the split keeps no
                // phantom). A run made of nothing but such pieces is a sector the face is not
                // over; a run with some of each is a producer inconsistency, as before.
                None => missing += 1,
            }
            cur = nxt;
        }
        if missing > 0 {
            return if run.is_empty() {
                Err(RunFail::AllMissing)
            } else {
                Err(RunFail::Other)
            };
        }
        // `a == b` means the sector's ends land on one rim node: no arc separates them, and the
        // caller has nothing to read here.
        if run.is_empty() {
            return Err(RunFail::Other);
        }
        Ok(run)
    }

    /// **Read one cell off the horizontal lines** — what the chart's emitter reads every cell
    /// with.
    ///
    /// No sign is derived here: `band_is_above` is `planes::plus_t_is_above` at the low end and
    /// its negation at the high end, and the bits come out through [`read_bits`]. The ruling
    /// labels (the
    /// vertical lines) are not consulted: a face that is here has a circle at both its ends, so
    /// the horizontal lines always speak, and `census` asserts that (`src0_present`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn read_cell<'a>(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        kind: BoolKind,
        side: SolidSide,
        cell: &Cell,
        curved: &'a Curved,
        lines: &Lines,
        rows: &[crate::bands::CylRow],
    ) -> Result<CellRead<'a>, BoolError> {
        let t = [
            self.z_lines[cell.interval].t,
            self.z_lines[cell.interval + 1].t,
        ];
        let mut ends: Vec<End<'a>> = Vec::with_capacity(2);
        let mut above: [bool; 2] = [false; 2];
        for e in 0..2 {
            let mut end = End::NoCircle;
            let mut saw_disk = false;
            let mut saw_arc = false;
            for &c in lines.classes(t[e]) {
                // The band leaves the low line toward `+t` and arrives at the high line from
                // `−t`.
                // ★ Bound to the class whose label is read, not to whichever class the loop
                // visited last: one `t` is one plane and so one class today, but the sign and
                // the label must come from the same row the day that stops being true.
                let world = jd.planes[c]
                    .world
                    .as_ref()
                    .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
                let up = crate::planes::plus_t_is_above(world, def);
                let band_is_above = if e == 0 { up } else { !up };
                if let Some(l) = curved.disk_labels.get(&(k, c)) {
                    saw_disk = true;
                    above[e] = band_is_above;
                    end = End::Disk(*l);
                    continue;
                }
                let (Some(arcs), Some(rim)) =
                    (curved.arc_labels.get(&(k, c)), curved.cut_rims.get(&(k, c)))
                else {
                    continue;
                };
                saw_arc = true;
                above[e] = band_is_above;
                end = match cell.walls {
                    None => {
                        // The whole circle is one cell here; every arc must read the same for
                        // *this* interval's side (the other half of a label is the neighbour's).
                        let mut bits: Option<(bool, bool)> = None;
                        let mut same = true;
                        for a in arcs {
                            let b = read_bits(&a.label, side, above[e]);
                            same &= bits.replace(b).is_none_or(|p| p == b);
                        }
                        match (same, arcs.first()) {
                            (true, Some(a)) => End::Disk(a.label),
                            (false, _) | (true, None) => End::Other,
                        }
                    }
                    Some([x, y]) => {
                        let nodes = (self.node_on(x, t[e]), self.node_on(y, t[e]));
                        let m = rim.nodes.len();
                        let adjacent = |nx, ny| {
                            (0..m).any(|j| rim.nodes[j] == nx && rim.nodes[(j + 1) % m] == ny)
                        };
                        match nodes {
                            (Some(nx), Some(ny)) if adjacent(nx, ny) => {
                                // ★ Adjacent rim nodes with no arc between them: the piece no face
                                // covers — the face ends before this sector.
                                match arcs.iter().find(|a| a.ends == [nx, ny]) {
                                    Some(arc) => End::Exact(vec![arc]),
                                    None => End::Uncovered,
                                }
                            }
                            _ => match self.arc_around(
                                jd,
                                k,
                                def,
                                c,
                                t[e],
                                x,
                                y,
                                rim,
                                arcs,
                                &curved.aliases,
                            ) {
                                Ok(run) => End::Exact(run),
                                Err(RunFail::AllMissing) => End::Uncovered,
                                Err(RunFail::Other) => End::Other,
                            },
                        }
                    }
                };
            }
            // A whole-disk label and arcs of one cylinder on one line is a producer inconsistency
            // (`ArcBoundNotYet`'s own sentence), not a chart shape.
            if saw_disk && saw_arc {
                return Err(reject(RejectReason::ArcBoundNotYet));
            }
            ends.push(end);
        }
        let ends: [End<'a>; 2] = match <[End<'a>; 2]>::try_from(ends) {
            Ok(a) => a,
            Err(_) => unreachable!("two ends are pushed, one per line"),
        };

        // ── Membership: every end that speaks must say the same. ──
        let mut chamber: Option<(bool, bool)> = None;
        let mut src2_disagree = false;
        for e in 0..2 {
            let Some(bits) = ends[e].chamber(side, above[e]) else {
                continue;
            };
            match chamber {
                None => chamber = Some(bits),
                Some(prev) if prev != bits => src2_disagree = true,
                Some(_) => {}
            }
        }
        if src2_disagree {
            chamber = None;
        }

        // ── The vertical answer: the chart's other axis, read the same way. ──
        //
        // ★★★★★ **A `Label` carries the material on both sides of its own plane**, so a
        // **ruling**'s label answers the cells beside it exactly as a rim's answers the cells
        // above and below — the chart's two axes are symmetric and only one of them was being
        // read. A cut end the reader cannot pair with its rim (`End::Other`) leaves a cell with
        // no horizontal answer at all, and that is the population this opens.
        //
        // ★★ **It fills in; it does not yet overrule.** Where the horizontal lines speak they
        // stay the answer, and the disagreements are *counted* instead (`probe::shadow`):
        // measured over the boss corpus, the two roads agree on 4,036 cells and disagree on 218
        // — **all of them one family** (`corner Common`, whose axis lies on *two* walls), where
        // both horizontal ends are arcs that agree with each other and the ruling dissents, in
        // the `own` bit only. That family already carries a recorded arrangement label defect
        // (the `(0, 0)` corner), so making the disagreement refuse would regress a
        // defect that is already named rather than fix it. Cross-checking is what this becomes
        // the day that corner is fixed.
        {
            let vertical = |i: usize, starts_here: bool| -> Option<(bool, bool)> {
                let seg = self.theta.get(i)?;
                // A tangent station has no chamber behind it to read: the face ends
                // there, and its `label` is `None` by design — said by the side, not inferred.
                if seg.side == 0 {
                    return None;
                }
                // ★ A station two walls state answers each side from the wall that bounds the
                // chamber there (`StationTwin::bounds_plus`); `starts_here` is the `+θ` side.
                // ☑ Measured: over the edge-on-a-ruling families this never decides a cell — the
                // horizontal lines speak there, and swapping the two walls leaves every boolean
                // as it was — so the rule stands on its derivation here, and on `slit_at`, which
                // reads the same `bounds_plus` and does go red when it is swapped.
                let (l, (wall, sd)) = match seg.twin {
                    Some(tw) if tw.bounds_plus == starts_here => (
                        tw.label?,
                        (tw.wall, ruling_side_of(jd, k, def, tw.wall, seg.end[0])?),
                    ),
                    _ => (seg.label?, self.ruling_name(jd, k, def, i)?),
                };
                // The sector leaves `x` counter-clockwise and arrives at `y`, so the two walls
                // are read from opposite sides of their own rulings.
                let above = crate::arrangement::plus_theta_is_above(jd, wall, sd)? == starts_here;
                Some(read_bits(&l, side, above))
            };
            // ★ A cell whose two walls are **one** ruling (a circle opened at a single point) lies
            // on both sides of that wall, so the vertical line says nothing about it — the two
            // readings would contradict by construction.
            let vert: Vec<(bool, bool)> = match cell.walls {
                Some([x, y]) if x != y => [vertical(x, true), vertical(y, false)]
                    .into_iter()
                    .flatten()
                    .collect(),
                _ => Vec::new(),
            };
            #[cfg(test)]
            probe::shadow::record(chamber, &vert, cell.walls.is_some_and(|[x, y]| x == y));
            // The two walls must agree with each other before either may answer.
            if chamber.is_none() && !src2_disagree {
                chamber = match vert.as_slice() {
                    [a] => Some(*a),
                    [a, b] if a == b => Some(*a),
                    _ => None,
                };
            }
        }

        // ── Existence: the trace where an end is cut, the span where none is. ──
        let (lo, hi) = (t[0].min(t[1]), t[0].max(t[1]));
        let by_span = rows.iter().any(|r| {
            r.class == k && r.span[0].min(r.span[1]) <= lo && hi <= r.span[0].max(r.span[1])
        });
        let mut by_marks: Option<bool> = None;
        for e in 0..2 {
            let End::Exact(arcs) = &ends[e] else { continue };
            // ★ **Existence asks the same run membership does.** A sector spanning several arcs
            // is present only if every one of them says so; arcs that disagree mean the cell
            // straddles a boundary the chart has no line for, and that is the refusal below.
            //
            // ☑ This is the half that moves results, and what it answers is **absent**: every run
            // the crossing corpus builds is over arcs whose `marks` are empty, which is the
            // geometry (the quadrant the chart could not name lies on the side the panel does not
            // reach). Forcing the run's answer to `true` here trips the record-site assertion
            // below with such an arc as its witness; forcing it to `false` changes nothing, which
            // is the same fact from the other side. An axial span knows no θ, and so would claim
            // every sector present.
            let mut span: Option<bool> = None;
            for arc in arcs {
                // `face_spans` refuses by name (two of this solid's faces disagreeing on one arc), and
                // two cut ends disagreeing with each other is `CylinderFaceUndecided`: the rims of
                // a hole are band boundaries, so a sector exists over its whole height or
                // not at all, and picking an end to believe is the guess this kernel does not make.
                let v = face_spans(arc, side, above[e])?;
                if span.replace(v).is_some_and(|p| p != v) {
                    return Err(reject(RejectReason::CylinderFaceUndecided));
                }
            }
            let Some(v) = span else { continue };
            if by_marks.replace(v).is_some_and(|p| p != v) {
                return Err(reject(RejectReason::CylinderFaceUndecided));
            }
        }
        // ★★★★★ **Two silent ends and a span that says «present» is not read as present**.
        // The span knows no θ, so over a cell whose both cut ends could not be paired with
        // their rims it would claim the face is there — the guess that puts a face over a hole
        // in the corner boss × mid slab. With stations placed by name (`station_on`) no cell in
        // the suite reaches
        // here (`src0_present` 0, the census's record-site assertion made structural), so
        // this is a guard on the reader's premise and the name is the chart's own for a cell it
        // cannot read.
        if by_marks.is_none() && by_span && ends.iter().all(|e| matches!(e, End::Other)) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
        // ★ The span is the coarser truth: a holed lateral's row spans the hole, and the trace is
        // what says the face is *not* there (the `exist_marks_false` population). So the
        // disagreement that would be a defect is only the other direction: the trace claiming a
        // face where no row spans.
        // ★ An uncovered rim piece decides existence outright: the face's own rim does
        // not run along this sector at that line, so the face is not over this cell — whatever
        // the row's axial span says.
        let uncovered = ends.iter().any(|e| matches!(e, End::Uncovered));
        let (present, exist_disagree) = match by_marks {
            _ if uncovered => (false, false),
            Some(v) => (v, v && !by_span),
            None => (by_span, false),
        };

        // A cell the face is not in emits nothing, chamber or no chamber — a cell outside every
        // span has no circle at either end and no question to answer.
        let emit = if !present {
            Some(false)
        } else {
            chamber.map(|(own, other)| {
                let keep = |in_own: bool| keep_for(kind, side, in_own, other);
                keep(own) != keep(!own)
            })
        };
        Ok(CellRead {
            ends,
            chamber,
            src2_disagree,
            present,
            exist_disagree,
            emit,
        })
    }
}

/// **Does this row's lateral face reach into the interval, at this arc?** — the *existence*
/// question, which no label answers.
///
/// A label states where material is. It is written about a **cell**, and it stays the same whether
/// the cylinder's own face bounds that cell or some other face does. That is enough while a
/// lateral marks a class over its whole circle; it stops being enough the moment a **holed**
/// lateral comes back in as an operand, because then two sectors of one rim can carry a literally
/// identical label and differ only in whether the face is there at all.
///
/// The trace already said which. One rule, no shapes counted:
///
/// | this face's mark on the arc | reaches into the interval |
/// |---|---|
/// | [`SegKind::Transversal`] | **yes** — the face passes through the plane here |
/// | [`SegKind::Graze`] | only when `body_above` names the interval's side |
/// | nothing | **no** — the face stops short of this arc |
///
/// `band_is_above` is the interval's side of this rim's plane in that plane's **stored** frame —
/// the very argument [`read_bits`] takes, as the caller already holds it. There is no second
/// derivation of "which way is the band" here, deliberately: a sign re-derived one call away from
/// its twin is where the two drift apart.
///
/// ★ [`SegKind::Seated`] is skipped because it is a *planar* face's word — see [`ArcLabel::marks`],
/// where the type is what rules a lateral out, not a convention.
///
/// ★★ **The rule is not about holes.** An outer rim grazes too (`circle_on_class` says
/// `Grazes { body_above: up }` at the face's own end), and there the interval on the face's side
/// gets `true` from this same test — one sentence about every rim. What is out of reach is only
/// an **uncut** outer rim: it never becomes an [`ArcLabels`] entry at all (a whole circle goes to
/// `DiskLabels` and the band road), so today only cut rims ask here.
///
/// ☑ **How often each row actually fires** (whole binary, production calls only):
/// `Transversal` **265** · `Graze` **13**, of which **1** reaches and **12**
/// do not · nothing at all **0**. The chart reads it once per cut end (`cyl_chart::census`'s
/// `arcs_read`, 756 over the lib suite).
/// The twelve are the six dropped sectors read at both rims, which is the cross-check that the
/// sector census and this one describe the same events.
pub(crate) fn face_spans(
    r: &ArcLabel,
    side: SolidSide,
    band_is_above: bool,
) -> Result<bool, BoolError> {
    let mut answer: Option<bool> = None;
    for (_, kind) in r.marks.iter().filter(|(s, _)| *s == side) {
        let reaches = match *kind {
            // A planar face's word — see `ArcLabel::marks`; a tangent ruling is a seated face's
            // too, and never a rim arc's mark.
            SegKind::Seated { .. } | SegKind::Tangent { .. } => continue,
            SegKind::Transversal { .. } => true,
            SegKind::Graze { body_above } => body_above == band_is_above,
        };
        // ★★★ **Two of this solid's faces meeting at one arc and *disagreeing*: refused, and the
        // refusal is a placeholder.** It takes one cylinder class carrying several faces that
        // share a rim — reachable geometry (a split bore whose two bands touch), just not
        // reachable today. The answer that day is almost certainly `.any()`: if any of the
        // solid's faces reaches into the interval, a face is there. It is not written that way
        // because it cannot be measured, and an unmeasured rule is a guess. Marks that **agree**
        // decide nothing by themselves, so they are
        // taken: refusing there would be refusing a case with no guess in it.
        if answer.replace(reaches).is_some_and(|prev| prev != reaches) {
            return Err(reject(RejectReason::CylinderFaceUndecided));
        }
    }
    Ok(answer.unwrap_or(false))
}

/// The one spelling of "which two bits of a label are this cell's chamber": the row's own
/// solid's bit and the counterpart's, on the side of the plane the cell occupies. Read for a disk
/// end and an arc end alike (`Chart::read_cell`), so the two cannot drift.
///
/// ★★ **The chamber is read off the arrangement, not measured.** The plane arrangement already
/// makes a **disk cell** for every circle a cylinder leaves on a class, and `label_cells` writes
/// four bits on it: which solid's material lies immediately above and below that plane *inside the
/// circle*. That is the band's chamber, stated by the engine that decided it.
///
/// ★ This replaced a witness ray (`point_in_faces_rat`): a rational point on the axis, cast
/// through the counterpart's faces, counting crossings. It gave the same answers, but it answered
/// a **3D containment** question — the shape the *component* probe asks — when the band's question
/// is the same shape as "does this face survive": two chambers either side of a boundary. Reading
/// the label needs no coordinates, no ray direction, no abstention retry, and has no width
/// ceiling, and it is why a cylinder may now stand on either side of the boolean.
pub(crate) fn read_bits(l: &Label, side: SolidSide, band_is_above: bool) -> (bool, bool) {
    let (cyl_bit, other_bit) = match side {
        SolidSide::A => (0usize, 2usize), // [A above, A below, B above, B below]
        SolidSide::B => (2usize, 0usize),
    };
    let i = usize::from(!band_is_above);
    (l[cyl_bit + i], l[other_bit + i])
}

/// The one spelling of the band's keep decision — the wall is a boundary face of its own solid,
/// so that solid's membership flips across it while the counterpart's does not.
pub(crate) fn keep_for(kind: BoolKind, side: SolidSide, in_own: bool, in_other: bool) -> bool {
    match side {
        SolidSide::A => crate::draft::keep(kind, in_own, in_other),
        SolidSide::B => crate::draft::keep(kind, in_other, in_own),
    }
}
