use super::*;
/// **What one line says to one cell** — the line's answer, as the cell reads it.
///
/// The plane arrangement wrote a label on every section a cylinder leaves on a class it traced,
/// whether or not a face was kept there (`emit_faces`' `disk_labels`/`arc_labels`, collected outside the keep
/// filter). So a cell's chamber is not *decided* here; it is **read** off the line at either end,
/// a disk or an arc ([`Chart::read_cell`] is the one reader).
pub(crate) enum End<'a> {
    /// The circle on this line is uncut: one label for the whole disk, whatever the slab.
    Disk(&'a DiskLabel),
    /// The circle is cut, and the cell's slab lies under one of its arcs — one, because every
    /// node the rim is cut at is an event of the chart ([`Chart::events`]), so no slab runs past
    /// one.
    Arc(&'a ArcLabel),
    /// The circle is cut, and no arc covers the cell's slab: the piece of the circle no face of
    /// this class has as a rim (the three quarters a fillet's quarter leaves, the half a slot's end
    /// leaves; the split drops such a piece rather than keep a phantom edge). The lateral face is
    /// **not here**: this end decides existence, not a chamber.
    Uncovered,
    /// No label of this cylinder on this line: a cut rim with no arc kept. The chart's lines are
    /// the classes that traced this cylinder ([`chart_of`]), so a present cell reading this is
    /// the stages disagreeing.
    NoCircle,
}

/// **What several readings of one chamber say together** — the one spelling of «unanimous or
/// nothing». Between a cell's sides a [`Agreement::Split`] is two statements of one chamber that
/// contradict each other.
enum Agreement<T> {
    /// No reading at all.
    Silent,
    /// Every reading the same.
    One(T),
    /// Two readings differ.
    Split,
}

fn agreement<T: PartialEq>(readings: impl IntoIterator<Item = T>) -> Agreement<T> {
    let mut out = Agreement::Silent;
    for r in readings {
        out = match out {
            Agreement::Silent => Agreement::One(r),
            Agreement::One(p) if p == r => Agreement::One(p),
            Agreement::One(_) | Agreement::Split => return Agreement::Split,
        };
    }
    out
}

impl End<'_> {
    /// **What this end says about the cell's chamber**, or `None` when it says nothing. The read
    /// is [`read_bits`], the one spelling.
    fn chamber(&self, side: SolidSide, above: bool) -> Option<(bool, bool)> {
        match self {
            End::Disk(d) => Some(read_bits(&d.label, side, above)),
            End::Arc(a) => Some(read_bits(&a.label, side, above)),
            End::Uncovered | End::NoCircle => None,
        }
    }
}

/// One cell, read.
/// `ends` is the census's reader; production reads `chamber`,
/// `present`, `emit` and `unread` (the emitter) and pays for the rest only as a copy.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct CellRead<'a> {
    pub(crate) ends: [End<'a>; 2],
    /// `(in_own, in_other)` on the cylinder's inside, agreed by every side that could speak — the
    /// two ends' rims and the two walls' rulings. `None` when no side could, or when two sides
    /// disagreed.
    pub(crate) chamber: Option<(bool, bool)>,
    /// Why `chamber` is empty, as the refusal a present cell gets: two sides contradicting each
    /// other is labels the arrangement wrote disagreeing (`LabelConflict`); no side speaking over
    /// a present cell is the stages disagreeing. `None` when there is a chamber.
    pub(crate) unread: Option<RejectReason>,
    /// Is this lateral face here at all — the existence question, answered by the trace's marks on
    /// the lines at the cell's ends ([`face_spans`]).
    pub(crate) present: bool,
    /// Would the cell be a result face: present, with a chamber, and `keep` differing across the
    /// wall. `None` when the chamber is unknown.
    pub(crate) emit: Option<bool>,
}

impl Chart {
    /// **Whether the lateral ends at station `i` on both sides** — where the material the lateral
    /// bounds does not run round the line: a sector between two neighbouring spokes of the line's
    /// star — the secant walls' chords ([`Station::chords`]) and the other solid's tangent wall
    /// where it ends on the line ([`Station::tangent`]) — **on the lateral's material side**, is not
    /// in the result. All kept, the material runs across the line and so does the lateral; one not
    /// kept, what lies on the two sides meets on the line alone. `false` for a tangent ruling
    /// (`side == 0`, whose `station` is `None` — no chord runs through it), for a station one chord
    /// alone runs through with no tangent wall ending there — nothing lies between the lateral's
    /// halves but that chord — and for one with a wall that
    /// draws no piece there ([`ChordSays::NoPiece`]): the line is then the lateral's own edge, and
    /// the lateral ends at it already (the walk refuses a station where it does not —
    /// `regions::walk`).
    ///
    /// **A tangent wall ending on the line** lies outside the cylinder, along `±τ`: between the
    /// lateral's half on that side, which bends to the axis, and every chord's outer half. Where the
    /// lateral keeps the outside, it splits the sector beside the lateral into the lateral's own
    /// (kept, as below) and the one next to the nearest outer half, which that chord's outside cell
    /// states — the last chord's on its `τ` side beside `+τ`, the first chord's on the other side
    /// beside `−τ`. A notched block less a box whose corner stands on the notch's lateral, one side
    /// tangent there and one across it, is two bodies meeting on that line only where this sector
    /// is read (`edge_on_a_ruling`).
    ///
    /// **Each chamber between is read off both its walls**: chord `k`'s inside cell on its far side
    /// from `τ` and chord `k + 1`'s on its near side, and outside, the vertically opposite chamber
    /// between their outer halves, off their outside cells the other way round. Two witnesses that
    /// disagree are labels the arrangement wrote disagreeing (`LabelConflict`). The material side
    /// is the one the first chord's inside cell answers beside `+θ`: kept, the lateral bounds
    /// material inside; not, outside. That reading only has to be true where the answer is used —
    /// the walk cuts a lateral at a slit only between two cells it emits, and an emitted `+θ` cell
    /// keeps one side of the lateral and not the other. One boolean keeps a lateral's material on
    /// one side all round a line: `keep` is monotone in each operand, and an operand, being a
    /// manifold, has its own material on one side of the cylinder on both sides of the line
    /// (measured: no station with its two cells emitted keeps different sides, in the suite or the
    /// census).
    ///
    /// Both halves carry an answer (`edge_on_a_ruling`): the inside chamber read where the
    /// material is outside cuts the lateral of a box less the wedge fused on a rim, across the
    /// pinch — the chamber between the chords is the wedge's own — and no chamber read there
    /// builds a box with a V-notch whose apex lies on the lateral, less the cylinder, as a body
    /// touching itself along the line.
    ///
    /// **What a station cannot answer**: a result face on the line from the lateral's *other*
    /// side — the walls' halves where the lateral keeps nothing (the inside chamber kept where the
    /// material is outside, the outside chamber kept where it is inside). The line is then an
    /// edge of those faces while the lateral runs on across it, which is right when they belong
    /// to another body touching along the line and wrong when the body is the lateral's own,
    /// joined elsewhere. Only the shells can tell, and cutting the lateral here breaks the first;
    /// the result's self-touch verdict asks it once the solids are known — the walls' edge inside
    /// the solid's own lateral (`assembly::self_touch`).
    ///
    /// `Err` when a label or a side could not be read: the labels are formed for every secant
    /// ruling piece, so their absence is the stages disagreeing.
    pub(crate) fn slit_at(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        i: usize,
        side: SolidSide,
        kind: BoolKind,
    ) -> Result<bool, RejectReason> {
        let stages = RejectReason::CylinderStagesDisagree;
        let Some(st) = &self.theta[i].station else {
            return Ok(false);
        };
        if st
            .chords
            .iter()
            .any(|c| matches!(c.says, ChordSays::NoPiece))
        {
            return Ok(false);
        }
        // Each wall's inside and outside labels and which side of it `τ` is on.
        let mut walls: Vec<(Label, Option<Label>, bool)> = Vec::with_capacity(st.chords.len());
        for c in &st.chords {
            let ChordSays::Piece { label, outer } = c.says else {
                unreachable!("a station with a wall that draws no piece returned above");
            };
            let near = crate::arrangement::plus_theta_is_above(jd, c.wall, c.side).ok_or(stages)?;
            walls.push((label.ok_or(stages)?, outer, near));
        }
        let keeps = |b: (bool, bool)| keep_for(kind, side, b.0, b.1);
        let inside = keeps(read_bits(&walls[0].0, side, walls[0].2));
        let mut between: Vec<bool> = Vec::with_capacity(walls.len());
        for pair in walls.windows(2) {
            let ((l0, o0, n0), (l1, o1, n1)) = (pair[0], pair[1]);
            let (a, b) = if inside {
                (read_bits(&l0, side, !n0), read_bits(&l1, side, n1))
            } else {
                (
                    read_bits(&o0.ok_or(stages)?, side, n0),
                    read_bits(&o1.ok_or(stages)?, side, !n1),
                )
            };
            if a != b {
                return Err(RejectReason::LabelConflict);
            }
            between.push(keeps(a));
        }
        // The lateral's two halves and the chords between them on its material side, in chord
        // order. The chambers next to the lateral are kept and the far side is not — the premise
        // above: the answer is read only between two cells the walk emits, each keeping one side
        // of the lateral.
        let lateral = 0;
        let n = walls.len();
        // Each spoke with the sector counter-clockwise after it.
        let mut star: Vec<(usize, bool)> = vec![(lateral, true)];
        if inside {
            // Inside, nothing but chords: the tangent plane lies outside the cylinder. The sector
            // after the last chord is the lateral's own.
            star.extend((1..=n).zip(between.into_iter().chain([true])));
        } else {
            // Outside, counter-clockwise from the lateral's `+τ` half: the tangent wall's `+τ` half
            // where it ends on the line, then the chords' outer halves from the last to the first,
            // then the `−τ` half. Beyond a tangent half lies the sector next to the nearest outer
            // half — the last chord's on its `τ` side, the first chord's on the other — which that
            // chord's outside cell states; the sector between the lateral and the tangent half is
            // the lateral's own, kept by the premise above.
            let tangent = |at: usize| n + 1 + at;
            let outside = |w: usize, toward_tau: bool| -> Result<bool, RejectReason> {
                let (_, o, near) = walls[w];
                Ok(keeps(read_bits(
                    &o.ok_or(stages)?,
                    side,
                    near == toward_tau,
                )))
            };
            if st.tangent[0] {
                star.push((tangent(0), outside(n - 1, true)?));
            }
            let beyond_first = if st.tangent[1] {
                outside(0, false)?
            } else {
                true
            };
            star.extend(
                (1..=n)
                    .rev()
                    .zip(between.into_iter().rev().chain([beyond_first])),
            );
            if st.tangent[1] {
                star.push((tangent(1), true));
            }
        }
        star.push((lateral, false));
        let star = crate::draft::LineStar::new(star);
        Ok(star.continues(lateral) == Some(false))
    }

    /// **The arc of class `c`'s cut rim over slab `u`** — the arc from the rim's last node at or
    /// before the slab's first event to the next node round, counter-clockwise: every rim node is
    /// an event, so that is the one arc the slab lies under. A rim cut at **one** node has one arc,
    /// from that node back to itself (a polygon vertex cut the circle there and nothing else did —
    /// `combinatorics::ring_floor`), and it covers every slab. `None` when the split kept no arc
    /// there — no contribution covers the piece.
    fn arc_over<'a>(
        &self,
        c: usize,
        u: usize,
        rim: &crate::draft::CutRim,
        arcs: &'a [ArcLabel],
    ) -> Result<Option<&'a ArcLabel>, BoolError> {
        let evs = self
            .rim_events
            .get(&c)
            .ok_or_else(|| reject(RejectReason::CylinderStagesDisagree))?;
        let m = rim.nodes.len();
        if m == 0 || evs.len() != m {
            return Err(reject(RejectReason::CylinderStagesDisagree));
        }
        let p = evs.iter().rposition(|&e| e <= u).unwrap_or(m - 1);
        let ends = [rim.nodes[p], rim.nodes[(p + 1) % m]];
        Ok(arcs.iter().find(|a| a.ends == ends))
    }

    /// **Read one cell off the lines around it** — what the chart's emitter reads every cell
    /// with.
    ///
    /// No sign is derived here: `band_is_above` is `planes::plus_t_is_above` at the low end and
    /// its negation at the high end, a wall's side is `plus_theta_is_above`, and the bits come out
    /// through [`read_bits`]. A face that is here has a circle at both its ends, so the horizontal
    /// lines always speak (`census` asserts `src0_present == 0`); the walls' rulings check them.
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
    ) -> Result<CellRead<'a>, BoolError> {
        let mut ends: Vec<End<'a>> = Vec::with_capacity(2);
        let mut above: [bool; 2] = [false; 2];
        for (e, &line) in cell.lines.iter().enumerate() {
            let c = self.lines[line];
            // The band leaves the low line toward `+t` and arrives at the high line from `−t`.
            let world = jd.planes[c]
                .world
                .as_ref()
                .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
            let up = crate::planes::plus_t_is_above(world, def);
            above[e] = if e == 0 { up } else { !up };
            let disk = curved.disk_labels.get(&(k, c));
            let cut = curved.split_rims.get(&(k, c));
            let end = match (disk, cut) {
                // A whole-disk label and a cut of one circle is the split and the labels disagreeing
                // about whether that circle was cut, not a chart shape.
                (Some(_), Some(_)) => return Err(reject(RejectReason::CylinderStagesDisagree)),
                (Some(d), None) => End::Disk(d),
                (None, Some(rim)) => match curved.arc_labels.get(&(k, c)) {
                    Some(arcs) => match self.arc_over(c, cell.slab, rim, arcs)? {
                        Some(arc) => End::Arc(arc),
                        None => End::Uncovered,
                    },
                    None => End::NoCircle,
                },
                (None, None) => End::NoCircle,
            };
            ends.push(end);
        }
        let ends: [End<'a>; 2] = match <[End<'a>; 2]>::try_from(ends) {
            Ok(a) => a,
            Err(_) => unreachable!("two ends are pushed, one per line"),
        };

        // ── Membership: every side that speaks must say the same. ──
        //
        // **A `Label` carries the material on both sides of its own plane**, so a
        // **ruling**'s label answers the cells beside it exactly as a rim's answers the cells
        // above and below — the chart's two axes state one chamber, and the rule is one for all
        // four sides: the chamber stands when every side that speaks says the same, and two sides
        // that contradict each other leave none (the emitter then refuses a present cell by name).
        // A contradiction is an arrangement label defect, and reading one side over the other would
        // hide it. Measured over the suite, the ignored sweep and the census:
        // on present cells both axes speak on 21,006 / 27,300 / 1,026 and the horizontal alone on
        // 4,048 / 3,260 / 978; no side has contradicted another anywhere.
        //
        // The symmetry is the chamber's: existence below is the rims' marks, and a ruling's own
        // marks are an instrument's field, so a wall says nothing there.
        let vertical = |i: usize, starts_here: bool| -> Option<(bool, bool)> {
            let seg = self.theta.get(i)?;
            // A tangent station has no chamber behind it to read: the face ends
            // there, and its `label` is `None` by design — said by the side, not inferred.
            if seg.side == 0 {
                return None;
            }
            // A ruling's station answers each side from the wall that closes the
            // chamber there — the first chord beside `+θ`, the last beside `−θ`
            // (`Station::chords`); `starts_here` is the `+θ` side. Held here as well as in
            // `slit_at`, which reads the same order (and the walls' outside cells too): swapping
            // the ends turns the corner and inward-wedge families of `edge_on_a_ruling` red, where
            // the rims speak too and the swapped wall contradicts them. A closing wall with no
            // piece there has no label to read, and the station's own wall would answer for a
            // chamber between chords — so that side is silent.
            let st = seg.station.as_ref()?;
            let c = if starts_here {
                st.chords.first()?
            } else {
                st.chords.last()?
            };
            let (l, (wall, sd)) = match c.says {
                ChordSays::Piece { label, .. } => (label?, (c.wall, c.side)),
                ChordSays::NoPiece => return None,
            };
            // The slab leaves its first event counter-clockwise and arrives at its second, so the
            // two walls are read from opposite sides of their own rulings.
            let above = crate::arrangement::plus_theta_is_above(jd, wall, sd)? == starts_here;
            Some(read_bits(&l, side, above))
        };
        // A cell whose two walls are **one** ruling (a circle opened at a single point) lies
        // on both sides of that wall, so the vertical line says nothing about it — the two
        // readings would contradict by construction. An event with no ruling across the gap
        // has no wall to read, and a wall that cannot be read (a tangent station, a side it cannot
        // name) is silent: it checks the chamber, it does not decide it — a present cell always
        // has a horizontal answer (the census asserts `src0_present == 0`).
        let walls = match cell.walls {
            [Some(x), Some(y)] if x == y => [None, None],
            [x, y] => [
                x.and_then(|x| vertical(x, true)),
                y.and_then(|y| vertical(y, false)),
            ],
        };
        let (chamber, disagree) = match agreement(
            (0..2)
                .filter_map(|e| ends[e].chamber(side, above[e]))
                .chain(walls.into_iter().flatten()),
        ) {
            Agreement::One(b) => (Some(b), false),
            Agreement::Silent => (None, false),
            Agreement::Split => (None, true),
        };
        let unread = match (chamber, disagree) {
            (Some(_), _) => None,
            (None, true) => Some(RejectReason::LabelConflict),
            // No side spoke over a present cell: both its lines carry no circle of this cylinder,
            // against «a face that is here has a circle at both its ends».
            (None, false) => Some(RejectReason::CylinderStagesDisagree),
        };

        // ── Existence: what the trace left on the lines at the cell's two ends. ──
        // A line carries this solid's lateral marks over the cell's slab — an arc's when its circle
        // is cut, the whole circle's when it is not — and `face_spans` reads them: through the
        // plane, on the side a rim bounds, or not here at all. `face_spans` refuses by name (two
        // of this solid's faces disagreeing on one line), and two ends disagreeing with each other
        // is `CylinderFaceUndecided`: a slab exists over its whole height or not at all, and
        // picking an end to believe is the guess this kernel does not make.
        let mut by_marks: Option<bool> = None;
        for e in 0..2 {
            let marks: &[(SolidSide, SegKind)] = match &ends[e] {
                End::Arc(arc) => &arc.marks,
                End::Disk(disk) => &disk.marks,
                End::Uncovered | End::NoCircle => continue,
            };
            let v = face_spans(marks, side, above[e])?;
            if by_marks.replace(v).is_some_and(|p| p != v) {
                return Err(reject(RejectReason::CylinderFaceUndecided));
            }
        }
        // An uncovered rim piece decides existence outright: the face's own rim does not run along
        // this slab at that line, so the face is not over this cell. A cell neither end speaks for
        // — both lines a cut rim with no arc kept — has no trace to read, and a face is never read
        // off nothing: the stages disagree.
        let uncovered = ends.iter().any(|e| matches!(e, End::Uncovered));
        let present = match by_marks {
            _ if uncovered => false,
            Some(v) => v,
            None => return Err(reject(RejectReason::CylinderStagesDisagree)),
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
            unread,
            present,
            emit,
        })
    }
}

/// **Does this solid's lateral face reach into the interval, at this end line?** — the *existence*
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
/// [`SegKind::Seated`] is skipped because it is a *planar* face's word — see [`ArcLabel::marks`],
/// where the type is what rules a lateral out, not a convention.
///
/// **The rule is not about holes.** An outer rim grazes too (`circle_on_class` says
/// `Grazes { body_above: up }` at the face's own end), and there the interval on the face's side
/// gets `true` from this same test — one sentence about every rim. A circle no split cut asks with
/// its whole marks ([`DiskLabel`]), an arc with its own.
///
/// **How often each row fires on cut ends** (the lib suite, an arc's marks read): `Transversal`
/// **23,744** · `Graze` **9,878**, of which **9,267** reach and **611** do not · nothing at all
/// **3**. The chart reads it once per end that carries marks.
pub(crate) fn face_spans(
    marks: &[(SolidSide, SegKind)],
    side: SolidSide,
    band_is_above: bool,
) -> Result<bool, BoolError> {
    let mut answer: Option<bool> = None;
    for (_, kind) in marks.iter().filter(|(s, _)| *s == side) {
        let reaches = match *kind {
            // A planar face's word — see `ArcLabel::marks`; a tangent ruling is a seated face's
            // too, and never a rim arc's mark.
            SegKind::Seated { .. } | SegKind::Tangent { .. } => continue,
            SegKind::Transversal { .. } => true,
            SegKind::Graze { body_above } => body_above == band_is_above,
        };
        // **Two of this solid's faces meeting at one arc and *disagreeing*: refused, and the
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
/// **The chamber is read off the arrangement, not measured.** The plane arrangement already
/// makes a **disk cell** for every circle a cylinder leaves on a class, and `label_cells` writes
/// four bits on it: which solid's material lies immediately above and below that plane *inside the
/// circle*. That is the band's chamber, stated by the engine that decided it.
///
/// A witness ray (`point_in_faces_rat`) — a rational point on the axis, cast through the
/// counterpart's faces, counting crossings — would give the same answers, but it answers a **3D
/// containment** question — the shape the *component* probe asks — when the band's question is the
/// same shape as "does this face survive": two chambers either side of a boundary. Reading the
/// label needs no coordinates, no ray direction, no abstention retry, and has no width ceiling, and
/// it is why a cylinder may stand on either side of the boolean.
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
