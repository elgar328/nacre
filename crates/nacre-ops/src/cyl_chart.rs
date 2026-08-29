//! **A cylinder's own chart, as an arrangement's line set** (capability D, first rung).
//!
//! The plane side long ago stopped walking faces by hand: each plane class gets a **cell complex**
//! (`arrangement`'s `walk_cells → nest_cells → label_cells → emit_faces`) and the case-work went
//! with it. The cylinder side did not. It still has three hand-written walks — `bands.rs`' band and
//! panel roads, and `boolean.rs`' `band_loop` slit-leg and `merge_curved_group` — and the last four
//! cells were each a falsehood inside one of them.
//!
//! `docs/design.md` names the way out: the lateral has an **isometric chart** `(z, r·θ)`, so the
//! same engine can run there — the seam is only the chart's cut line, and notches, holes and
//! θ-panels all become cells. This module is that road's **first** piece: the arrangement's two
//! axes, and nothing else.
//!
//! ## Why the lines are orthogonal — a derivation, not a measurement
//!
//! The M6-2a population gate (`planes.rs`) sorts every `(plane class, cylinder class)` pair into
//! five, and only two leave a mark on the lateral:
//!
//! | plane vs axis | distance | on the lateral | ☑ suite |
//! |---|---|---|---|
//! | `n ∥ m` | — | a **circle** — the chart's *horizontal* | 1128 |
//! | `n ⊥ m` | `d = 0` | **two rulings** — the chart's *vertical* | 75 |
//! | `n ⊥ m` | `0 < d ≤ r` | refused, `WallMeetsLateral` (capability B) | 7 |
//! | `n ⊥ m` | `d > r` | nothing | 1494 |
//! | oblique | — | refused, `ObliqueCylinderCut` (M6-3) | 2 |
//!
//! So within this milestone the chart carries **axis-parallel segments only**. Not a grid, though:
//! a ruling exists over a finite axis interval and a cut circle is arcs, so the cells are what
//! rectilinear *segments* cut out — which is why building them is its own rung.
//!
//! ★ The gate is the reason, so the day the gate changes this does too. Capability B opening
//! `0 < d ≤ r` adds vertical lines at **irrational** θ — and the chart never reads a θ *value*,
//! only an order, so nothing here would break. That is plausible and **unmeasured**; it is not
//! a reason this module leans on.
//!
//! ## What is here, and what is not
//!
//! The two axes (D1a) and the **cells** they cut (D1b), beside a census that counts how today's
//! emitted lateral faces cover those cells. The **labels** and any cutover are not: the two roads
//! this will one day replace are untouched, because the census has to measure what they answer
//! *today*.

use crate::arrangement::{Curved, RulingExtent};
use crate::boolean::{Bound, LocalFace};
use crate::planes::{ClassIx, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, combinatorics};
use nacre_scalar::Rat;

/// One **horizontal** line: a ⊥ plane class, and where it crosses the axis.
///
/// ★ `t` is exact (`Rat`) — a ⊥ class's axis parameter is a rational, and `bands::param` is the
/// one spelling of it. A class whose parameter is too wide for `Rat` refuses rather than being
/// skipped, exactly as `bands_of` already does: a missing band boundary merges two regions whose
/// membership differs, which is a closed and wrong answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ZLine {
    pub(crate) class: usize,
    pub(crate) t: Rat,
}

/// One **vertical** segment: a ruling, and the axis interval it spans.
///
/// Carried straight from [`RulingExtent`] — the plane pass decides this and the chart must not
/// decide it again.
///
/// ★ No `side`: a ruling is named `(wall class, root)` and the side derives from that — see
/// [`RulingExtent`], which dropped its own copy for the same reason.
#[derive(Clone, Debug)]
pub(crate) struct ThetaSeg {
    pub(crate) wall: usize,
    /// The piece's own two branch nodes: the θ **order** is asked of these
    /// (`circular_order_about_seam` via `arrangement::circular_order`), never of a coordinate.
    pub(crate) end: [combinatorics::NodeId; 2],
    /// `[lower, upper]` on the axis.
    pub(crate) z: [Rat; 2],
    /// **The chart's vertical answer** — the label of the cell this ruling borders inside the
    /// cylinder, carried straight from [`RulingExtent::label`]. The horizontal lines' answers are
    /// `DiskLabels`/`ArcLabels`; this is the half they never had.
    pub(crate) label: Option<crate::arrangement::Label>,
    /// Who traced this line — [`RulingExtent::marks`]. Membership is `label`'s question; whether
    /// this lateral face is even here is this one's.
    pub(crate) marks: Vec<(crate::planes::SolidSide, crate::arrangement::SegKind)>,
}

/// One cylinder class's chart: the two axes, and the cells they cut.
#[derive(Clone, Debug, Default)]
pub(crate) struct Chart {
    pub(crate) z_lines: Vec<ZLine>,
    pub(crate) theta: Vec<ThetaSeg>,
}

/// A ruling's **name**: the wall class that cut it, and which of the two roots it is.
///
/// ★ One name, not two. `RulingExtent` also carries a `side`, but that is *derived* from the same
/// pair (`arrangement::ruling_side`), and a second spelling of one identity is how this ladder
/// has been bitten before. The join against `panel_faces`' rings uses this key because that is
/// the key `name_on` already builds there.
type RulingName = (usize, nacre_topo::QuadRoot);

/// **One cell of the chart** — an axis interval crossed with a θ-sector.
///
/// ★★★ Addressed by `(interval, sector)` rather than kept as a flat list, because the next rung
/// asks for **neighbours**: crossing a horizontal line moves to the interval above or below,
/// crossing a ruling to the sector beside. That adjacency is *derivable* from this address, and
/// deriving it is D2's job — building it here would be a field with no consumer.
///
/// ☑ The adjacency is not 1:1 and the address already says so: an interval carrying one cell sits
/// against every sector of a neighbour that carries six.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Cell {
    /// The interval between `z_lines[interval]` and `z_lines[interval + 1]`.
    pub(crate) interval: usize,
    /// Which sector of that interval, in the interval's own θ order.
    pub(crate) sector: usize,
    /// The sector's two bounding rulings **as indices into [`Chart::theta`]**, `None` when the
    /// interval carries none and the cell is the whole circle. For a single cut the circle opens
    /// at one point and both are that ruling.
    ///
    /// ★★★ **Indices, not names.** The first spelling kept the `(wall, root)` names here because
    /// that is the key `panel_faces` joins on — but then the ruling's *label* (D2a) could not be
    /// reached from a cell without searching for it, and a second copy of one identity is what
    /// this ladder keeps being bitten by. The name is one call away ([`Chart::ruling_name`]).
    pub(crate) walls: Option<[usize; 2]>,
}

impl Chart {
    /// **The cells the two axes cut** — `None` if any θ order could not be formed.
    ///
    /// ★★★★★ **Every ruling's axis endpoints are themselves z-lines** — measured across the suite
    /// (263 charts, 0 exceptions) and structural besides: a ruling's node is named
    /// `[wall, other plane]`, and a plane that meets the lateral while the wall holds the axis can
    /// only be ⊥ to it, which is exactly what `chart_of` collects. So the arrangement **splits at
    /// every axis interval**, and inside one interval a cell is just a θ-sector: no face walk, no
    /// angular order at a vertex.
    ///
    /// ★★★★ **A refusal loses the whole chart, not one interval.** Skipping the interval that
    /// refused would leave a *partial* chart that the census then compares as if it were whole —
    /// the exact shape in which a counting claim goes blind. ☑ Measured 0 refusals; the shape is
    /// not left to that.
    pub(crate) fn cells(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
    ) -> Option<Vec<Cell>> {
        let mut out = Vec::new();
        for i in 0..self.z_lines.len().saturating_sub(1) {
            let (lo, hi) = (self.z_lines[i].t, self.z_lines[i + 1].t);
            // Alive over the *whole* interval: a ruling that stopped inside one would mean an
            // endpoint that is not an axis line, which is the premise above.
            let alive: Vec<usize> = (0..self.theta.len())
                .filter(|&j| self.theta[j].z[0] <= lo && hi <= self.theta[j].z[1])
                .collect();
            if alive.is_empty() {
                out.push(Cell {
                    interval: i,
                    sector: 0,
                    walls: None,
                });
                continue;
            }
            // ★★★★★ **One node per ruling.** A ruling is vertical, so its two ends share a θ —
            // handing both to `circular_order` asks it to rank one point twice, and it honestly
            // refuses that as the two-names-for-one-point it checks for. The order is a property
            // of the *segment*, and one endpoint states it.
            let nodes: Vec<combinatorics::NodeId> =
                alive.iter().map(|&j| self.theta[j].end[0]).collect();
            let (order, _) = crate::arrangement::circular_order(jd, k, def, &nodes).ok()?;
            let n = order.len();
            for s in 0..n {
                // `k` cuts make `k` sectors, and one cut makes one: the circle opens at that point
                // and closes at it again.
                out.push(Cell {
                    interval: i,
                    sector: s,
                    walls: Some([alive[order[s]], alive[order[(s + 1) % n]]]),
                });
            }
        }
        Some(out)
    }

    /// A ruling's **name** — `(wall class, root)`, the key `panel_faces::name_on` builds on the
    /// other side of the census's join. Derived from the segment rather than stored beside it.
    pub(crate) fn ruling_name(&self, i: usize) -> Option<RulingName> {
        let t = self.theta.get(i)?;
        let (_, _, root) = combinatorics::branch_name(t.end[0])?;
        Some((t.wall, root))
    }
}

/// **Every ⊥ class that cuts this cylinder, and every ruling on it** — per cylinder *class*.
///
/// ★★★★★ **The sources are `bands_of`'s, and one thing it does is deliberately *not* done here.**
/// That function ends with `keyed.retain(|(t, _)| t ∈ row.span)`: it clips to the **face** it was
/// asked about, because its row is a face and not a class. A class has no such span — one lateral
/// surface can carry several faces with gaps between them, which is exactly why `CylRow` became
/// per-face. So the clip is dropped, and that is a *repair* rather than a loss: `CylRow`'s own doc
/// warns that merging the spans into one `min..max` would **invent a band where the solid has no
/// face at all**, and on a chart that gap is simply a cell that keeps nothing — the argument the
/// plane side already makes for its unbounded cells.
///
/// ★ The cylinder class itself is one solid's surface, not both operands': `planes.rs` interns
/// cylinders by `Handle<Surface>` (planes by geometry), and two solids share no handles. That is
/// sound here only because a **coincident pair is refused** — `cylinders_clear` turns coaxial
/// cylinders away as `CylinderPairContact`. What does see both operands is what the chart is made
/// *of*: the plane classes, which are interned geometrically.
pub(crate) fn chart_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    k: usize,
    plane_faces: &[LocalFace],
    curved: &Curved,
) -> Result<Chart, BoolError> {
    let def = &cyls[k].def;
    let mut classes: Vec<usize> = Vec::new();
    let push = |c: usize, v: &mut Vec<usize>| {
        if !v.contains(&c) {
            v.push(c);
        }
    };
    // (i) the circles the plane arrangement actually emitted, and (ii) the cut rims, which emit
    // arcs instead and so are invisible to (i). `bands_of` reads both for the same reason.
    for lf in plane_faces {
        let ClassIx::Plane(c) = lf.surf else { continue };
        if std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, Bound::Circle { cyl } if *cyl == k))
        {
            push(c, &mut classes);
        }
    }
    for &(kk, c) in curved.cut_rims.keys() {
        if kk == k {
            push(c, &mut classes);
        }
    }
    // (iii) every ⊥ class, full stop. `bands_of` only takes the ones that land on a face's rim
    // because it is answering about a face; a chart covers the whole axis.
    for c in 0..jd.planes.len() {
        let Some(coeffs) = combinatorics::class_coeffs_rat(jd, c) else {
            continue;
        };
        if nacre_scalar::parallel_rat(&[coeffs[0], coeffs[1], coeffs[2]], &def.dir()) {
            push(c, &mut classes);
        }
    }
    let mut z_lines: Vec<ZLine> = classes
        .into_iter()
        .map(|class| crate::bands::axis_param(jd, class, def).map(|t| ZLine { class, t }))
        .collect::<Result<_, BoolError>>()?;
    z_lines.sort_by_key(|l| l.t);
    // ★ Deduped by **parameter**, and that is exact rather than lucky: two classes ⊥ to this axis
    // at one `t` are the *same plane*, which is `bands_of`'s own words for why it matches span ends
    // against `t` without a tolerance. So one of the two indices survives arbitrarily — which is
    // why everything downstream compares `t`, never the class index.
    z_lines.dedup_by_key(|l| l.t);

    let theta: Vec<ThetaSeg> = curved
        .rulings
        .get(&k)
        .map(|v| {
            v.iter()
                .map(|r: &RulingExtent| {
                    let (lo, hi) = if r.z[0] <= r.z[1] {
                        (r.z[0], r.z[1])
                    } else {
                        (r.z[1], r.z[0])
                    };
                    ThetaSeg {
                        wall: r.wall,
                        end: r.end,
                        z: [lo, hi],
                        label: r.label,
                        marks: r.marks.clone(),
                    }
                })
                .collect()
        })
        .unwrap_or_default();
    Ok(Chart { z_lines, theta })
}

/// **What the chart's line set looks like, beside what the hand-written roads answer today.**
///
/// The first rung of capability D ships no capability: it builds the two axes and *measures* them.
/// The numbers decide the next rung's design — above all whether the cells are worth walking or
/// worth overlaying — and until they are in, that design would rest on an unmeasured premise.
///
/// ★ Called from the one place the plane arrangement's faces and [`Curved`] are both in hand. It
/// reads; nothing downstream reads it.
pub(crate) fn census(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
    lateral: Option<&[LocalFace]>,
) {
    for (k, _) in cyls.iter().enumerate() {
        // ★ A class whose axis parameter cannot be stated is skipped, not recorded: `bands_of`
        // refuses on the very same `axis_param` a few lines further down, so this is not an
        // outcome the chart introduces. ☑ Measured **0** across the suite either way.
        let Ok(chart) = chart_of(jd, cyls, k, plane_faces, curved) else {
            continue;
        };
        // How today's road sees the same cylinder: one row per lateral *face*, each with its own
        // clipped span. The chart has one line set for the whole class.
        let mine: Vec<&crate::bands::CylRow> = rows.iter().filter(|r| r.class == k).collect();
        // ★★★★★ **The claim is asserted here, where the fact is made.** A test that reads the
        // ledger afterwards sees only what ran before it — this session measured that exact hole
        // (a census reported 764 of 2174 because it sat mid-suite). Asserting at the record makes
        // the coverage total and names the offending test in the panic.
        //
        // ★ The relation is **⊇, not equality**: a class's z-lines include every ⊥ plane, whether
        // or not a face reaches it, while a row's band boundaries are clipped to that face. So
        // what is checked is *cover*, and the surplus is what the chart gains — the gaps a row
        // cannot describe at all.
        //
        // ★★★★★ **The comparison asks `bands_of` itself, not the row's two ends.** The first
        // spelling checked `r.span[0]`/`r.span[1]` and the message still claimed "a band
        // boundary" — a weak check wearing a strong claim, which is the shape this session kept
        // catching. A row's boundaries come from three sources (emitted circle carriers, cut rims,
        // span ends), so only that function can say what they are.
        //
        // ★★ **Compared by `t`, never by class index.** The chart dedups its z-lines by axis
        // parameter, which is sound because two ⊥ classes at one `t` *are* one plane (`bands_of`
        // says so in its own words) — but it means the surviving line may carry the other class's
        // index. Matching on the index would fail for a reason that is not a lost boundary.
        //
        // ★★★★ **A refusal here is skipped, never a panic.** `bands_of` runs for real a few lines
        // below and may honestly refuse (a class whose axis parameter is too wide for `Rat`);
        // this census sits *before* it, so panicking would turn a legitimate reject into a crash
        // — an instrument changing the answer it came to watch. ☑ Measured 0 either way.
        let mut intervals = 0usize;
        let mut intervals_split = 0usize;
        for r in &mine {
            let Ok(bands) = crate::bands::bands_of(r, plane_faces, jd, &curved.cut_rims) else {
                continue;
            };
            for (lo, hi) in bands {
                for c in [lo, hi] {
                    let Ok(t) = crate::bands::axis_param(jd, c, &r.def) else {
                        continue;
                    };
                    assert!(
                        chart.z_lines.iter().any(|l| l.t == t),
                        "a chart lost a band boundary its own class's row had: cyl {k}, \
                         class {c}, z-lines {:?}",
                        chart.z_lines.iter().map(|l| l.class).collect::<Vec<_>>()
                    );
                }
                // How much finer the chart sees this same interval — the `545 → 1/2/3` split the
                // rung's design was chosen on, recounted here in the committed code.
                intervals += 1;
                if let (Ok(a), Ok(b)) = (
                    crate::bands::axis_param(jd, lo, &r.def),
                    crate::bands::axis_param(jd, hi, &r.def),
                ) {
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    let inside = chart.z_lines.iter().filter(|l| l.t > a && l.t < b).count();
                    if inside > 0 {
                        intervals_split += 1;
                    }
                }
            }
        }

        // ★★★★★ **Every vertical line carries an answer** (capability D, D2a). A ruling without
        // one is a chart that can state where a wall crosses but not what changes across it — and
        // the census below would then be comparing a partial chart. ☑ Measured 862/862 across the
        // suite; `world_rat_sense` declines none of them.
        for t in &chart.theta {
            assert!(
                t.label.is_some(),
                "a ruling reached the chart with no label: cyl {k}, wall {}",
                t.wall
            );
        }

        // ★★★★★ **The premise `cells()` is built on, checked where the fact is made.** That
        // function keeps a ruling only where it spans an axis interval *whole* — so a ruling
        // ending strictly inside one would be **silently dropped**, and the cells would come out
        // wrong with nothing red. The premise is that a ruling's ends are always z-lines, which is
        // structural (a ruling node is `[wall, other plane]` and that other plane can only be ⊥ to
        // the axis, which `chart_of` collects in full) — but it joins two *different* spellings of
        // the axis parameter (`node_axis_param` for the ends, `bands::axis_param` for the lines),
        // and the day those disagree the drop is what would happen. ☑ Measured 0 across the suite.
        for t in &chart.theta {
            for e in t.z {
                assert!(
                    chart.z_lines.iter().any(|l| l.t == e),
                    "a ruling ends off the chart's own z-lines, so a cell would lose it: \
                     cyl {k}, wall {}",
                    t.wall
                );
            }
        }

        let def = &cyls[k].def;
        let Some(cells) = chart.cells(jd, k, def) else {
            probe::push(probe::Row {
                z_lines: chart.z_lines.len(),
                theta: chart.theta.len(),
                rows: mine.len(),
                refused: true,
                ..probe::Row::empty()
            });
            continue;
        };

        // How many rulings live over each interval — `k = 1` opens the circle at one point and
        // makes **one** cell, and an odd `k` is the only way to reach that branch.
        let mut whole_circle = 0usize;
        let mut odd_k = 0usize;
        // ── The vertical answer, read (D2a). ──
        //
        // ★★★★★ **The consistency check is «sign-free», and deliberately so.** Asking "which side
        // of the wall is this sector" would need a third sign beside the two D2a-a already
        // composes, and this ladder's defect shape is a sign spelled a second time. A ruling's
        // label states the chamber on **both** sides of its wall, so what each ruling contributes
        // to a walk around the interval is just *whether the two differ* — and going all the way
        // round must come back to where it started. That is `label_cells`' final verification,
        // stated on the chart with no orientation at all.
        let mut rulings = 0usize;
        let mut wall_flips = 0usize;
        let mut intervals_with_flip = 0usize;
        let mut closes = 0usize;
        let mut does_not_close = 0usize;
        let mut grazing_rulings = 0usize;
        for i in 0..chart.z_lines.len().saturating_sub(1) {
            let (lo, hi) = (chart.z_lines[i].t, chart.z_lines[i + 1].t);
            let alive: Vec<&ThetaSeg> = chart
                .theta
                .iter()
                .filter(|t| t.z[0] <= lo && hi <= t.z[1])
                .collect();
            let n = alive.len();
            if n == 0 {
                whole_circle += 1;
                continue;
            }
            if n % 2 == 1 {
                odd_k += 1;
            }
            rulings += n;
            // Does crossing this wall change the material at the lateral? `Label` is
            // `[A_above, A_below, B_above, B_below]`, so the two sides are the even/odd halves.
            let flip = |l: crate::arrangement::Label| [l[0] != l[1], l[2] != l[3]];
            let mut acc = [false; 2];
            let mut all_known = true;
            let mut any_flip = false;
            for t in &alive {
                // ★★ **Existence before membership** — cell ㉒'s split, on the vertical axis. The
                // chart holds every wall's ruling, whether or not *this* lateral face reaches it,
                // and a label only ever answers membership. `face_spans` is the one spelling of
                // the rule and it reads exactly this list.
                // ★★★★★ **The first spelling of this was vacuous.** It asked whether the mark
                // list was *empty*, and a `MergedRuling` exists only because a lateral traced it —
                // so the answer was `0` because nothing was being looked at, not because the
                // population is empty. The signal `face_spans` actually reads is the **kind**: a
                // face running through says `Transversal`, one whose boundary stops at the line
                // says `Graze` and may not reach the interval at all.
                // ☑ Measured over the suite: 412 `Transversal`, **12 `Graze`** — not empty.
                if t.marks
                    .iter()
                    .all(|(_, k)| matches!(k, crate::arrangement::SegKind::Graze { .. }))
                {
                    grazing_rulings += 1;
                }
                match t.label {
                    Some(l) => {
                        let f = flip(l);
                        any_flip |= f[0] || f[1];
                        acc = [acc[0] ^ f[0], acc[1] ^ f[1]];
                    }
                    None => all_known = false,
                }
            }
            if any_flip {
                intervals_with_flip += 1;
            }
            wall_flips += alive
                .iter()
                .filter(|t| t.label.is_some_and(|l| flip(l)[0] || flip(l)[1]))
                .count();
            if all_known {
                if acc == [false; 2] {
                    closes += 1;
                } else {
                    does_not_close += 1;
                }
                // ★★★★★ **`label_cells`' final verification, on the cylinder chart.** The plane
                // side ends its labelling by checking every edge's flip relation and calling a
                // failure `LabelConflict`; this is that check, stated where no orientation is
                // needed — walk the circle, XOR what each wall changes, and come back to where you
                // started. ☑ Measured 201 intervals, **none** failing to close.
                // ★ It is asserted (not merely counted) because a failure would mean the chart's
                // own labels contradict each other, which is a defect on this side and not a
                // disagreement with today's road.
                //
                // ★★★★★ **What it cannot see, measured rather than assumed.** A walk that XORs is
                // blind to any error appearing an **even** number of times around the circle — and
                // every interval here carries an even count of rulings (`odd_k`, measured 0), so a
                // *per-ruling* systematic flip cancels itself and passes. ☑ Probed both ways:
                // adding one flip per ruling stays green, seeding the accumulator wrong goes red.
                // So this holds the labels against each other; what holds their absolute sense is
                // `ruling_probe::SIDE_CHECK`, which compares against content and is not a walk.
                assert!(
                    acc == [false; 2],
                    "the chart's labels do not close around an interval: cyl {k}, interval {i}"
                );
            }
        }

        // ── The comparison: today's **emitted** lateral faces against the chart's cells. ──
        //
        // ★ `None` means the band pass declined for this operand pair. The chart is still recorded
        // — it is the arrangement's property, not that pass's — but there is nothing to compare it
        // against, and the row says so by leaving the comparison counters at zero.
        let Some(lateral) = lateral else {
            probe::push(probe::Row {
                z_lines: chart.z_lines.len(),
                theta: chart.theta.len(),
                rows: mine.len(),
                cells: cells.len(),
                whole_circle,
                odd_k,
                intervals,
                intervals_split,
                no_faces: true,
                ..probe::Row::empty()
            });
            continue;
        };
        let mut claims = vec![0usize; cells.len()];
        let mut unnamed = 0usize;
        let mut reversed = 0usize;
        let mut band_multi = 0usize;
        let mut band_over_rulings = 0usize;
        // The 80's three answers: the ruling is not on this face · it is and the chamber changes
        // across it (a boundary the band road misses) · it is and nothing changes (harmless).
        let mut band_split_absent = 0usize;
        let mut band_split_missed = 0usize;
        let mut band_split_harmless = 0usize;
        // ★★★★ Counted because a green `unnamed == 0` says nothing if an arm was never walked —
        // the "is the zero a population or an untravelled road" question this ladder has already
        // been caught by once.
        let mut bands_seen = 0usize;
        let mut panels_seen = 0usize;
        let line_at = |t: Rat| chart.z_lines.iter().position(|l| l.t == t);
        for lf in lateral
            .iter()
            .filter(|f| matches!(f.surf, ClassIx::Cyl(c) if c == k))
        {
            match &lf.outer {
                Bound::Band { lo, hi } => {
                    bands_seen += 1;
                    // A band names its two rim classes; its cells are every interval between them,
                    // whole circle by whole circle.
                    let (Ok(a), Ok(b)) = (
                        crate::bands::axis_param(jd, *lo, def),
                        crate::bands::axis_param(jd, *hi, def),
                    ) else {
                        unnamed += 1;
                        continue;
                    };
                    let (Some(a), Some(b)) = (line_at(a), line_at(b)) else {
                        unnamed += 1;
                        continue;
                    };
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    // A band whose two rims land on one z-line spans no interval, so it names no
                    // cell — which is what `unnamed` says. Without this it would claim nothing and
                    // be counted as nothing, the one way a face could go missing quietly.
                    // ☑ Measured 0 across the suite.
                    if a == b {
                        unnamed += 1;
                        continue;
                    }
                    let mut hit = 0usize;
                    let mut sectors = 0usize;
                    for (ci, c) in cells.iter().enumerate() {
                        if c.interval >= a && c.interval < b {
                            claims[ci] += 1;
                            hit += 1;
                            sectors = sectors.max(c.sector + 1);
                        }
                    }
                    if hit > 1 {
                        band_multi += 1;
                    }
                    // ★★★★★ **The one thing this census can newly find.** `band_faces` takes the
                    // band road whenever the interval's two rims are not both cut — and then emits
                    // a **whole circle**, even where a wall holding the axis has cut rulings across
                    // it. Its own comment leans on "the gate keeps its boundary faces clear of the
                    // lateral", which the record-and-pass arm (`d = 0`) does not satisfy. So this
                    // counts the bands that span more than one sector: not predicted, because a
                    // zero means the population is empty and a non-zero means today's road is
                    // calling a wall-crossed strip uniform.
                    if sectors > 1 {
                        band_over_rulings += 1;
                        // ★★★★★ **The three-way split D1b could not make.** That rung asked
                        // whether the chart cuts finer than the band and had no way to say what
                        // the extra line *was*. With a vertical answer there are three, and the
                        // third is cell ㉒'s: the ruling may not be on this face at all.
                        // ★ "May not be on this face" is the **graze** kind, not an empty mark
                        // list — a `MergedRuling` is never traceless, so asking that measured
                        // nothing (caught while running this rung).
                        let mut absent = true;
                        let mut flips = false;
                        for t in chart.theta.iter().filter(|t| {
                            line_at(t.z[0]).is_some_and(|x| x < b)
                                && line_at(t.z[1]).is_some_and(|y| y > a)
                        }) {
                            if !t.marks.iter().all(|(_, k)| {
                                matches!(k, crate::arrangement::SegKind::Graze { .. })
                            }) {
                                absent = false;
                            }
                            if t.label.is_some_and(|l| l[0] != l[1] || l[2] != l[3]) {
                                flips = true;
                            }
                        }
                        if absent {
                            band_split_absent += 1;
                        } else if flips {
                            band_split_missed += 1;
                        } else {
                            band_split_harmless += 1;
                        }
                    }
                }
                Bound::Ring(ring) => {
                    panels_seen += 1;
                    // A panel ring is `[a_lo, b_lo, b_hi, a_hi]`: the two rim classes and the two
                    // ruling names come out of the nodes' own `Branch` names, which is the key
                    // `panel_faces::name_on` already builds.
                    let nodes = &ring.nodes;
                    let Some(cell) = (|| {
                        let [a_lo, b_lo, _, a_hi] =
                            *<&[combinatorics::NodeId; 4]>::try_from(nodes.get(..4)?).ok()?;
                        let (pa, _, ra) = combinatorics::branch_name(a_lo)?;
                        let (pb, _, rb) = combinatorics::branch_name(b_lo)?;
                        let (ph, _, _) = combinatorics::branch_name(a_hi)?;
                        // The wall is the class `a_lo` and `a_hi` share; the rims are the others.
                        let wall = *pa.iter().find(|c| ph.contains(c))?;
                        let rim_lo = *pa.iter().find(|c| **c != wall)?;
                        let rim_hi = *ph.iter().find(|c| **c != wall)?;
                        let wall_b = *pb.iter().find(|c| **c != rim_lo)?;
                        let i = line_at(crate::bands::axis_param(jd, rim_lo, def).ok()?)?;
                        let j = line_at(crate::bands::axis_param(jd, rim_hi, def).ok()?)?;
                        let lo = i.min(j);
                        let (na, nb) = ((wall, ra), (wall_b, rb));
                        Some((lo, na, nb))
                    })() else {
                        unnamed += 1;
                        continue;
                    };
                    let (i, na, nb) = cell;
                    // ★ The cell's rulings are indices; the panel names them. `ruling_name` is the
                    // one derivation of that name, so the join stays a single spelling.
                    let named = |c: &Cell| -> Option<[RulingName; 2]> {
                        let [x, y] = c.walls?;
                        Some([chart.ruling_name(x)?, chart.ruling_name(y)?])
                    };
                    let exact = cells
                        .iter()
                        .position(|c| c.interval == i && named(c) == Some([na, nb]));
                    let ci = match exact {
                        Some(ci) => Some(ci),
                        None => {
                            let r = cells
                                .iter()
                                .position(|c| c.interval == i && named(c) == Some([nb, na]));
                            if r.is_some() {
                                reversed += 1;
                            }
                            r
                        }
                    };
                    match ci {
                        Some(ci) => claims[ci] += 1,
                        None => unnamed += 1,
                    }
                }
                Bound::Circle { .. } => {}
            }
        }

        // ★★★★★ **The absolute claim comes first, and the counting one is only support.**
        // "every cell is claimed exactly once" stays true under a *consistently* wrong mapping —
        // `3c1fa56` measured exactly that (a global `axis_up` inversion satisfied its relative
        // locks). So what is asserted is that every emitted face's **own name** finds a cell in
        // the chart: the chart is what has to contain it, and a chart missing a cell goes red.
        assert_eq!(
            unnamed,
            0,
            "an emitted lateral face names a cell the chart does not have: cyl {k}, \
             {} z-lines, {} rulings, {} cells",
            chart.z_lines.len(),
            chart.theta.len(),
            cells.len()
        );
        // Two faces claiming one cell is a defect on either side; a cell nobody claims is a gap,
        // and counting those is this rung's output.
        assert!(
            claims.iter().all(|&n| n <= 1),
            "two emitted faces claim one cell: cyl {k}, claims {claims:?}"
        );
        probe::push(probe::Row {
            z_lines: chart.z_lines.len(),
            theta: chart.theta.len(),
            rows: mine.len(),
            refused: false,
            cells: cells.len(),
            whole_circle,
            odd_k,
            claimed: claims.iter().filter(|&&n| n == 1).count(),
            unclaimed: claims.iter().filter(|&&n| n == 0).count(),
            reversed,
            band_multi,
            band_over_rulings,
            intervals,
            intervals_split,
            no_faces: false,
            bands_seen,
            panels_seen,
        });
        probe::d2::push(probe::d2::Row {
            rulings,
            wall_flips,
            intervals_with_flip,
            closes,
            does_not_close,
            grazing_rulings,
            band_split_absent,
            band_split_missed,
            band_split_harmless,
        });
    }
}

/// The census's ledger — the same shape as `ruling_probe`/`panel_probe`: filled where the fact is
/// made, read by one test that reports it.
pub(crate) mod probe {
    use std::sync::Mutex;

    /// One chart's shape, and how today's emitted lateral faces cover its cells.
    #[derive(Clone, Copy, Debug, Default)]
    pub(crate) struct Row {
        pub(crate) z_lines: usize,
        pub(crate) theta: usize,
        /// How many **rows** the hand-written road split the same class into.
        pub(crate) rows: usize,
        /// The θ order could not be formed, so this chart has **no** cells — recorded rather than
        /// silently dropped, because a partial chart compared as a whole one is how a counting
        /// claim goes blind.
        pub(crate) refused: bool,
        pub(crate) cells: usize,
        /// Intervals with no ruling over them: the cell is the whole circle.
        pub(crate) whole_circle: usize,
        /// Intervals crossed by an **odd** number of rulings — the only way to reach the
        /// one-cut-one-cell branch. ★ A zero here means that branch is *unexercised*, not verified.
        pub(crate) odd_k: usize,
        pub(crate) claimed: usize,
        pub(crate) unclaimed: usize,
        /// Panel rings whose sector matched only with its two rulings **swapped**.
        pub(crate) reversed: usize,
        /// Emitted bands covering more than one cell.
        pub(crate) band_multi: usize,
        /// Emitted bands spanning more than one θ-**sector** — a whole circle where the chart sees
        /// a wall crossing. ★ Not predicted; see the census.
        pub(crate) band_over_rulings: usize,
        /// `bands_of` intervals, and how many of them the chart cuts further.
        pub(crate) intervals: usize,
        pub(crate) intervals_split: usize,
        /// The band pass declined, so this chart has no emitted faces to be compared against —
        /// recorded rather than dropped, because a chart belongs to the arrangement.
        pub(crate) no_faces: bool,
        /// How many emitted faces of each shape the comparison actually walked.
        pub(crate) bands_seen: usize,
        pub(crate) panels_seen: usize,
    }

    impl Row {
        pub(crate) fn empty() -> Self {
            Self::default()
        }
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    /// **Capability D's third rung has its own ledger**, deliberately not more fields on [`Row`].
    /// That one is already fourteen wide, and every field of it needs a reader or `dead_code`
    /// stops the build — which is how a table nobody can read grows contrived invariants to feed
    /// it. One instrument per rung.
    pub(crate) mod d2 {
        use std::sync::Mutex;

        /// One chart's vertical answers.
        #[derive(Clone, Copy, Debug, Default)]
        pub(crate) struct Row {
            /// Alive-ruling observations across all of this chart's intervals — the denominator
            /// the two counts below are subsets of. Kept so a ratio is never read against a
            /// denominator that lives in another table.
            pub(crate) rulings: usize,
            /// Rulings whose wall changes the material at the lateral.
            pub(crate) wall_flips: usize,
            pub(crate) intervals_with_flip: usize,
            /// Intervals whose walk around the circle returns to where it started, and those whose
            /// does not — `label_cells`' final verification, stated on the chart.
            pub(crate) closes: usize,
            pub(crate) does_not_close: usize,
            /// Rulings whose every mark is a **graze** — the face stops at the line rather than
            /// crossing it, so whether it reaches the interval is `face_spans`' question and not a
            /// label's. Existence, not membership (cell ㉒'s split, on the vertical axis).
            pub(crate) grazing_rulings: usize,
            /// D1b's 80, split three ways.
            pub(crate) band_split_absent: usize,
            pub(crate) band_split_missed: usize,
            pub(crate) band_split_harmless: usize,
        }

        pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

        pub(crate) fn push(r: Row) {
            ROWS.lock()
                .expect("the probe's lock is never held across a panic")
                .push(r);
        }
    }

    pub(crate) fn push(r: Row) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
}

#[cfg(test)]
mod tests {
    use super::probe::{ROWS, Row};
    use crate::BoolKind;
    use nacre_math::{Point3, Vector3};
    use nacre_scalar::Rat;
    use nacre_topo::Model;

    /// **The chart census is live, and every chart it records has the shape a chart must have.**
    ///
    /// The real claim — "a chart covers its own class's rows" — is asserted in `census`, where the
    /// fact is made, because a test that reads the ledger sees only what ran before it.
    ///
    /// ★★★★★ **And it may not read its *own* row out of that ledger.** The first spelling did
    /// (`rows.last()`), which is only its own under `--test-threads=1`; the workspace gate runs
    /// **parallel**, and a filtered run caught it immediately — the last row belonged to a
    /// `bands::tests` fixture (θ = 6 where a bore has 0). So what is asserted here is **universal
    /// over every recorded chart**, which no interleaving can break:
    ///
    /// * `rows >= 1` — a cylinder class exists because a face made it (`cyl_rows` refuses
    ///   otherwise), so the hand-written road always has at least one row for it.
    /// * `z_lines >= 2` — a band needs two boundaries; a lateral face's own two rims are ⊥ classes
    ///   and both are in the set. ☑ Measured minimum across the suite: exactly 2.
    /// * `theta` is **even** — a plane holding the axis cuts the lateral in *two* rulings
    ///   (`QuadRoot::{Lo, Hi}`), and nothing else contributes a vertical line. ☑ Measured: 0, 4,
    ///   6, 8.
    ///
    /// And the same for the cells (D1b), again universally:
    ///
    /// * `refused` is never set — a chart's θ orders can always be formed. ☑ Measured 0 across the
    ///   suite; asserted rather than merely counted, so the day a population appears it says so.
    /// * `cells >= 1` — `z_lines >= 2` is one interval, and an interval with no ruling is still one
    ///   cell (the whole circle).
    /// * every cell is claimed at most once, and `claimed + unclaimed == cells` — the partition is
    ///   total, so a cell can neither be lost nor counted twice by the comparison.
    /// * `reversed` is never set: a panel ring's two rulings always match the sector's own order.
    ///   ☑ Measured 0 over 65 panels, which is what makes the exact-order join sound.
    #[test]
    fn the_chart_census_is_running() {
        let before = ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .len();
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a through bore");
        let rows = ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        assert!(
            rows.len() > before,
            "a boolean with a cylinder recorded no chart"
        );
        for r in &rows {
            let Row {
                z_lines,
                theta,
                rows: n,
                refused,
                cells,
                claimed,
                unclaimed,
                reversed,
                no_faces,
                ..
            } = *r;
            assert!(n >= 1, "a class with no row: {r:?}");
            assert!(z_lines >= 2, "a band needs two boundaries: {r:?}");
            assert!(
                theta % 2 == 0,
                "a wall cuts two rulings, not {theta}: {r:?}"
            );
            assert!(!refused, "a chart's theta order could not be formed: {r:?}");
            assert!(cells >= 1, "two z-lines are one interval: {r:?}");
            assert_eq!(reversed, 0, "a panel matched only reversed: {r:?}");
            // A chart has `z_lines - 1` intervals, and an interval carries either no ruling or
            // some number of them — so these two counts are disjoint subsets of that many.
            assert!(
                r.whole_circle + r.odd_k < z_lines,
                "more intervals than the chart has: {r:?}"
            );
            assert!(
                r.intervals_split <= r.intervals,
                "an interval the chart cuts further is one of them: {r:?}"
            );
            // Spanning two sectors means spanning two cells, so the finding is a subset of the
            // coarser count — and if it ever were not, one of the two is being counted wrong.
            assert!(
                r.band_over_rulings <= r.band_multi,
                "a band over a wall covers more than one cell: {r:?}"
            );
            if !no_faces {
                assert_eq!(
                    claimed + unclaimed,
                    cells,
                    "the claim count must partition the cells: {r:?}"
                );
                // ★★★★ **Is the zero a population, or a road nobody walked?** A green
                // `unnamed == 0` in `census` says nothing about an arm that was never entered, so
                // the two are tied together here: cells get claimed only by faces the comparison
                // actually saw. ☑ Measured over the suite: 284 bands and 65 panels walked.
                assert_eq!(
                    r.bands_seen + r.panels_seen == 0,
                    claimed == 0,
                    "cells are claimed exactly when a face was walked: {r:?}"
                );
            }
        }
    }

    /// **The chart's vertical answers are read, and they close** (capability D, D2a).
    ///
    /// The real claims are asserted in `census`, where the facts are made — every ruling reaches
    /// the chart with a label, and the walk around each interval returns to where it started.
    /// What this holds is that the reading is **live** and that its counters are not vacuous.
    ///
    /// ★★★★★ **What the vertical answer settled, which no label alone could.** D1b handed on 82
    /// emitted bands that span more than one θ-sector and could not say what the extra line was.
    /// With a label on the vertical lines there are three answers, and the third is cell ㉒'s —
    /// the ruling may not be on this face at all. ☑ Measured over the suite: **absent 0 · a
    /// boundary the band road misses 0 · harmless 84**. The band road is calling those strips
    /// uniform and the chart agrees they *are*: the walls crossing them change nothing at the
    /// lateral there.
    ///
    /// ☑ Beside it, over **450 alive-ruling observations** — a piece counts once per interval it
    /// spans, so this is **not** the 862-piece denominator above: **122 whose wall does change the
    /// material** (so the counter has a population), **61 intervals** carrying at least one, **12
    /// whose every mark is a graze** — the existence question is not empty either — and **201
    /// intervals closing, 0 not**.
    ///
    /// ★★ **The first spelling of the existence counter was vacuous**: it asked whether a ruling's
    /// mark list was *empty*, and a `MergedRuling` exists only because a lateral traced it. It
    /// measured `0` because nothing was being looked at. The signal `face_spans` actually reads is
    /// the **kind** (`Graze` = the face stops here), and that has 12.
    #[test]
    fn the_charts_vertical_answers_close() {
        let before = super::probe::d2::ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .len();
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([12.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([2.0, 0.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the wall boss builds");
        let rows = super::probe::d2::ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        assert!(
            rows.len() > before,
            "a boolean recorded no vertical answers"
        );
        // ★ Universal over every recorded chart, so no interleaving can break it.
        for r in &rows {
            assert_eq!(r.does_not_close, 0, "an interval did not close: {r:?}");
            // Both counts are subsets of the same denominator, in the same table.
            assert!(r.wall_flips <= r.rulings, "more flips than rulings: {r:?}");
            assert!(
                r.grazing_rulings <= r.rulings,
                "more grazes than rulings: {r:?}"
            );
            assert!(
                r.intervals_with_flip <= r.closes + r.does_not_close,
                "more flipping intervals than walked ones: {r:?}"
            );
            assert_eq!(
                r.band_split_missed, 0,
                "a band spans a wall that changes the material: {r:?}"
            );
            assert_eq!(
                r.band_split_absent, 0,
                "a band spans a ruling no face traced: {r:?}"
            );
        }
        // ★★★★ And the counters are not vacuous — a zero above must mean "the population is
        // empty here", never "nothing was looked at". That is the mistake this rung made once.
        let sum = |f: fn(&super::probe::d2::Row) -> usize| rows.iter().map(f).sum::<usize>();
        assert!(sum(|r| r.closes) > 0, "no interval was ever walked");
        assert!(
            sum(|r| r.wall_flips) > 0,
            "no wall ever changed the material, so the check saw nothing"
        );
        assert!(
            sum(|r| r.band_split_harmless) > 0,
            "no band ever spanned a sector, so the three-way split saw nothing"
        );
    }

    /// **The derivation, on an axis that is not `+z`** — both roads, on a genuinely tilted one.
    ///
    /// ★★★★★ **Why this fixture had to be built.** Everything above is derived from the gate's
    /// branch table and depends on no axis direction, but every number beside it was measured on
    /// a corpus that is effectively axis-aligned: the kit builder raises cylinders on world XY
    /// (`+z` only), the ops corpus is `+z`/`+x`/`+y`, and all three *tilted* cylinders in it are
    /// **reject** fixtures — so no tilted axis had ever reached the chart. On `+z` against an
    /// axis-aligned box a class's axis parameter is just a `z` coordinate and the division in
    /// `bands::axis_param` is trivial; tilted, it is a real rational one.
    ///
    /// ★★★★ **Both halves are stated on the exact road, and that is the whole trick.**
    /// A *Pythagorean frame* — `u = (0.6, 0.8, 0)`, `v = (−0.48, 0.36, 0.8)` — lifts to axes that
    /// are exactly orthonormal in the rationals, so `SketchPlane::from_axes` takes the
    /// world-rational path and every face of the prism raised on it states narrow coefficients
    /// (asserted below, small integers). Its normal `(0.64, −0.48, 0.6)` is the cylinder's axis,
    /// and the frame's own `u` is a **unit rational perpendicular** to it — exactly the `ref_dir`
    /// `Model::add_cylinder_exact` requires.
    ///
    /// ★★★★★ **`Model::add_cylinder` cannot state this, and that is not a kernel limit.** That
    /// entry normalizes an `f64` axis and lifts its cap points back out of `any_perpendicular`'s
    /// computed floats, so on a tilted axis the caps land outside the narrow window and the
    /// population gate honestly declines (`CylinderGateUndecided` — measured while building this).
    /// Its own doc calls it *the test entry*; the production road is `add_cylinder_exact`, and on
    /// that road the tilted case is not the irrational case.
    ///
    /// ☑ What this establishes: the chart's two axes are built, and D1a's universal claims hold,
    /// on an axis with no zero component — through-bore (circles only) **and** a wall holding the
    /// axis (rulings).
    #[test]
    fn the_chart_stands_on_a_tilted_axis() {
        // Two runs of one fixture: the bore centred inside the prism (⊥ classes only), then
        // centred **on** a wall, which puts that wall through the axis and cuts rulings.
        // The prism is 2 x 2 x 1 and the bore has r = 1/2, so the removed volume is `pi/4`
        // whole and half of that when the axis lies in a wall.
        let bore = std::f64::consts::FRAC_PI_4;
        for (name, base, want_rulings, want_vol) in [
            ("through bore", [-0.52, 1.64, 0.2], false, 4.0 - bore),
            (
                "half bore on a wall",
                [0.08, 2.44, 0.2],
                true,
                4.0 - bore / 2.0,
            ),
        ] {
            let before = ROWS
                .lock()
                .expect("the probe's lock is never held across a panic")
                .len();
            let mut m = Model::new();
            let plane = crate::SketchPlane::from_axes(
                Point3::from_array([0.0; 3]),
                Vector3::from_array([0.6, 0.8, 0.0]),
                Vector3::from_array([-0.48, 0.36, 0.8]),
            );
            let frame = match crate::apply(
                &mut m,
                &crate::Operation::DatumPlane {
                    def: crate::DatumDef::Stated(plane),
                },
            ) {
                Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
                other => panic!("stating the tilted plane: {other:?}"),
            };
            let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
            let crate::OpOutput::Extrude { solid, faces, .. } = crate::apply(
                &mut m,
                &crate::Operation::Extrude {
                    frame,
                    profile: crate::Profile2d::polygon(vec![
                        p(0.0, 0.0),
                        p(2.0, 0.0),
                        p(2.0, 2.0),
                        p(0.0, 2.0),
                    ])
                    .expect("a square"),
                    dist: 1.0,
                },
            )
            .expect("the tilted prism") else {
                unreachable!()
            };
            // ★★★★ **The fixture qualifies its own population** (the S2 census lesson, which this
            // file's neighbours keep citing): every face must state narrow world coefficients, and
            // each must be either ⊥ or ∥ to the axis **exactly** — otherwise the gate would be
            // answering the oblique branch and the chart would never be reached.
            let axis = [0.64, -0.48, 0.6].map(|x| Rat::from_decimal(x).expect("a short decimal"));
            for &f in &faces {
                let n = m
                    .surface_name
                    .get(&m.faces.get(f).surface)
                    .and_then(|nm| nm.narrow())
                    .copied()
                    .expect("a tilted prism's face states itself");
                let n = [n[0], n[1], n[2]];
                let dot = nacre_scalar::dot_sign_rat(&n, &axis);
                assert!(
                    dot == nacre_scalar::Orient::Zero || nacre_scalar::parallel_rat(&n, &axis),
                    "a face neither ⊥ nor ∥ to the axis would be the oblique branch"
                );
            }
            let q = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            let (cyl, _) = m
                .add_cylinder_exact(
                    base.map(q),
                    axis,
                    // ★ The frame's own `u`: unit, and `u · axis = 0.384 − 0.384 = 0` exactly.
                    [q(0.6), q(0.8), q(0.0)],
                    q(0.5),
                    q(3.0),
                    None,
                )
                .expect("the exact road states a tilted cylinder");
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Cut, solid, cyl)
                .unwrap_or_else(|e| panic!("{name}: {e:?}"));
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{name}: the tilted result must be sound"
            );
            // ★★★★★ **The oracle is on this model, derived from the inputs.** Green plus a clean
            // `validate` would say only that *some* solid was built; the volume says the tilted
            // boolean removed the right material, and it is computed from the prism's and bore's
            // own dimensions rather than from anything the kernel answered.
            let [s] = out[..] else {
                panic!("{name}: one solid")
            };
            let got = nacre_props::mass_props(&m, s).expect("mass props").volume;
            assert!(
                (got - want_vol).abs() < 1e-9,
                "{name}: volume {got} vs {want_vol}"
            );
            let rows = ROWS
                .lock()
                .expect("the probe's lock is never held across a panic")
                .clone();
            assert!(rows.len() > before, "{name}: no chart was recorded");
            // ★★★ Read as a **set difference**, not `last()` (a `rows.last()` spelling had to be
            // corrected in D1a — the gate runs this suite in parallel). ★ And it stays an
            // *existential* claim over a **shared** ledger, so it can be diluted by a concurrent
            // test but never carries the weight alone: what proves this fixture right is the
            // volume above, on this model.
            let want = if want_rulings { 6 } else { 0 };
            assert!(
                rows[before..]
                    .iter()
                    .any(|r| r.z_lines == 4 && r.theta == want),
                "{name}: expected a chart with {want} rulings, got {:?}",
                &rows[before..]
            );
        }
    }
}
