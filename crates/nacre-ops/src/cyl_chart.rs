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
//! The two axes (D1a), the **cells** they cut (D1b), the vertical lines' answers (D2a), and the
//! **cell reader** (D2b-0, [`Chart::read_cell`]) — the function the cutover will call, which reads
//! a cell's chamber and existence off the *horizontal* lines alone (`DiskLabels`/`ArcLabels`
//! through `bands::read_bits`/`face_spans`) and is measured, cell by cell, against what the two
//! hand-written roads emit today. The cutover itself is not: those roads are untouched, because
//! the census has to measure what they answer *today*.

use crate::arrangement::{ArcLabel, Curved, Label, RulingExtent};
use crate::boolean::{Bound, LocalFace};
use crate::planes::{ClassIx, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, combinatorics, reject};
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
/// ★ No `side` stored: a ruling is named `(wall class, side)` by [`Chart::ruling_name`] — see
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
    #[cfg(test)]
    pub(crate) label: Option<crate::arrangement::Label>,
    /// Who traced this line — [`RulingExtent::marks`]. Membership is `label`'s question; whether
    /// this lateral face is even here is this one's.
    #[cfg(test)]
    pub(crate) marks: Vec<(crate::planes::SolidSide, crate::arrangement::SegKind)>,
}

/// One cylinder class's chart: the two axes, and the cells they cut.
#[derive(Clone, Debug, Default)]
pub(crate) struct Chart {
    pub(crate) z_lines: Vec<ZLine>,
    pub(crate) theta: Vec<ThetaSeg>,
    /// How many rulings arrived with `z` descending, so that [`chart_of`] had to swap `end`
    /// alongside `z` to keep `end[e] ↔ z[e]` true. ☑ Predicted 0 (`MergedRuling::end` is stated
    /// in ascending axis order); counted rather than assumed, and promoted to an assertion the
    /// day the count is in.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) end_swapped: usize,
}

/// A ruling's **name**: the wall class that cut it, and which of the two roots it is.
///
/// ★ One name, not two. `RulingExtent` also carries a `side`, but that is *derived* from the same
/// pair (`arrangement::ruling_side`), and a second spelling of one identity is how this ladder
/// has been bitten before. The join against `panel_faces`' rings uses this key because that is
/// the key `name_on` already builds there.
type RulingName = (usize, i8);

/// A partial run as the census records it: its first and last ruling, and its chamber.
#[cfg(test)]
type Run = (RulingName, RulingName, Option<(bool, bool)>);

/// **Which of a wall's two rulings a branch node lies on** — `arrangement::ruling_side`, the one
/// spelling the assembly's `Wall::Ruling { side }` and `panel_faces`' `side_at` already read.
///
/// ★★★★★ **Not the node's `QuadRoot`.** A `Branch` name's root is `Lo`/`Hi` along
/// `ℓ = n₁ × n₂` of *its own* plane pair, so the same physical ruling reads `Hi` where its node
/// pairs the wall with the plate's top and `Lo` where it pairs it with the boss's cap (the pair
/// order and the ⊥ normal's sense both reverse `ℓ` — `QuadRoot::canonical`'s doc). ☑ Measured on
/// the corner boss while building D2b-0: two different rulings of one wall carried one name and
/// the panel join claimed the wrong sectors. The side is a fact about the *point*, so every piece
/// of one ruling answers alike, whatever its end nodes pair it with.
fn ruling_side_of(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    wall: usize,
    n: combinatorics::NodeId,
) -> Option<i8> {
    let w = combinatorics::class_coeffs_rat(jd, wall)?;
    let (line, s) = combinatorics::branch_meet(jd, k, def, n)?;
    crate::arrangement::ruling_side(&w, def, (&line, &s))
}

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
    /// ★★★ **Indices, not names.** The first spelling kept the `(wall, side)` names here because
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

    /// A ruling's **name** — `(wall class, side)`, the key the panel join on the other side of
    /// the census builds the same way ([`ruling_side_of`]). Derived from the segment rather than
    /// stored beside it.
    pub(crate) fn ruling_name(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        i: usize,
    ) -> Option<RulingName> {
        let t = self.theta.get(i)?;
        Some((t.wall, ruling_side_of(jd, k, def, t.wall, t.end[0])?))
    }
}

/// **The plane classes where this cylinder's rim is a boundary the arrangement already made**:
/// (i) the circles a plane face emitted as a `Bound::Circle`, and (ii) the cut rims, which emit
/// arcs instead and so are invisible to (i). `bands_of` reads both for the same reason, and the
/// chart's line set, its boundary rule and its census all read this one list.
pub(crate) fn rim_classes(k: usize, plane_faces: &[LocalFace], curved: &Curved) -> Vec<usize> {
    let mut classes: Vec<usize> = Vec::new();
    let mut push = |c: usize| {
        if !classes.contains(&c) {
            classes.push(c);
        }
    };
    for lf in plane_faces {
        let ClassIx::Plane(c) = lf.surf else { continue };
        if std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, Bound::Circle { cyl } if *cyl == k))
        {
            push(c);
        }
    }
    let mut cut: Vec<usize> = curved
        .cut_rims
        .keys()
        .filter(|&&(kk, _)| kk == k)
        .map(|&(_, c)| c)
        .collect();
    cut.sort_unstable(); // a HashMap's order decides nothing downstream, but the list is stated once
    for c in cut {
        push(c);
    }
    classes
}

/// **The lines a lateral face's band may not run across** — `bands_of`'s three boundary sources,
/// as axis parameters: the rim classes ([`rim_classes`]) and every row's own span ends. A
/// band-shaped emission that continues across any other line is what `bands_of` never cut, so
/// the chart merges there and nowhere else.
pub(crate) fn boundary_lines(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
) -> Result<Vec<Rat>, BoolError> {
    let mut out: Vec<Rat> = rim_classes(k, plane_faces, curved)
        .into_iter()
        .map(|c| crate::bands::axis_param(jd, c, def))
        .collect::<Result<_, _>>()?;
    for r in rows.iter().filter(|r| r.class == k) {
        out.extend(r.span);
    }
    out.sort();
    out.dedup();
    Ok(out)
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
    let mut classes: Vec<usize> = rim_classes(k, plane_faces, curved);
    let push = |c: usize, v: &mut Vec<usize>| {
        if !v.contains(&c) {
            v.push(c);
        }
    };
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

    // ★★ **`end` travels with `z`.** `per_class` states `z[e]` as the axis parameter of `end[e]`,
    // so sorting `z` alone would silently break that pairing — and the cell reader (D2b-0) asks
    // "this ruling's node *on this z-line*" of exactly that pairing. Swapping both keeps it a
    // carried fact rather than a rule spelled a second time (`node_axis_param` is not re-asked).
    let mut end_swapped = 0usize;
    let theta: Vec<ThetaSeg> = curved
        .rulings
        .get(&k)
        .map(|v| {
            v.iter()
                .map(|r: &RulingExtent| {
                    let (z, end) = if r.z[0] <= r.z[1] {
                        (r.z, r.end)
                    } else {
                        end_swapped += 1;
                        ([r.z[1], r.z[0]], [r.end[1], r.end[0]])
                    };
                    ThetaSeg {
                        wall: r.wall,
                        end,
                        z,
                        #[cfg(test)]
                        label: r.label,
                        #[cfg(test)]
                        marks: r.marks.clone(),
                    }
                })
                .collect()
        })
        .unwrap_or_default();
    Ok(Chart {
        z_lines,
        theta,
        end_swapped,
    })
}

/// **What one z-line says to one cell** — the horizontal line's answer, as the cell reads it.
///
/// The plane arrangement wrote a label on every circle a cylinder leaves on a ⊥ class, whether or
/// not a face was kept there (`emit_faces`' `disk_labels`/`arc_labels`, collected outside the keep
/// filter). So a cell's chamber is not *decided* here; it is **read** off the line at either end,
/// exactly as `bands::chamber` reads a disk and `bands::panel_faces` reads an arc.
pub(crate) enum End<'a> {
    /// The circle on this line is uncut (or cut, but every arc reads the same for this cell): one
    /// label for the whole disk, whatever the sector.
    Disk(Label),
    /// The circle is cut, and the cell's two rulings are an **adjacent pair** of its branch nodes:
    /// this is the arc between them, the very row `panel_faces::arc_at` joins on.
    Exact(&'a ArcLabel),
    /// The circle is cut but the cell's rulings are not an adjacent node pair on it (the sector
    /// lies inside one arc, or spans several, or a ruling has no node here) — a shape that needs
    /// a θ placement this rung does not build. ☑ Counted; the population decides whether it is
    /// built (`end_other`).
    Other,
    /// No circle of this cylinder on this line at all: the line is a ⊥ class outside every
    /// lateral face's span (`circle_on_class` leaves a circle on every class *within* a span).
    NoCircle,
}

impl End<'_> {
    fn label(&self) -> Option<Label> {
        match self {
            End::Disk(l) => Some(*l),
            End::Exact(a) => Some(a.label),
            End::Other | End::NoCircle => None,
        }
    }
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
    /// Is this lateral face here at all — the existence question (cell ㉒), answered by the trace
    /// where an end is cut (`bands::face_spans`) and by the row's span where none is.
    pub(crate) present: bool,
    /// The trace says the face is here while no row's span covers the cell — the one direction
    /// of disagreement that is a defect (the other is a hole in a spanned face).
    pub(crate) exist_disagree: bool,
    /// Would the cell be a result face: present, with a chamber, and `keep` differing across the
    /// wall. `None` when the chamber is unknown.
    pub(crate) emit: Option<bool>,
}

/// The lines of one chart, looked up **by axis parameter** — never by class index. Two ⊥ classes
/// at one `t` are one plane, and `chart_of` keeps an arbitrary one of their indices, so the label
/// tables (keyed by class) are joined through `t`. Built once per chart; the cells then pay no
/// rational arithmetic per lookup.
pub(crate) struct Lines {
    by_t: std::collections::HashMap<Rat, Vec<usize>>,
}

impl Lines {
    pub(crate) fn of(
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        curved: &Curved,
    ) -> Result<Self, BoolError> {
        let mut by_t: std::collections::HashMap<Rat, Vec<usize>> = std::collections::HashMap::new();
        let keys = curved
            .disk_labels
            .keys()
            .chain(curved.arc_labels.keys())
            .chain(curved.cut_rims.keys());
        for &(kk, c) in keys {
            if kk != k {
                continue;
            }
            // ★ Not swallowed: `axis_param` refuses through `reject()`, which the reject census
            // records at the raise whether or not the caller looks. `chart_of` (iii) already asked
            // this of every ⊥ class, so a refusal here cannot be reached after a chart was built —
            // the `?` removes a silent arm rather than adding a population.
            let t = crate::bands::axis_param(jd, c, def)?;
            let v = by_t.entry(t).or_default();
            if !v.contains(&c) {
                v.push(c);
            }
        }
        Ok(Self { by_t })
    }

    /// The classes carrying a label of this cylinder at `t`.
    fn classes(&self, t: Rat) -> &[usize] {
        self.by_t.get(&t).map(Vec::as_slice).unwrap_or(&[])
    }
}

impl Chart {
    /// The node of ruling `i` **on** the z-line `t`, from the carried `end ↔ z` pairing — or
    /// `None` when the ruling merely passes `t` without a node there.
    fn node_on(&self, i: usize, t: Rat) -> Option<combinatorics::NodeId> {
        let s = &self.theta[i];
        (0..2).find(|&e| s.z[e] == t).map(|e| s.end[e])
    }

    /// **The arc of a cut rim that contains the sector `[x, y)`** when the sector's rulings are not
    /// themselves an adjacent node pair of that rim — a wall whose *face* stops short of this
    /// line leaves a ruling on the chart but no node on the circle (☑ the corner boss: two chords
    /// end at the plate's corner inside the footprint, so the rim has two nodes while the chart
    /// has four rulings, and three sectors lie inside the long arc).
    ///
    /// The θ order is asked of `arrangement::circular_order`, the one spelling the cells' own
    /// order comes from — never of a coordinate. A ruling with no node here is placed by its end
    /// node on another line (a ruling is vertical, so the θ is the same); a ruling *with* a node
    /// here is that node, so no point is handed to the order twice. `None` when a rim node lies
    /// strictly inside the sector (the cell would span two arcs — no single label), when the
    /// sector is the whole circle less one ruling, or when the order cannot be formed.
    #[allow(clippy::too_many_arguments)]
    fn arc_around<'a>(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        t: Rat,
        x: usize,
        y: usize,
        rim: &crate::arrangement::CutRim,
        arcs: &'a [ArcLabel],
    ) -> Option<&'a ArcLabel> {
        let m = rim.nodes.len();
        if x == y && m >= 2 {
            return None;
        }
        let mut list: Vec<combinatorics::NodeId> = rim.nodes.clone();
        let mut index_of = |i: usize| -> Option<usize> {
            match self.node_on(i, t) {
                Some(n) => rim.nodes.iter().position(|&r| r == n),
                None => {
                    list.push(self.theta[i].end[0]);
                    Some(list.len() - 1)
                }
            }
        };
        let ix = index_of(x)?;
        let iy = if x == y { ix } else { index_of(y)? };
        let (order, _) = crate::arrangement::circular_order(jd, k, def, &list).ok()?;
        let n = order.len();
        let pos = |li: usize| order.iter().position(|&o| o == li);
        let (px, py) = (pos(ix)?, pos(iy)?);
        let is_rim = |p: usize| order[p] < m;
        // No rim node strictly inside the sector, walking CCW from x to y.
        let mut q = (px + 1) % n;
        while q != py {
            if is_rim(q) {
                return None;
            }
            q = (q + 1) % n;
        }
        // The arc's ends: the nearest rim node at or before `x` (clockwise), and at or after `y`.
        let (mut a, mut b) = (px, py);
        while !is_rim(a) {
            a = (a + n - 1) % n;
        }
        while !is_rim(b) {
            b = (b + 1) % n;
        }
        let (na, nb) = (list[order[a]], list[order[b]]);
        arcs.iter().find(|arc| arc.ends == [na, nb])
    }

    /// **Read one cell off the horizontal lines** — the function the cutover will call, measured
    /// first against what the hand-written roads emit today.
    ///
    /// No sign is derived here: `band_is_above` is `planes::plus_t_is_above` at the low end and
    /// its negation at the high end — the one spelling `bands::chamber` and `bands::panel_faces`
    /// both use — and the bits come out through `bands::read_bits`. The ruling labels (the
    /// vertical lines) are not consulted: a face that is here has a circle at both its ends, so
    /// the horizontal lines always speak, and `census` asserts that (`src0_present`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn read_cell<'a>(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        kind: crate::BoolKind,
        side: crate::planes::SolidSide,
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
                // `−t` — `chamber`'s `toward_hi(lo)` / `!toward_hi(hi)`.
                // ★ Bound to the class whose label is read, not to whichever class the loop
                // visited last: one `t` is one plane and so one class today, but the sign and
                // the label must come from the same row the day that stops being true.
                let up = crate::planes::plus_t_is_above(&jd.planes[c], def);
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
                            let b = crate::bands::read_bits(&a.label, side, above[e]);
                            same &= bits.replace(b).is_none_or(|p| p == b);
                        }
                        match (same, arcs.first()) {
                            (true, Some(a)) => End::Disk(a.label),
                            _ => End::Other,
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
                                // The split's arc order and the rim's node order are one table
                                // read twice (`emit_faces` copies `ma.end`); a missing row is the
                                // "sector labels missing" this name states (`panel_faces::arc_at`).
                                let Some(arc) = arcs.iter().find(|a| a.ends == [nx, ny]) else {
                                    return Err(reject(RejectReason::RulingBoundNotYet));
                                };
                                End::Exact(arc)
                            }
                            _ => self
                                .arc_around(jd, k, def, t[e], x, y, rim, arcs)
                                .map_or(End::Other, End::Exact),
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
            let Some(l) = ends[e].label() else { continue };
            let bits = crate::bands::read_bits(&l, side, above[e]);
            match chamber {
                None => chamber = Some(bits),
                Some(prev) if prev != bits => src2_disagree = true,
                Some(_) => {}
            }
        }
        if src2_disagree {
            chamber = None;
        }

        // ── Existence: the trace where an end is cut, the span where none is. ──
        let (lo, hi) = (t[0].min(t[1]), t[0].max(t[1]));
        let by_span = rows.iter().any(|r| {
            r.class == k && r.span[0].min(r.span[1]) <= lo && hi <= r.span[0].max(r.span[1])
        });
        let mut by_marks: Option<bool> = None;
        for e in 0..2 {
            let End::Exact(arc) = &ends[e] else { continue };
            // `face_spans` refuses by name (two of this solid's faces disagreeing on one arc), and
            // two cut ends disagreeing with each other is `panel_faces`' `CylinderFaceUndecided`:
            // the rims of a hole are band boundaries, so a sector exists over its whole height or
            // not at all, and picking an end to believe is the guess this kernel does not make.
            let v = crate::bands::face_spans(arc, side, above[e])?;
            if by_marks.replace(v).is_some_and(|p| p != v) {
                return Err(reject(RejectReason::CylinderFaceUndecided));
            }
        }
        // ★ The span is the coarser truth: a holed lateral's row spans the hole, and the trace is
        // what says the face is *not* there (cell ㉒'s population — `exist_marks_false`). So the
        // disagreement that would be a defect is only the other direction: the trace claiming a
        // face where no row spans.
        let (present, exist_disagree) = match by_marks {
            Some(v) => (v, v && !by_span),
            None => (by_span, false),
        };

        // A cell the face is not in emits nothing, chamber or no chamber — a cell outside every
        // span has no circle at either end and no question to answer.
        let emit = if !present {
            Some(false)
        } else {
            chamber.map(|(own, other)| {
                let keep = |in_own: bool| crate::bands::keep_for(kind, side, in_own, other);
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

/// **The lateral faces of every cylinder class, emitted from the chart** — the cutover's road.
///
/// The vocabulary is the band road's own: a band-shaped interval (every cell kept, one chamber)
/// is a [`Bound::Band`] between two lines, and a partial run of kept sectors is a panel ring
/// `[a_lo, b_lo, b_hi, a_hi]` walled `[Arc ccw, Ruling up, Arc cw, Ruling down]` — the ring
/// `bands::panel_faces` emits, node for node. Four rules, each the band road's, restated on the
/// chart and measured against it before this was written (D2b-1):
///
/// * **`Band` iff band-shaped and not both rims cut** — `band_faces`' dispatch. A whole-circle
///   cell with both rims cut is still a `Band`; `band_loop` refuses that one by name.
/// * **A band runs across a line only where `bands_of` never cut**: not a rim class (an emitted
///   circle or a cut rim) and not a row's span end ([`boundary_lines`]). ☑ 47 = 47 merges.
/// * **Runs split at a ruling with a node on a cut rim** — one ring per adjacent node pair, as
///   the panel road walks `CutRim.nodes`, so every ring arc is one edge. ☑ 0 changed today.
/// * **Runs are emitted by first sector** — `lo_rim.nodes` order, since both come from
///   `circular_order` seam-first; a wrap run sorts last, as today's `j = m − 1`.
///
/// No sign is derived here: chambers come from [`Chart::read_cell`], the keep rule is
/// `bands::keep_for`, the ruling identity is [`Chart::ruling_name`] (`ruling_side`), and the
/// nodes are the rulings' own ([`Chart::node_on`]).
pub(crate) fn emit_lateral(
    kind: crate::BoolKind,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
) -> Result<Vec<LocalFace>, BoolError> {
    use crate::boolean::{Ring, Wall};
    let mut out: Vec<LocalFace> = Vec::new();
    for k in 0..cyls.len() {
        let def = &cyls[k].def;
        // One class is one solid's surface (`cyl_rows` refuses a class with no row; two solids
        // share no handles) — stated as the refusal `cyl_rows` names for the wiring failure.
        let mut sides = rows.iter().filter(|r| r.class == k).map(|r| r.side);
        let Some(side) = sides.next() else {
            return Err(reject(RejectReason::CylinderGateUndecided));
        };
        if sides.any(|s| s != side) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
        let chart = chart_of(jd, cyls, k, plane_faces, curved)?;
        let lines = Lines::of(jd, k, def, curved)?;
        // The θ order the split already formed cannot fail to form again; if it does, the point's
        // own `(line, s)` could not be stated — that name's sentence.
        let cells = chart
            .cells(jd, k, def)
            .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
        let reads: Vec<CellRead<'_>> = cells
            .iter()
            .map(|c| chart.read_cell(jd, k, def, kind, side, c, curved, &lines, rows))
            .collect::<Result<_, _>>()?;
        let boundary = boundary_lines(jd, k, def, plane_faces, curved, rows)?;
        let is_cut = |t: Rat| {
            lines
                .classes(t)
                .iter()
                .any(|&c| curved.cut_rims.contains_key(&(k, c)))
        };
        let is_boundary = |t: Rat| boundary.binary_search(&t).is_ok();
        let ladder = || reject(RejectReason::RulingBoundNotYet);
        let keep = |own: bool, other: bool| crate::bands::keep_for(kind, side, own, other);

        // A band-shaped stretch not yet flushed: `(first line, last line, chamber, flip)`.
        let mut pending: Option<(usize, usize, (bool, bool), bool)> = None;
        let flush = |pending: &mut Option<(usize, usize, (bool, bool), bool)>,
                     out: &mut Vec<LocalFace>| {
            if let Some((lo, hi, _, flip)) = pending.take() {
                out.push(LocalFace {
                    surf: ClassIx::Cyl(k),
                    outer: Bound::Band {
                        lo: chart.z_lines[lo].class,
                        hi: chart.z_lines[hi].class,
                    },
                    inner: Vec::new(),
                    flip,
                });
            }
        };
        let n_int = chart.z_lines.len().saturating_sub(1);
        for i in 0..n_int {
            let idx: Vec<usize> = (0..cells.len())
                .filter(|&ci| cells[ci].interval == i)
                .collect();
            let (t_lo, t_hi) = (chart.z_lines[i].t, chart.z_lines[i + 1].t);
            // A present cell whose chamber could not be read: the two ends disagreed, which is
            // `chamber`'s own refusal.
            if idx
                .iter()
                .any(|&ci| reads[ci].present && reads[ci].emit.is_none())
            {
                return Err(reject(RejectReason::CylinderGateUndecided));
            }
            let em: Vec<bool> = idx.iter().map(|&ci| reads[ci].emit == Some(true)).collect();
            if !em.iter().any(|&e| e) {
                flush(&mut pending, &mut out);
                continue;
            }
            let n = idx.len();
            let chamber0 = reads[idx[0]].chamber;
            let all_kept = em.iter().all(|&e| e);
            let one_chamber = idx.iter().all(|&ci| reads[ci].chamber == chamber0);
            let whole = n == 1 && cells[idx[0]].walls.is_none();
            // `band_faces`' dispatch: a both-cut interval is the panel road's, whatever its
            // chambers say — except the whole circle, which has no sectors to be panels of and
            // goes to `band_loop` to be refused by name there.
            let bandlike = all_kept && one_chamber && (whole || !(is_cut(t_lo) && is_cut(t_hi)));
            if bandlike {
                let (own, other) =
                    chamber0.ok_or_else(|| reject(RejectReason::CylinderGateUndecided))?;
                let flip = !keep(own, other);
                match pending {
                    Some((lo, hi, ch, fl))
                        if hi == i && ch == (own, other) && fl == flip && !is_boundary(t_lo) =>
                    {
                        pending = Some((lo, i + 1, ch, fl));
                    }
                    _ => {
                        flush(&mut pending, &mut out);
                        pending = Some((i, i + 1, (own, other), flip));
                    }
                }
                continue;
            }
            flush(&mut pending, &mut out);
            // Sectors: both rims must be cut, or the chamber has no sector answer for an end.
            if !(is_cut(t_lo) && is_cut(t_hi)) {
                return Err(ladder());
            }
            // A run breaks between sector `s − 1` and `s` where either is dropped, the chamber
            // changes, or the shared ruling has a node on a rim (today's panel boundary).
            let breaks = |s: usize| -> bool {
                let (p, c) = (idx[(s + n - 1) % n], idx[s]);
                !em[s]
                    || !em[(s + n - 1) % n]
                    || reads[p].chamber != reads[c].chamber
                    || cells[c].walls.is_none_or(|[x, _]| {
                        chart.node_on(x, t_lo).is_some() || chart.node_on(x, t_hi).is_some()
                    })
            };
            let Some(start) = (0..n).find(|&s| breaks(s)) else {
                // Every sector kept, one chamber, no ruling with a rim node — but both rims are
                // cut, so a node exists on some ruling; a circle that breaks nowhere has none.
                return Err(ladder());
            };
            // Runs as `(first sector, first cell, last cell)`, then in first-sector order.
            let mut runs: Vec<(usize, usize, usize)> = Vec::new();
            let mut s = 0;
            while s < n {
                let at = (start + s) % n;
                if !em[at] {
                    s += 1;
                    continue;
                }
                let mut len = 1;
                while s + len < n && !breaks((start + s + len) % n) {
                    len += 1;
                }
                runs.push((
                    cells[idx[at]].sector,
                    idx[at],
                    idx[(start + s + len - 1) % n],
                ));
                s += len;
            }
            runs.sort_unstable_by_key(|&(sector, _, _)| sector);
            for (_, f, l) in runs {
                let (Some([a, _]), Some([_, b])) = (cells[f].walls, cells[l].walls) else {
                    return Err(ladder());
                };
                let node = |r: usize, t: Rat| chart.node_on(r, t).ok_or_else(ladder);
                let (a_lo, a_hi) = (node(a, t_lo)?, node(a, t_hi)?);
                let (b_lo, b_hi) = (node(b, t_lo)?, node(b, t_hi)?);
                let side_a = chart.ruling_name(jd, k, def, a).ok_or_else(ladder)?.1;
                let side_b = chart.ruling_name(jd, k, def, b).ok_or_else(ladder)?.1;
                let (own, other) = reads[f]
                    .chamber
                    .ok_or_else(|| reject(RejectReason::CylinderGateUndecided))?;
                out.push(LocalFace {
                    surf: ClassIx::Cyl(k),
                    outer: Bound::Ring(Ring::new(
                        vec![a_lo, b_lo, b_hi, a_hi],
                        vec![
                            Wall::Arc { cyl: k, ccw: true },
                            Wall::Ruling {
                                cyl: k,
                                side: side_b,
                                up: true,
                            },
                            Wall::Arc { cyl: k, ccw: false },
                            Wall::Ruling {
                                cyl: k,
                                side: side_a,
                                up: false,
                            },
                        ],
                    )),
                    inner: Vec::new(),
                    flip: !keep(own, other),
                });
            }
        }
        flush(&mut pending, &mut out);
    }
    Ok(out)
}

/// **What the chart's line set looks like, beside what the hand-written roads answer today.**
///
/// The first rung of capability D ships no capability: it builds the two axes and *measures* them.
/// The numbers decide the next rung's design — above all whether the cells are worth walking or
/// worth overlaying — and until they are in, that design would rest on an unmeasured premise.
///
/// ★ Called from the one place the plane arrangement's faces and [`Curved`] are both in hand. It
/// reads; nothing downstream reads it.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn census(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    kind: crate::BoolKind,
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
    lateral: Option<&[LocalFace]>,
    emitted: &Result<Vec<LocalFace>, BoolError>,
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
        // The cell reader's inputs (D2b-0): the lines by parameter, and whose solid this class
        // is. ★ One class is one solid's surface — `planes.rs` interns cylinders by handle and
        // two solids share none — asserted where it is relied on rather than assumed.
        let Ok(lines) = Lines::of(jd, k, &cyls[k].def, curved) else {
            continue;
        };
        // `bands_of`'s boundaries, for the merge counters; a refusal is a skip, as for the chart.
        let Ok(boundary) = boundary_lines(jd, k, &cyls[k].def, plane_faces, curved, rows) else {
            continue;
        };
        // ── The cutover's emitter against the reference road, face by face (D2b shadow). ──
        //
        // The zip is the claim: same faces, same order (handles are minted in face order — the
        // replay contract), bands compared by their lines' parameters and rings node for node.
        // A refusal on the chart side where the reference emitted is the one outcome that must
        // never happen; a refusal on the reference side where the chart emits is a capability
        // opening, counted.
        let of_k = |faces: &[LocalFace]| -> Vec<LocalFace> {
            faces
                .iter()
                .filter(|f| matches!(f.surf, ClassIx::Cyl(c) if c == k))
                .cloned()
                .collect()
        };
        let mut shadow = (0usize, 0usize, 0usize); // (compared, opened, both_refused)
        match (lateral.map(of_k), emitted.as_ref().ok().map(|v| of_k(v))) {
            (Some(r), Some(n)) => {
                assert_eq!(
                    r.len(),
                    n.len(),
                    "the chart emits a different number of lateral faces than the band road: \
                     cyl {k}"
                );
                for (fi, (a, b)) in r.iter().zip(n.iter()).enumerate() {
                    assert_eq!(a.flip, b.flip, "face {fi} of cyl {k}: flip differs");
                    assert!(
                        a.inner.is_empty() && b.inner.is_empty(),
                        "face {fi} of cyl {k}: a lateral face with inner bounds"
                    );
                    let t_of = |c: usize| crate::bands::axis_param(jd, c, &cyls[k].def).ok();
                    match (&a.outer, &b.outer) {
                        (Bound::Band { lo, hi }, Bound::Band { lo: l2, hi: h2 }) => {
                            assert!(
                                t_of(*lo) == t_of(*l2) && t_of(*hi) == t_of(*h2),
                                "face {fi} of cyl {k}: bands on different lines"
                            );
                        }
                        (Bound::Ring(x), Bound::Ring(y)) => {
                            assert_eq!(x.nodes, y.nodes, "face {fi} of cyl {k}: ring nodes differ");
                            assert_eq!(x.walls, y.walls, "face {fi} of cyl {k}: ring walls differ");
                        }
                        (x, y) => panic!("face {fi} of cyl {k}: shapes differ: {x:?} vs {y:?}"),
                    }
                }
                shadow.0 = r.len();
            }
            (None, Some(_)) => shadow.1 = 1,
            (Some(_), None) => panic!(
                "the chart road refuses what the band road emits: cyl {k}: {:?}",
                emitted.as_ref().err()
            ),
            (None, None) => shadow.2 = 1,
        }
        let side = mine.first().map(|r| r.side);
        assert!(
            mine.iter().all(|r| Some(r.side) == side),
            "one cylinder class carries rows of both solids: cyl {k}"
        );
        let mut split_flip = 0usize;
        let mut split_nocircle = 0usize;
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
                    // ★★★★★ **D1b's 44, judged (D2b-0).** A line the chart sees inside a band
                    // is not asked what *kind* it is but whether **crossing it changes the
                    // chamber** — the label's two halves differ for either solid. A line with a
                    // circle but no flip is a harmless subdivision; one with no circle at all is
                    // a line inside a face's span that `circle_on_class` did not mark, which its
                    // own contract says cannot happen (☑ predicted 0, counted).
                    for l in chart.z_lines.iter().filter(|l| l.t > a && l.t < b) {
                        let mut seen = false;
                        let mut flips = false;
                        for &c in lines.classes(l.t) {
                            let flip = |l: &crate::arrangement::Label| l[0] != l[1] || l[2] != l[3];
                            if let Some(d) = curved.disk_labels.get(&(k, c)) {
                                seen = true;
                                flips |= flip(d);
                            }
                            if let Some(arcs) = curved.arc_labels.get(&(k, c)) {
                                seen = true;
                                flips |= arcs.iter().any(|a| flip(&a.label));
                            }
                        }
                        if !seen {
                            split_nocircle += 1;
                        } else if flips {
                            split_flip += 1;
                        }
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

        // ── The cells, read off the horizontal lines (D2b-0). ──
        //
        // ★ `cyl_rows` refuses a class with no row before this census runs, so a class here always
        // has a side; stated as a check rather than an `unwrap` so the message names the class.
        let Some(side) = side else {
            panic!("a cylinder class reached the chart with no row: cyl {k}");
        };
        // A cell the reader refuses is a refusal of the whole class (the emitter's `?`); the
        // chart is still recorded, with the comparison empty and the refusal counted.
        let reads: Vec<CellRead<'_>> = match cells
            .iter()
            .map(|c| chart.read_cell(jd, k, def, kind, side, c, curved, &lines, rows))
            .collect::<Result<Vec<_>, BoolError>>()
        {
            Ok(v) => v,
            Err(_) => {
                probe::d2b::push(probe::d2b::Row {
                    cells: cells.len(),
                    no_faces: true,
                    read_refused: 1,
                    ..probe::d2b::Row::default()
                });
                continue;
            }
        };
        let mut d2b = probe::d2b::Row {
            cells: cells.len(),
            compared: shadow.0,
            opened: shadow.1,
            both_refused: shadow.2,
            end_swapped: chart.end_swapped,
            split_flip,
            split_nocircle,
            ..probe::d2b::Row::default()
        };
        for r in &reads {
            for e in &r.ends {
                match e {
                    End::Disk(_) => d2b.end_disk += 1,
                    End::Exact(_) => d2b.end_exact += 1,
                    End::Other => d2b.end_other += 1,
                    End::NoCircle => d2b.end_nocircle += 1,
                }
            }
            let sources = r
                .ends
                .iter()
                .filter(|e| matches!(e, End::Disk(_) | End::Exact(_)))
                .count();
            if sources == 0 && r.present {
                d2b.src0_present += 1;
            }
            if r.present && r.ends.iter().any(|e| matches!(e, End::Other)) {
                d2b.other_present += 1;
            }
            d2b.src2_disagree += usize::from(r.src2_disagree);
            d2b.exist_disagree += usize::from(r.exist_disagree);
            // Both ends cut is exactly the panel road's population, where today `face_spans`
            // drops a sector for existence — the one count with a twin on the other side.
            if r.ends.iter().all(|e| matches!(e, End::Exact(_))) && !r.present {
                d2b.exist_marks_false += 1;
            }
            match r.emit {
                Some(true) => d2b.emit += 1,
                None => d2b.emit_unknown += 1,
                Some(false) => {}
            }
        }
        assert_eq!(
            d2b.end_disk + d2b.end_exact + d2b.end_other + d2b.end_nocircle,
            2 * cells.len(),
            "every cell has two ends: cyl {k}"
        );
        // ★★★★★ **The horizontal lines always speak for a face that is there** — asserted where
        // the fact is made, and structural rather than measured: `circle_on_class` leaves a circle
        // on every ⊥ class within a lateral face's span (Crosses inside, Grazes at the rims), and
        // `emit_faces` labels every circle cell and arc outside the keep filter. So a present cell
        // with no label at either end is a broken premise upstream, not a chart shape — and it is
        // why no "which side of the wall" sign is needed to read a cell.
        assert_eq!(
            d2b.src0_present, 0,
            "a cell with a face has no label at either end: cyl {k}"
        );
        // ★ Promoted from counts to record-site assertions once the suite measured them 0
        // (D1b's `unnamed == 0` discipline): a reporting test sees only the rows recorded before
        // it, an assertion here sees every chart and names the offending test.
        assert_eq!(
            d2b.end_swapped, 0,
            "a ruling arrived with z descending: cyl {k}"
        );
        assert_eq!(
            d2b.src2_disagree, 0,
            "a cell's two ends disagree about its chamber: cyl {k}"
        );
        assert_eq!(
            d2b.split_nocircle, 0,
            "a line inside a band has no circle of this cylinder: cyl {k}"
        );

        // ── Neighbours with one chamber: where the cutover would have to merge. ──
        let n_int = chart.z_lines.len().saturating_sub(1);
        let idx_of = |i: usize| -> Vec<usize> {
            (0..cells.len())
                .filter(|&ci| cells[ci].interval == i)
                .collect()
        };
        let emitted = |ci: usize| reads[ci].emit == Some(true);
        let same_chamber = |a: usize, b: usize| reads[a].chamber == reads[b].chamber;
        for i in 0..n_int {
            let idx = idx_of(i);
            let n = idx.len();
            // θ-neighbours: sectors are in circular order, and two sectors make one pair.
            if n >= 2 {
                let pairs = if n == 2 { 1 } else { n };
                for s in 0..pairs {
                    let (a, b) = (idx[s], idx[(s + 1) % n]);
                    if emitted(a) && emitted(b) && same_chamber(a, b) {
                        d2b.theta_merge_pairs += 1;
                    }
                }
            }
            // A disk-ended interval whose emitted sectors do not close into one whole circle
            // would need a ring with arcs on an *uncut* circle — no `Bound` spells that today.
            let disk_ended = idx
                .iter()
                .any(|&ci| reads[ci].ends.iter().any(|e| matches!(e, End::Disk(_))));
            if disk_ended && n >= 2 {
                let em: Vec<usize> = idx.iter().copied().filter(|&ci| emitted(ci)).collect();
                if let Some(&first) = em.first()
                    && (em.len() != n || em.iter().any(|&ci| !same_chamber(ci, first)))
                {
                    d2b.partial_theta_in_disk_interval += 1;
                }
            }
        }
        // ── What the cutover would emit per interval, measured before it is built (D2b-1 census). ──
        //
        // A run = a maximal circular stretch of emitted sectors with one chamber. An interval that
        // is one chamber all the way round is **band-shaped** (today's `Band`); any other run is a
        // panel ring `[a_lo, b_lo, b_hi, a_hi]`, which needs the run's two boundary rulings to
        // have a **node on each cut rim** — a cut rim without that node has no vertex to close the
        // ring on. "Cut" is one spelling throughout: a `cut_rims` key at that line.
        let is_cut = |t: Rat| {
            lines
                .classes(t)
                .iter()
                .any(|&c| curved.cut_rims.contains_key(&(k, c)))
        };
        let is_boundary = |t: Rat| boundary.binary_search(&t).is_ok();
        // Per interval: the chamber if band-shaped (every cell emitted, one chamber), else `None`.
        let mut bandlike: Vec<Option<(bool, bool)>> = vec![None; n_int];
        // Per interval: its partial runs as `(first ruling name, last ruling name, chamber)`.
        let mut runs: Vec<Vec<Run>> = vec![Vec::new(); n_int];
        for i in 0..n_int {
            let idx = idx_of(i);
            let n = idx.len();
            if n == 0 {
                continue;
            }
            let cut_end = |ci: usize, e: usize| matches!(reads[ci].ends[e], End::Exact(_));
            let (t_lo, t_hi) = (chart.z_lines[i].t, chart.z_lines[i + 1].t);
            // A whole-circle cell is a `Band`, which `band_loop` refuses when both rims are cut.
            if n == 1 && cells[idx[0]].walls.is_none() {
                if emitted(idx[0]) {
                    d2b.whole_emitted += 1;
                    bandlike[i] = reads[idx[0]].chamber;
                    if is_cut(t_lo) && is_cut(t_hi) {
                        d2b.whole_both_cut += 1;
                    }
                }
                continue;
            }
            let em: Vec<bool> = idx.iter().map(|&ci| emitted(ci)).collect();
            if em.iter().all(|&e| e) && (1..n).all(|s| same_chamber(idx[0], idx[s])) {
                d2b.full_runs += 1; // one chamber all the way round: a `Band`
                bandlike[i] = reads[idx[0]].chamber;
                // …unless both rims are cut: today's road sends that interval to the panel road
                // (`band_faces`' dispatch), and a `Band` there would be refused by `band_loop`.
                if is_cut(t_lo) && is_cut(t_hi) {
                    d2b.full_run_both_cut += 1;
                }
                continue;
            }
            // Runs, circularly, starting where a run cannot begin mid-way.
            let start = (0..n)
                .find(|&s| !em[s] || !same_chamber(idx[s], idx[(s + n - 1) % n]))
                .unwrap_or(0);
            let mut s = 0;
            while s < n {
                let at = (start + s) % n;
                if !em[at] {
                    s += 1;
                    continue;
                }
                let mut len = 1;
                while s + len < n {
                    let nx = (start + s + len) % n;
                    if em[nx] && same_chamber(idx[at], idx[nx]) {
                        len += 1;
                    } else {
                        break;
                    }
                }
                d2b.partial_runs += 1;
                // The same run, counted again as the panel road would cut it: a ruling with a node
                // on a cut rim ends a panel there even when the chamber continues past it.
                for j in 1..len {
                    let ci = idx[(start + s + j) % n];
                    if let Some([x, _]) = cells[ci].walls
                        && [(0, t_lo), (1, t_hi)]
                            .iter()
                            .any(|&(e, t)| cut_end(ci, e) && chart.node_on(x, t).is_some())
                    {
                        d2b.run_split_at_node += 1;
                    }
                }
                let first = idx[at];
                let last = idx[(start + s + len - 1) % n];
                match (cells[first].walls, cells[last].walls) {
                    (Some([x, _]), Some([_, y])) => {
                        for (ci, r) in [(first, x), (last, y)] {
                            for (e, t) in [(0, t_lo), (1, t_hi)] {
                                if cut_end(ci, e) && chart.node_on(r, t).is_none() {
                                    d2b.run_boundary_no_node += 1;
                                }
                            }
                        }
                        if let (Some(na), Some(nb)) = (
                            chart.ruling_name(jd, k, def, x),
                            chart.ruling_name(jd, k, def, y),
                        ) {
                            runs[i].push((na, nb, reads[first].chamber));
                        }
                    }
                    _ => d2b.run_boundary_no_node += 1,
                }
                s += len;
            }
            if runs[i].len() > 1 {
                d2b.intervals_multi_run += 1;
            }
        }
        // ── Across a z-line: what the cutover would merge at the chart, and what it may not. ──
        //
        // Two band-shaped neighbours with one chamber merge across their shared line **unless** it
        // is one of `bands_of`'s boundaries (`boundary_lines`). The sub-count across a line a
        // plane face emitted as a `Bound::Circle` is asserted 0 where it is made: a whole circle
        // through a present band with the same chamber on both sides is a T-junction, which the
        // arrangement does not produce.
        let mut z_merge_bandlike_circle = 0usize;
        for i in 1..n_int {
            let t = chart.z_lines[i].t;
            if let (Some(p), Some(q)) = (bandlike[i - 1], bandlike[i])
                && p == q
            {
                if !is_boundary(t) {
                    d2b.z_merge_bandlike += 1;
                } else if rim_classes(k, plane_faces, curved).iter().any(|&c| {
                    lines.classes(t).contains(&c) && !curved.cut_rims.contains_key(&(k, c))
                }) {
                    z_merge_bandlike_circle += 1;
                }
            }
            // Two partial runs with the same boundary rulings and chamber across a non-boundary
            // line: a ring the cutover would have to stretch across z — no arm is written for it
            // until this counts one.
            if !is_boundary(t) {
                for a in &runs[i - 1] {
                    for b in &runs[i] {
                        if a.0 == b.0 && a.1 == b.1 && a.2 == b.2 {
                            d2b.run_run_nonboundary += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(
            z_merge_bandlike_circle, 0,
            "two band-shaped cells with one chamber meet across an emitted circle: cyl {k}"
        );
        // z-neighbours: the same sector (by ruling *names* — indices differ per interval) in
        // consecutive intervals — the cell-level count, kept as the finer view.
        let sector_name = |ci: usize| -> Option<[Option<RulingName>; 2]> {
            match cells[ci].walls {
                None => Some([None, None]),
                Some([x, y]) => Some([
                    Some(chart.ruling_name(jd, k, def, x)?),
                    Some(chart.ruling_name(jd, k, def, y)?),
                ]),
            }
        };
        for i in 1..n_int {
            for &a in &idx_of(i - 1) {
                for &b in &idx_of(i) {
                    let same = match (sector_name(a), sector_name(b)) {
                        (Some(p), Some(q)) => p == q,
                        _ => false,
                    };
                    if same && emitted(a) && emitted(b) && same_chamber(a, b) {
                        d2b.z_merge_pairs += 1;
                    }
                }
            }
        }

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
            // The reads stand on their own (they need no emitted face); only the comparison
            // against today's emission is empty here, and the row says so.
            d2b.no_faces = true;
            probe::d2b::push(d2b);
            continue;
        };
        let mut claims = vec![0usize; cells.len()];
        let mut ref_band_merges = 0usize;
        let mut ref_merge_over_boundary = 0usize;
        // The order today's road emits faces in, as the first cell each claims — the chart's
        // cells are in `(interval, sector)` order, so a descent here is an order the cutover
        // would change (handles are minted in face order: the replay contract).
        let mut first_claim: Vec<usize> = Vec::new();
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
                    // The merges today's road made inside this band: one per interior line. Each
                    // must be a line the chart's boundary rule lets a band-shaped pair cross.
                    ref_band_merges += b - a - 1;
                    for i in (a + 1)..b {
                        let t = chart.z_lines[i].t;
                        if is_boundary(t) {
                            ref_merge_over_boundary += 1;
                        }
                    }
                    let mut hit = 0usize;
                    let mut sectors = 0usize;
                    for (ci, c) in cells.iter().enumerate() {
                        if c.interval >= a && c.interval < b {
                            if hit == 0 {
                                first_claim.push(ci);
                            }
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
                        let (pa, _, _) = combinatorics::branch_name(a_lo)?;
                        let (pb, _, _) = combinatorics::branch_name(b_lo)?;
                        let (ph, _, _) = combinatorics::branch_name(a_hi)?;
                        // The wall is the class `a_lo` and `a_hi` share; the rims are the others.
                        let wall = *pa.iter().find(|c| ph.contains(c))?;
                        let rim_lo = *pa.iter().find(|c| **c != wall)?;
                        let rim_hi = *ph.iter().find(|c| **c != wall)?;
                        let wall_b = *pb.iter().find(|c| **c != rim_lo)?;
                        let i = line_at(crate::bands::axis_param(jd, rim_lo, def).ok()?)?;
                        let j = line_at(crate::bands::axis_param(jd, rim_hi, def).ok()?)?;
                        let lo = i.min(j);
                        let na = (wall, ruling_side_of(jd, k, def, wall, a_lo)?);
                        let nb = (wall_b, ruling_side_of(jd, k, def, wall_b, b_lo)?);
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
                        Some([
                            chart.ruling_name(jd, k, def, x)?,
                            chart.ruling_name(jd, k, def, y)?,
                        ])
                    };
                    // ★★★★ **A panel claims every sector it covers, not one** (corrected in
                    // D2b-0). A panel's arc runs between two *rim nodes*, but the chart's sectors
                    // are cut by every ruling alive over the interval — including rulings of walls
                    // whose face stops short of this rim and so leave no node on it. The corner
                    // boss is the population: its long arc is one panel over **three** chart
                    // sectors. D1b's join asked for the one cell whose two rulings *are* the
                    // panel's, claimed it, and counted the other two as unclaimed gaps. The sweep
                    // below claims from the sector leaving `na` round to the one arriving at `nb`.
                    let mut sweep =
                        |claims: &mut Vec<usize>, from: RulingName, to: RulingName| -> bool {
                            let ring: Vec<usize> = (0..cells.len())
                                .filter(|&ci| cells[ci].interval == i)
                                .collect();
                            let n = ring.len();
                            let Some(start) = ring
                                .iter()
                                .position(|&ci| named(&cells[ci]).is_some_and(|w| w[0] == from))
                            else {
                                return false;
                            };
                            let mut hit = Vec::new();
                            for step in 0..n {
                                let ci = ring[(start + step) % n];
                                hit.push(ci);
                                if named(&cells[ci]).is_some_and(|w| w[1] == to) {
                                    // The run's *first* sector — the statistic the emitter
                                    // orders runs by (a wrap run's minimum is the adjacent one).
                                    first_claim.push(ring[start]);
                                    for ci in hit {
                                        claims[ci] += 1;
                                    }
                                    return true;
                                }
                            }
                            false
                        };
                    if sweep(&mut claims, na, nb) {
                    } else if sweep(&mut claims, nb, na) {
                        reversed += 1;
                    } else {
                        unnamed += 1;
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
        // ── The headline (D2b-0): does the cell reader agree with today's roads, cell by cell? ──
        //
        // `claims[ci] == 1` is "an emitted band or panel names this cell"; `emit` is what the
        // reader would emit. Every disagreement is a cell the cutover would change — a forecast
        // the vertical answer (D2a) could not make without a third sign, and the horizontal one
        // makes without any.
        for (ci, r) in reads.iter().enumerate() {
            if let Some(e) = r.emit
                && e != (claims[ci] == 1)
            {
                d2b.emit_mismatch += 1;
            }
        }
        // ★ The headline, asserted where it is made (☑ measured 0 over 1056 cells first).
        assert_eq!(
            d2b.emit_mismatch, 0,
            "the cell reader disagrees with today's emitted lateral faces: cyl {k}"
        );
        d2b.ref_band_merges = ref_band_merges;
        d2b.ref_merge_over_boundary = ref_merge_over_boundary;
        // ★★★★★ **The chart's boundary rule is today's, asserted where it is made.** Every line a
        // band of the reference road runs across is one the chart lets a band-shaped pair cross,
        // and the chart makes exactly the merges the reference made — so the cutover's raw
        // emission is today's, band for band (☑ measured 47 = 47 over the suite first).
        assert_eq!(
            ref_merge_over_boundary, 0,
            "a reference band runs across a line the chart calls a boundary: cyl {k}"
        );
        assert_eq!(
            d2b.z_merge_bandlike, ref_band_merges,
            "the chart would merge more or fewer band-shaped intervals than the reference road \
             did: cyl {k}"
        );
        d2b.order_descents = first_claim.windows(2).filter(|w| w[1] < w[0]).count();
        // ★ Asserted where it is made: today's emission order is the chart's cell order, so the
        // cutover renumbers no handle (☑ measured 0 first).
        assert_eq!(
            d2b.order_descents, 0,
            "today's lateral faces are not emitted in the chart's cell order: cyl {k}"
        );
        d2b.band_over_rulings = band_over_rulings;
        probe::d2b::push(d2b);
    }
}

/// The census's ledger — the same shape as `ruling_probe`/`panel_probe`: filled where the fact is
/// made, read by one test that reports it.
#[cfg(test)]
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

    /// **The fourth rung's ledger (D2b-0)**: what the cell reader read, and how it compares to
    /// what the hand-written roads emitted. One row per chart with cells, `no_faces` included —
    /// the reads need no emitted face, only the comparison does.
    ///
    /// ★ Compared **within a row** (`cells`, `band_over_rulings` are copied in), never row-by-row
    /// against another ledger: tests run in parallel and the ledgers interleave independently.
    pub(crate) mod d2b {
        use std::sync::Mutex;

        #[derive(Clone, Copy, Debug, Default)]
        pub(crate) struct Row {
            pub(crate) cells: usize,
            pub(crate) no_faces: bool,
            /// Copied from the D1b row of the same chart, for the merge check below.
            pub(crate) band_over_rulings: usize,
            /// Rulings whose `end` had to be swapped alongside `z` — predicted 0.
            pub(crate) end_swapped: usize,
            /// Every cell's two ends, by what they said.
            pub(crate) end_disk: usize,
            pub(crate) end_exact: usize,
            pub(crate) end_other: usize,
            pub(crate) end_nocircle: usize,
            /// Present cells with an `Other` end — read from their other end alone. The
            /// population a θ-placement finer than [`Chart::arc_around`] would serve.
            pub(crate) other_present: usize,
            /// Cells whose two speaking ends disagreed — `chamber`'s two-end refusal, on the chart.
            pub(crate) src2_disagree: usize,
            /// Cells with a face and no speaking end. Asserted 0 at the record; here as a count
            /// so the reporting test can say the assertion was live.
            pub(crate) src0_present: usize,
            pub(crate) exist_disagree: usize,
            /// The reader refused a cell of this chart by name (`face_spans`' or `panel_faces`'
            /// refusals, or an arrangement-table inconsistency) — the whole class is refused.
            pub(crate) read_refused: usize,
            /// Both-cut cells the trace says the face is not in — the panel road's dropped sectors.
            pub(crate) exist_marks_false: usize,
            /// Cells the reader would emit, cells it could not decide, and — the headline — cells
            /// where it disagrees with what was emitted today.
            pub(crate) emit: usize,
            pub(crate) emit_unknown: usize,
            pub(crate) emit_mismatch: usize,
            /// D1b's 44, judged: lines inside a band whose label changes across them, and lines
            /// inside a band with no circle at all.
            pub(crate) split_flip: usize,
            pub(crate) split_nocircle: usize,
            /// Emitted neighbours with one chamber, across a ruling and across a z-line.
            pub(crate) theta_merge_pairs: usize,
            pub(crate) z_merge_pairs: usize,
            /// Disk-ended intervals whose emitted sectors do not close into one whole circle.
            pub(crate) partial_theta_in_disk_interval: usize,
            /// D2b-1: what the cutover would emit. Intervals that are one chamber round (`Band`),
            /// partial runs (panel rings), runs whose boundary ruling has no node on a cut rim,
            /// emitted whole-circle cells with **both** rims cut (`band_loop` refuses those), and
            /// descents in today's emission order against the chart's cell order.
            pub(crate) full_runs: usize,
            pub(crate) partial_runs: usize,
            pub(crate) run_boundary_no_node: usize,
            pub(crate) whole_both_cut: usize,
            pub(crate) order_descents: usize,
            /// Emitted whole-circle cells (`Band`s over intervals with no ruling), and runs the
            /// panel road would cut further because an interior ruling has a node on a cut rim.
            pub(crate) whole_emitted: usize,
            pub(crate) run_split_at_node: usize,
            /// Band-shaped intervals with both rims cut (today's road sends those to the panel
            /// road); intervals with more than one partial run; and partial-run pairs that would
            /// have to merge across a non-boundary z-line.
            pub(crate) full_run_both_cut: usize,
            pub(crate) intervals_multi_run: usize,
            pub(crate) run_run_nonboundary: usize,
            /// Consecutive band-shaped intervals with one chamber across a line `bands_of` never
            /// cut — the merges the cutover makes at the chart.
            pub(crate) z_merge_bandlike: usize,
            /// The merges the reference road made (interior lines of its bands), and how many of
            /// those lines the chart's boundary rule would refuse to cross.
            pub(crate) ref_band_merges: usize,
            pub(crate) ref_merge_over_boundary: usize,
            /// The shadow comparison: faces compared equal to the reference road's; charts the
            /// chart road emits where the reference refused; charts both refused.
            pub(crate) compared: usize,
            pub(crate) opened: usize,
            pub(crate) both_refused: usize,
        }

        pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

        pub(crate) fn push(r: Row) {
            ROWS.lock()
                .expect("the probe's lock is never held across a panic")
                .push(r);
        }
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

    /// **The cells read their chamber off the horizontal lines, and the reader agrees with today's
    /// roads cell by cell** (capability D, D2b-0).
    ///
    /// The absolute claims are asserted in `census`, where the facts are made: a present cell
    /// always has a speaking end (`src0_present == 0`, the reason no wall-side sign is needed),
    /// and the ends' bookkeeping is total. What this holds, universally over every recorded row
    /// so no interleaving can break it, is the rest of the forecast:
    ///
    /// * `emit_mismatch == 0` — the headline: no cell where the reader and the emitted faces
    ///   disagree. The cutover would change **no** face of today's corpus.
    /// * `src2_disagree == 0` — a cell's two ends never contradict (the two-end agreement
    ///   `bands::chamber` demands, on every cell).
    /// * `exist_disagree == 0`, `read_refused == 0` — the trace and the span tell the same
    ///   existence story wherever both speak.
    /// * `split_flip == 0`, `split_nocircle == 0` — D1b's 44 extra lines are harmless: each has a
    ///   circle and nothing changes across it.
    /// * `partial_theta_in_disk_interval == 0` — a disk-ended interval's kept sectors always close
    ///   into a whole circle, so the cutover can emit it as today's `Band` verbatim.
    /// * `theta_merge_pairs >= band_over_rulings` — every band that spans a wall is, on the chart,
    ///   emitted sectors with one chamber (the merge the cutover must do at the chart).
    ///
    /// And the counters are not vacuous: the fixtures below put every end kind but `Other` on the
    /// ledger, the chained ones drop sectors for existence (**≥ 4**, the panel road's own count,
    /// `tests::a_chained_cylinder_bounded_by_the_first_builds`), and both merge directions are
    /// seen. `end_other` is **unmeasured** if it stays 0 — the θ-placement it would need is not
    /// built, and a zero here is a population claim only once a fixture reaches it.
    #[test]
    fn the_cells_read_their_chamber_from_the_horizontal_lines() {
        use super::probe::d2b::ROWS as D2B;
        let snapshot = || -> Vec<super::probe::d2b::Row> {
            D2B.lock()
                .expect("the probe's lock is never held across a panic")
                .clone()
        };
        let sum = |rows: &[super::probe::d2b::Row], f: fn(&super::probe::d2b::Row) -> usize| {
            rows.iter().map(f).sum::<usize>()
        };
        let before = snapshot();
        // A wall boss (exact arcs, θ- and z-merges), a corner boss and a boss on top (a line with
        // no circle of this cylinder, cells outside the face), and the chained pairs (sectors
        // dropped for existence).
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        for (base, h) in [
            ([2.0, 0.0, -1.0], 4.0),
            ([4.0, 4.0, -1.0], 4.0),
            ([2.0, 2.0, 2.0], 1.0),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(Point3::from_array(base), up, 0.5, h);
            m.rebuild_adjacency();
            crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        }
        let wall_boss = ([2.0, 0.0, -1.0], 4.0, BoolKind::Fuse);
        for second in [
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            // A short boss whose cap (z = 2.5) is a ⊥ line strictly inside the wall boss's span:
            // the wall boss's circle is traced there but bounds no face, so the line is one
            // `bands_of` never cuts — the band-shaped merge population.
            ([6.0, 2.0, 2.0], 0.5, BoolKind::Fuse),
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Fuse),
            ([6.0, 0.0, -1.0], 4.0, BoolKind::Cut),
        ] {
            let (_, r) = crate::tests::chained(wall_boss, second);
            r.expect("the chained pair builds");
        }
        let rows = snapshot();
        assert!(
            rows.len() > before.len(),
            "a boolean recorded no cell reads"
        );

        for r in &rows {
            assert_eq!(
                r.emit_mismatch, 0,
                "the reader disagrees with today's road: {r:?}"
            );
            assert_eq!(
                r.emit_unknown, 0,
                "a cell's chamber could not be read: {r:?}"
            );
            assert_eq!(r.src2_disagree, 0, "a cell's two ends disagree: {r:?}");
            assert_eq!(
                r.src0_present, 0,
                "the record-site assertion was not live: {r:?}"
            );
            assert_eq!(
                r.exist_disagree, 0,
                "trace and span disagree on existence: {r:?}"
            );
            assert_eq!(
                r.read_refused, 0,
                "the reader refused a cell by name: {r:?}"
            );
            assert!(r.opened <= 1, "opened is per chart: {r:?}");
            assert_eq!(r.both_refused, 0, "both roads refused a chart: {r:?}");
            assert_eq!(
                r.split_flip, 0,
                "a line inside a band changes the chamber: {r:?}"
            );
            assert_eq!(
                r.split_nocircle, 0,
                "a line inside a band has no circle: {r:?}"
            );
            assert_eq!(
                r.end_swapped, 0,
                "a ruling arrived with z descending: {r:?}"
            );
            assert_eq!(
                r.partial_theta_in_disk_interval, 0,
                "a disk-ended interval's sectors do not close into a circle: {r:?}"
            );
            assert_eq!(
                r.end_disk + r.end_exact + r.end_other + r.end_nocircle,
                2 * r.cells,
                "two ends per cell: {r:?}"
            );
            assert!(r.emit <= r.cells, "more emitted cells than cells: {r:?}");
            assert!(
                r.other_present <= r.cells,
                "more other-ended cells than cells: {r:?}"
            );
            assert!(
                r.exist_marks_false <= r.cells,
                "more dropped cells than cells: {r:?}"
            );
            if !r.no_faces {
                assert!(
                    r.theta_merge_pairs >= r.band_over_rulings,
                    "a band over a wall is not same-chamber sectors on the chart: {r:?}"
                );
            }
        }
        // ── D2b-1: what the cutover would emit, held as counts (their production twins are
        // honest rejects, so they are not asserted where the fact is made).
        for r in &rows {
            assert_eq!(
                r.run_boundary_no_node, 0,
                "a run's boundary ruling has no rim node: {r:?}"
            );
            assert_eq!(
                r.whole_both_cut, 0,
                "an emitted whole circle has both rims cut: {r:?}"
            );
            assert_eq!(
                r.run_split_at_node, 0,
                "a run passes a ruling with a rim node: {r:?}"
            );
            assert_eq!(
                r.full_run_both_cut, 0,
                "a band-shaped interval has both rims cut: {r:?}"
            );
            assert_eq!(
                r.run_run_nonboundary, 0,
                "a panel would have to stretch across z: {r:?}"
            );
        }
        // ★★★★ Non-vacuity — a zero above must mean "the population is empty", never "nothing
        // was looked at". `end_other` and `intervals_multi_run` are deliberately not in this list.
        assert!(sum(&rows, |r| r.end_disk) > 0, "no disk end was ever read");
        assert!(
            sum(&rows, |r| r.end_exact) > 0,
            "no exact arc end was ever read"
        );
        assert!(
            sum(&rows, |r| r.end_nocircle) > 0,
            "no line without a circle was ever seen"
        );
        assert!(
            sum(&rows, |r| r.emit) > 0,
            "the reader never emitted a cell"
        );
        assert!(
            sum(&rows, |r| r.theta_merge_pairs) > 0,
            "no θ-merge was ever seen"
        );
        assert!(
            sum(&rows, |r| r.z_merge_pairs) > 0,
            "no z-merge was ever seen"
        );
        assert!(
            sum(&rows, |r| r.full_runs) > 0,
            "no band-shaped interval with rulings"
        );
        assert!(
            sum(&rows, |r| r.partial_runs) > 0,
            "no partial run was ever seen"
        );
        assert!(
            sum(&rows, |r| r.whole_emitted) > 0,
            "no whole-circle cell was ever emitted"
        );
        assert!(
            sum(&rows, |r| r.z_merge_bandlike) > 0,
            "no band-shaped merge was ever seen"
        );
        assert!(
            sum(&rows, |r| r.compared) > 0,
            "the shadow comparison never ran"
        );
        // The chained pairs above drop exactly one sector per operation for existence (four),
        // and other tests may add theirs in between — so what is held is the growth.
        assert!(
            sum(&rows, |r| r.exist_marks_false) >= sum(&before, |r| r.exist_marks_false) + 4,
            "the chained pairs dropped fewer than four sectors for existence"
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
