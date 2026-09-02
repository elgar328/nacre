//! **A cylinder's own chart, as an arrangement's line set** (capability D, first rung).
//!
//! The plane side long ago stopped walking faces by hand: each plane class gets a **cell complex**
//! (`arrangement`'s `walk_cells → nest_cells → label_cells → emit_faces`) and the case-work went
//! with it. The cylinder side did not, when this module began: it had three hand-written walks —
//! `bands.rs`' band and panel roads, and `boolean.rs`' `band_loop` slit-leg and
//! `merge_curved_group` — and four cells in a row were each a falsehood inside one of them. The
//! band and panel roads are gone (D2b/D3: this module emits), and the merge reads one rule (D4:
//! a cycle's winding, `seam_step`); the slit-leg remains, generalized to a chain rim.
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
//! The two axes (D1a), the **cells** they cut (D1b), the vertical lines' answers (D2a), the
//! **cell reader** ([`Chart::read_cell`], D2b-0), which reads a cell's chamber and existence off
//! the *horizontal* lines alone (`DiskLabels`/`ArcLabels` through `bands::read_bits`/`face_spans`),
//! and — since the D2b cutover — the **emitter** ([`emit_lateral`]): production's lateral faces
//! come from here, in the band road's own vocabulary and by its own four rules (see the fn doc).
//! The band road it restates (`bands::{band_faces, bands_of, chamber, panel_faces}`) was the
//! **reference** the census zipped this emission against, face by face, through the cutover; D3
//! deleted it, and the census now holds the chart against its own rules ([`census`]).

use crate::arrangement::{ArcLabel, Curved, Label, RulingExtent};
use crate::boolean::{Bound, LocalFace, Rim};
use crate::planes::{ClassIx, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, combinatorics, reject};
use nacre_scalar::Rat;

/// One **horizontal** line: a ⊥ plane class, and where it crosses the axis.
///
/// ★ `t` is exact (`Rat`) — a ⊥ class's axis parameter is a rational, and `bands::axis_param` is
/// the one spelling of it. A class whose parameter is too wide for `Rat` refuses rather than being
/// skipped, as the band road did: a missing band boundary merges two regions whose membership
/// differs, which is a closed and wrong answer.
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
    /// `DiskLabels`/`ArcLabels`; this is the half they never had, and [`Chart::read_cell`] reads
    /// it where they are silent.
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
/// has been bitten before. `name_on` builds this key from the ruling's own pair, and it was the
/// key the band road's rings were joined on while that road was the reference.
type RulingName = (usize, i8);

/// A partial run as the census records it: its first and last ruling, and its chamber.
#[cfg(test)]
type Run = (RulingName, RulingName, Option<(bool, bool)>);

/// **Which of a wall's two rulings a branch node lies on** — `arrangement::node_ruling_side`, the
/// one spelling (`ruling_side` asked of a name) the assembly's `Wall::Ruling { side }` reads.
/// ★ This used to be a second copy of that function's body (D5, 1a): the same `branch_meet` +
/// `ruling_side`, spelled twice one module apart.
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
    // A ruling piece's node names this chart's own cylinder; a name of another cylinder has no
    // side here (M6b's pair, which no piece carries today).
    let (_, cyl, _) = combinatorics::branch_name(n)?;
    if cyl != k {
        return None;
    }
    let w = combinatorics::class_coeffs_rat(jd, wall)?;
    crate::arrangement::node_ruling_side(jd, def, &w, n)
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
    /// that was the key the band road's `panel_faces` joined on — but then the ruling's *label* (D2a) could not be
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

    /// **A station's canonical name on a z-line** — `crossing_on_ruling(c, wall, k, side)`: the
    /// very name `ruling_sweep` gives a piece ending on that line and `split_circles` gives a rim
    /// node there (`NodeId::branch` folds the pair and the root), so «does this rim carry this
    /// station» is **name equality**, and a station the rim does not carry joins the θ order under
    /// the name a piece ending there would have worn.
    ///
    /// ★★★★★ **Not a node of another line** (D5, 1a). The chart used to place such a station by
    /// its piece's `end[0]` — a node on a *different* z-line at the same θ — and where the rim did
    /// carry the station (a tool wall crossing the circle there), `circular_order` was handed one
    /// point under two names and refused it as the coincidence it checks for. Every one of those
    /// ends read `End::Other` (dev-log's Q2 «한 점 두 이름»; ledger `end_other` 321), and the
    /// crossing census's whole `RulingBoundNotYet` column read its rims that way.
    fn station_on(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        c: usize,
        i: usize,
    ) -> Option<combinatorics::NodeId> {
        let (wall, side) = self.ruling_name(jd, k, def, i)?;
        crate::arrangement::crossing_on_ruling(jd, def, c, wall, k, side).ok()
    }
}

/// **The plane classes where this cylinder's rim is a boundary the arrangement already made**:
/// (i) the circles a plane face emitted as a `Bound::Circle`, and (ii) the cut rims, which emit
/// arcs instead and so are invisible to (i). The band road's `bands_of` read both for the same
/// reason; the chart's line set, its boundary rule and its census all read this one list.
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

/// **The lines a lateral face's band may not run across** — the band road's three boundary
/// sources (`bands_of`, deleted in D3), as axis parameters: the rim classes ([`rim_classes`]) and
/// every row's own span ends. A band-shaped emission may continue across any other line, so
/// [`emit_lateral`] merges there and nowhere else, and the census asserts the chart's z-lines ⊇
/// this list where it is made.
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
/// ★★★★★ **The sources were `bands_of`'s (the band road, deleted in D3), and one thing it did is
/// deliberately *not* done here.** That function ended with `keyed.retain(|(t, _)| t ∈ row.span)`:
/// it clipped to the **face** it was asked about, because its row was a face and not a class. A
/// class has no such span — one lateral
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
    // (iii) every ⊥ class, full stop. The band road only took the ones that land on a face's rim
    // because it was answering about a face; a chart covers the whole axis.
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
    // at one `t` are the *same plane*, which is why span ends match against `t` without a
    // tolerance (`boundary_lines`). So one of the two indices survives arbitrarily — which is
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
/// exactly as the band road's `chamber` read a disk and its `panel_faces` read an arc (both deleted
/// in D3; [`Chart::read_cell`] is their one successor).
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
    /// (`end_other`); a sector that merely spans several arcs is no longer here — it is an
    /// [`End::Exact`] run.
    Other,
    /// No circle of this cylinder on this line at all: the line is a ⊥ class outside every
    /// lateral face's span (`circle_on_class` leaves a circle on every class *within* a span).
    NoCircle,
}

impl End<'_> {
    /// **What this end says about the cell's chamber**, or `None` when it says nothing — which
    /// now includes a run of arcs that do not agree. The read is `bands::read_bits`, the one
    /// spelling, and it lives here so the two consumers below cannot drift into two.
    ///
    /// ☑ **Measured: the run case here is not exercised by the suite.** Making this answer `None`
    /// whenever the run holds two or more arcs leaves all 340 lib tests green — because the cells
    /// a run decides today are decided **absent** by the existence read below, and an absent cell
    /// is never asked for its chamber. This stays the general rule rather than a written-out
    /// refusal because it *is* the whole-circle arm's rule (see [`End::Disk`]'s site) with the
    /// sector's own arcs in place of the circle's; narrowing it would be the second spelling.
    fn chamber(&self, side: crate::planes::SolidSide, above: bool) -> Option<(bool, bool)> {
        let one = |l: &Label| crate::bands::read_bits(l, side, above);
        match self {
            End::Disk(l) => Some(one(l)),
            End::Exact(arcs) => {
                let mut it = arcs.iter().map(|a| one(&a.label));
                let first = it.next()?;
                it.all(|b| b == first).then_some(first)
            }
            End::Other | End::NoCircle => None,
        }
    }
}

/// **Why an end reads `Other`** — one name per `None` site of [`Chart::arc_around`] and per
/// `Other` arm of the whole-circle read, so the ledger's `end_other` is a table of causes rather
/// than one number (D5, 1a). The enum lives outside the instrument because the sites that name a
/// cause are production code; the counting is `cfg(test)` (`probe::other`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
// The whole-circle arms record only under `cfg(test)`, so two variants are built nowhere in a
// release build — stated per build rather than blanket-allowed, the `RulingExtent` precedent.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum OtherWhy {
    /// Both walls one ruling: the sector is the whole circle less that ruling.
    SingleCut,
    /// A wall's station could not be named or placed in the order.
    Unplaced,
    /// `circular_order` refused the rim nodes plus the stations.
    OrderFailed,
    /// A run arc between two adjacent rim nodes is not among the split's arcs.
    RowMissing,
    /// The sector's two ends land on one rim node — no arc between them.
    EmptyRun,
    /// A whole-circle cell over a cut rim whose arcs disagree for this side.
    WholeDisagree,
    /// A whole-circle cell over a cut rim with no arcs at all.
    WholeNoArcs,
}

/// `arc_around`'s `None`, with its reason counted: every `End::Other` the suite still produces
/// is named at the site that produced it, so «what is left» is a table and not a guess.
fn other_none<'a>(why: OtherWhy) -> Option<Vec<&'a ArcLabel>> {
    #[cfg(test)]
    probe::other::record(why);
    #[cfg(not(test))]
    let _ = why;
    None
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
        rim: &crate::arrangement::CutRim,
        arcs: &'a [ArcLabel],
    ) -> Option<Vec<&'a ArcLabel>> {
        use OtherWhy as Why;
        let m = rim.nodes.len();
        if x == y && m >= 2 {
            return other_none(Why::SingleCut);
        }
        let mut list: Vec<combinatorics::NodeId> = rim.nodes.clone();
        let mut index_of = |i: usize| -> Option<usize> {
            // A station the rim's own arc decomposition does not hold — the wall's face stops
            // short of the rim — joins the order as itself, under its canonical name, and the
            // search below asks which arc contains it. ★ It used to join under a node of
            // *another* line ([`Chart::station_on`]'s note), which the order refused wherever
            // the rim did hold the station.
            let n = match self.node_on(i, t) {
                Some(n) => n,
                None => self.station_on(jd, k, def, c, i)?,
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
            return other_none(Why::Unplaced);
        };
        let iy = if x == y {
            ix
        } else {
            match index_of(y) {
                Some(i) => i,
                None => return other_none(Why::Unplaced),
            }
        };
        let Ok((order, _)) = crate::arrangement::circular_order(jd, k, def, &list) else {
            return other_none(Why::OrderFailed);
        };
        let n = order.len();
        let pos = |li: usize| order.iter().position(|&o| o == li);
        let (Some(px), Some(py)) = (pos(ix), pos(iy)) else {
            return other_none(Why::Unplaced);
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
            let Some(arc) = arcs.iter().find(|arc| arc.ends == [na, nb]) else {
                return other_none(Why::RowMissing);
            };
            run.push(arc);
            cur = nxt;
        }
        // `a == b` means the sector's ends land on one rim node: no arc separates them, and the
        // caller has nothing to read here.
        if run.is_empty() {
            return other_none(Why::EmptyRun);
        }
        Some(run)
    }

    /// **Read one cell off the horizontal lines** — the function the cutover will call, measured
    /// first against what the hand-written roads emit today.
    ///
    /// No sign is derived here: `band_is_above` is `planes::plus_t_is_above` at the low end and
    /// its negation at the high end — the one spelling the band road's `chamber` and `panel_faces`
    /// both used — and the bits come out through `bands::read_bits`. The ruling labels (the
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
                        #[cfg(test)]
                        let mut seen: Vec<(bool, bool)> = Vec::new();
                        for a in arcs {
                            let b = crate::bands::read_bits(&a.label, side, above[e]);
                            same &= bits.replace(b).is_none_or(|p| p == b);
                            #[cfg(test)]
                            seen.push(b);
                        }
                        match (same, arcs.first()) {
                            (true, Some(a)) => End::Disk(a.label),
                            (false, _) => {
                                #[cfg(test)]
                                probe::other::record_whole(k, t[e], e, above[e], seen);
                                End::Other
                            }
                            (true, None) => {
                                #[cfg(test)]
                                probe::other::record(OtherWhy::WholeNoArcs);
                                End::Other
                            }
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
                                // "sector labels missing" this name states (the band road's
                                // `arc_at` stated it first).
                                let Some(arc) = arcs.iter().find(|a| a.ends == [nx, ny]) else {
                                    return Err(reject(RejectReason::RulingBoundNotYet));
                                };
                                let arc = vec![arc];
                                End::Exact(arc)
                            }
                            _ => self
                                .arc_around(jd, k, def, c, t[e], x, y, rim, arcs)
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
        // (D2b's `(0, 0)`-corner finding), so making the disagreement refuse would regress a
        // defect that is already named rather than fix it. Cross-checking is what this becomes
        // the day that corner is fixed.
        {
            let vertical = |i: usize, starts_here: bool| -> Option<(bool, bool)> {
                let seg = self.theta.get(i)?;
                let l = seg.label?;
                let (wall, sd) = self.ruling_name(jd, k, def, i)?;
                // The sector leaves `x` counter-clockwise and arrives at `y`, so the two walls
                // are read from opposite sides of their own rulings.
                let above = crate::arrangement::plus_theta_is_above(jd, wall, sd)? == starts_here;
                Some(crate::bands::read_bits(&l, side, above))
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
            // is the same fact from the other side. What this replaced was the `by_span`
            // fallback — an axial span that knows no θ and so claimed every sector present.
            let mut span: Option<bool> = None;
            for arc in arcs {
                // `face_spans` refuses by name (two of this solid's faces disagreeing on one arc), and
                // two cut ends disagreeing with each other is `CylinderFaceUndecided`, as the band
                // road's `panel_faces` refused it:
                // the rims of a hole are band boundaries, so a sector exists over its whole height or
                // not at all, and picking an end to believe is the guess this kernel does not make.
                let v = crate::bands::face_spans(arc, side, above[e])?;
                if span.replace(v).is_some_and(|p| p != v) {
                    return Err(reject(RejectReason::CylinderFaceUndecided));
                }
            }
            let Some(v) = span else { continue };
            if by_marks.replace(v).is_some_and(|p| p != v) {
                return Err(reject(RejectReason::CylinderFaceUndecided));
            }
        }
        // ★★★★★ **Two silent ends and a span that says «present» is not read as present** (D5,
        // 1a). The span knows no θ, so over a cell whose both cut ends could not be paired with
        // their rims it used to claim the face was there — the guess that put a face over a hole
        // in the corner boss × mid slab, caught only because the emitter refused that class one
        // step later. With stations placed by name (`station_on`) no cell in the suite reaches
        // here (`src0_present` 14 → 0, the census's record-site assertion made structural), so
        // this is a guard on the reader's premise and the name is the chart's own for a cell it
        // cannot read.
        if by_marks.is_none() && by_span && ends.iter().all(|e| matches!(e, End::Other)) {
            return Err(reject(RejectReason::CylinderGateUndecided));
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
/// the band road's `panel_faces` emitted, node for node. Four rules, each the band road's, restated
/// on the chart and measured against it before this was written (D2b-1; the road itself was deleted
/// in D3, and the census's own walk now predicts this emitter's face count — `emitted_faces`):
///
/// * **`Band` iff band-shaped and not both rims cut** — the band road's dispatch. A whole-circle
///   cell with both rims cut is still a `Band`; `band_loop` refuses that one by name.
/// * **A band runs across a line only where the band road never cut**: not a rim class (an emitted
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
                        lo: Rim::Circle(chart.z_lines[lo].class),
                        hi: Rim::Circle(chart.z_lines[hi].class),
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
            // The band road's dispatch: a both-cut interval is the panel rings', whatever its
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
                // ★ Stage-0 instrument (D5): the one site every `RulingBoundNotYet` of the
                // crossing census dies at — a run-end ruling with no node on the interval's
                // z-line. Recorded with what that z-line says under this sector, so the plan's
                // diagnosis («the line is transversal there, so no class has a node») is
                // measured rather than argued. Production behaviour is unchanged.
                let node =
                    |r: usize, t: Rat, e: usize| -> Result<combinatorics::NodeId, BoolError> {
                        match chart.node_on(r, t) {
                            Some(n) => Ok(n),
                            None => {
                                #[cfg(test)]
                                probe::rbny::record(
                                    jd, k, def, &chart, &lines, curved, side, &reads[f], r, t, e,
                                );
                                #[cfg(not(test))]
                                let _ = e;
                                Err(ladder())
                            }
                        }
                    };
                let (a_lo, a_hi) = (node(a, t_lo, 0)?, node(a, t_hi, 1)?);
                let (b_lo, b_hi) = (node(b, t_lo, 0)?, node(b, t_hi, 1)?);
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
    emission: &Result<Vec<LocalFace>, BoolError>,
) {
    for (k, _) in cyls.iter().enumerate() {
        // ★ A class whose axis parameter cannot be stated is skipped, not recorded: `emit_lateral`
        // refuses on the very same `chart_of`, so this is a refusal the emitter already named, not
        // an outcome the census introduces. ☑ Measured **0** across the suite either way.
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
        // The boundary lines, for the merge counters; a refusal is a skip, as for the chart.
        let Ok(boundary) = boundary_lines(jd, k, &cyls[k].def, plane_faces, curved, rows) else {
            continue;
        };
        // ★★★★★ **The chart's lines cover its own boundary rule, asserted where the fact is made.**
        // The rim half is shared with `chart_of` by construction; the load-bearing half is the
        // rows' span ends — a *different spelling* of the axis parameter (the tracer's `t_range`)
        // than the classes' `axis_param` — so this says every lateral face's rim lands, with exact
        // `Rat` equality, on a ⊥ class the chart collected. Not a tautology: a rim on a class
        // `chart_of` skipped (no rational coefficients) or a span end the two spellings disagree on
        // would both fail here. (D3: the successor of the ⊇ check against `bands_of`.)
        for t in &boundary {
            assert!(
                chart.z_lines.iter().any(|l| l.t == *t),
                "a chart lost a line its own boundary rule names: cyl {k}, t {t:?}"
            );
        }
        let side = mine.first().map(|r| r.side);
        assert!(
            mine.iter().all(|r| Some(r.side) == side),
            "one cylinder class carries rows of both solids: cyl {k}"
        );

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
                    test: std::thread::current().name().unwrap_or("?").to_string(),
                    cells: cells.len(),
                    read_refused: 1,
                    ..probe::d2b::Row::default()
                });
                continue;
            }
        };
        let mut d2b = probe::d2b::Row {
            test: std::thread::current().name().unwrap_or("?").to_string(),
            emission_reason: match emission {
                Err(BoolError::Rejected { reason, .. }) => Some(*reason),
                _ => None,
            },
            cells: cells.len(),
            emitter_refused: emission.is_err(),
            end_swapped: chart.end_swapped,
            ..probe::d2b::Row::default()
        };
        // D5 stage 0 (P4): the premise of naming a station on a line — every (z-line, station)
        // pair of this chart asked for the station's canonical name there.
        {
            let mut names: Vec<(usize, i8)> = Vec::new();
            for j in 0..chart.theta.len() {
                if let Some(n) = chart.ruling_name(jd, k, def, j) {
                    if !names.contains(&n) {
                        names.push(n);
                    }
                }
            }
            for l in &chart.z_lines {
                for &(wall, sd) in &names {
                    d2b.station_pairs += 1;
                    if crate::arrangement::crossing_on_ruling(jd, def, l.class, wall, k, sd)
                        .is_err()
                    {
                        d2b.station_name_failures += 1;
                    }
                }
            }
        }
        // ── D5, 1b — the shadow: the class's faces as regions of its chart, beside what today's
        // road ships (its emission **after** the curved cleaning pass, which is the face set the
        // assembly actually receives). Compared by signature, never by order.
        {
            let (shadow, rep, refused) = match regions::emit_regions(
                jd, k, def, kind, side, &chart, &lines, &cells, &reads, curved,
            ) {
                Ok((f, r)) => (f, r, false),
                Err(_) => (Vec::new(), regions::Report::default(), true),
            };
            let today: Vec<LocalFace> = match emission {
                Ok(faces) => crate::boolean::unify_curved_faces(faces.clone(), &curved.cut_rims)
                    .into_iter()
                    .filter(|f| matches!(f.surf, ClassIx::Cyl(c) if c == k))
                    .collect(),
                Err(_) => Vec::new(),
            };
            // ★ Two spellings of one rim, read as one: today's emitter says `Rim::Circle(c)` of
            // a **cut** rim (the assembly's `rim_walk` then walks its arcs), the cleaning pass
            // and the shadow say the chain of those arcs. Same edges, same nodes — so a cut
            // circle is signed as the chain it stands for, nodes included.
            let sig = |f: &LocalFace| -> String {
                let rim = |r: &Rim| match r {
                    Rim::Circle(c) => match curved.cut_rims.get(&(k, *c)) {
                        Some(cr) => format!("K{}", cr.nodes.len()),
                        None => format!("C{c}"),
                    },
                    Rim::Chain(r) => format!("K{}", r.nodes.len()),
                };
                let bound = |b: &Bound| -> String {
                    match b {
                        Bound::Ring(r) => format!("ring{}", r.nodes.len()),
                        Bound::Circle { cyl } => format!("circle{cyl}"),
                        Bound::Band { lo, hi } => format!("band[{},{}]", rim(lo), rim(hi)),
                    }
                };
                let mut nodes: Vec<String> = f
                    .poly_rings()
                    .flat_map(|r| r.nodes.iter().map(|n| format!("{n:?}")))
                    .collect();
                if let Bound::Band { lo, hi } = &f.outer {
                    for r in [lo, hi] {
                        if let Rim::Circle(c) = r {
                            if let Some(cr) = curved.cut_rims.get(&(k, *c)) {
                                nodes.extend(cr.nodes.iter().map(|n| format!("{n:?}")));
                            }
                        }
                    }
                }
                nodes.sort();
                nodes.dedup();
                let mut inner: Vec<String> = f.inner.iter().map(bound).collect();
                inner.sort();
                format!(
                    "{} flip{} inner[{}] nodes[{}]",
                    bound(&f.outer),
                    f.flip,
                    inner.join(","),
                    nodes.join(",")
                )
            };
            let mut a: Vec<String> = shadow.iter().map(sig).collect();
            a.sort();
            let mut b: Vec<String> = today.iter().map(sig).collect();
            b.sort();
            if a != b && !emission.is_err() {
                eprintln!(
                    "1b-sig cyl {k}\n  shadow: {}\n  today:  {}",
                    a.join(" | "),
                    b.join(" | ")
                );
            }
            probe::regions::push(probe::regions::Row {
                test: std::thread::current().name().unwrap_or("?").to_string(),
                cyl: k,
                emitter_refused: emission.is_err(),
                shadow_refused: refused,
                shadow_faces: shadow.len(),
                today_faces: today.len(),
                agree: a == b,
                rep,
            });
        }
        for (r, cell) in reads.iter().zip(&cells) {
            for e in &r.ends {
                match e {
                    End::Disk(_) => d2b.end_disk += 1,
                    End::Exact(a) => {
                        d2b.end_exact += 1;
                        d2b.exact_run_arcs += a.len() - 1;
                    }
                    End::Other => {
                        d2b.end_other += 1;
                        // D5, 1a: with stations placed by name the only `Other` left should be
                        // the single-cut sector (both walls one ruling — `arc_around`'s first
                        // return); anything else is counted apart so it can be named.
                        if cell.walls.is_some_and(|[x, y]| x == y) {
                            d2b.end_other_single_cut += 1;
                        }
                    }
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
            // A present cell with a line that carries no circle of this cylinder: the
            // `circle_on_class` premise at one end instead of both (`src0_present`).
            if r.present && r.ends.iter().any(|e| matches!(e, End::NoCircle)) {
                d2b.nocircle_present += 1;
            }
            // ★ The one-mark contract the band road's `panel_probe` used to hold on its rim reads:
            // every cut end read carries exactly one lateral mark of its own solid (`Seated` is a
            // planar face's word — `face_spans` skips it for the same reason). Counted here, not
            // in `read_cell`, so the emitter's read and the census's read are not counted twice.
            for e in &r.ends {
                let End::Exact(arcs) = e else { continue };
                // ★ Each arc of a run is read on its own — the contract is per arc, and a run
                // that spans several is exactly where a violation would first show.
                for arc in arcs {
                    let lateral = arc
                        .marks
                        .iter()
                        .filter(|(s, kd)| {
                            *s == side && !matches!(kd, crate::arrangement::SegKind::Seated { .. })
                        })
                        .count();
                    d2b.arcs_read += 1;
                    match lateral {
                        0 => d2b.arcs_no_mark += 1,
                        1 => {}
                        _ => d2b.arcs_multi_mark += 1,
                    }
                    // ☑ Measured 0/0 over 756 reads (lib) and 828 (corpus) before promotion (D3).
                    // ★ **Refuted once the scan named crossings on rulings (E3):** a ⊥ cap through
                    // a wall boss's notch cuts the circle *inside the lateral's hole*, and the arc
                    // there carries **no** mark — the face is absent along it, and the cell reads
                    // absent. So the contract is «at most one, and none exactly where the cell is
                    // not there»; the population that measured 0/0 had no cut through a hole.
                    assert!(
                        lateral <= 1,
                        "a cut end read carries {lateral} lateral marks of its own solid: {arc:?}"
                    );
                    assert!(
                        lateral == 1 || !r.present,
                        "a cut end read with no lateral mark of its own solid reads present: {arc:?}"
                    );
                }
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
        // the fact is made: `circle_on_class` leaves a circle on every ⊥ class within a lateral
        // face's span (Crosses inside, Grazes at the rims), and `emit_faces` labels every circle
        // cell and arc outside the keep filter — which is why no "which side of the wall" sign is
        // needed to read a cell. ★ **Refuted as a bare zero by the corner boss × a mid slab
        // (E3):** the circle is there, but a cut end the reader cannot pair with its rim
        // (`End::Other`, E2-2's item) leaves *no* end speaking and `present` falls back to the
        // row's span — which cannot see the hole. The true proposition is the `src2_disagree`
        // guard's: over such a cell the emitter builds nothing, it refuses the class by name.
        assert!(
            d2b.src0_present == 0 || emission.is_err(),
            "a cell with a face has no label at either end and the emitter read it: cyl {k}"
        );
        // ★★★★★ **A cell with a face never has an end that says nothing** (D5, 1a) — promoted
        // from a count the suite measured 60 → **0** once stations were placed by name. What
        // `Other` still means is a whole-circle interval *beyond* a face's span whose cut rim's
        // arcs disagree about the far side (the other solid stands on one side of its own wall
        // there), and such a cell is absent. A present cell reading `Other` would be a chart with
        // a θ boundary it has no line for — named here, at the fact, not read from the other end.
        assert_eq!(
            d2b.other_present, 0,
            "a cell with a face has an end that says nothing: cyl {k}"
        );
        // ★ A run whose boundary ruling has no node on a cut rim is never built — it is the
        // `RulingBoundNotYet` population (D5's), refused by the emitter at the run end. Asserted
        // where the fact is made rather than only in the reporting test, whose rows depend on
        // what ran before it.
        assert!(
            d2b.run_boundary_no_node == 0 || emission.is_err(),
            "a run's boundary ruling has no rim node and the emitter built it: cyl {k}"
        );
        // ★ Promoted from counts to record-site assertions once the suite measured them 0
        // (D1b's `unnamed == 0` discipline): a reporting test sees only the rows recorded before
        // it, an assertion here sees every chart and names the offending test.
        assert_eq!(
            d2b.end_swapped, 0,
            "a ruling arrived with z descending: cyl {k}"
        );
        // ★★★★★ **Over a present cell whose chamber could not be read, the emitter emits
        // nothing** — the true proposition. The plain `src2_disagree == 0` was promoted on a
        // lib-suite zero that turned out to be the suite's, not the corpus's: the (0,0)-corner
        // boss (`wal corner-lo`) has plate-class disk cells with no B material while its caps
        // carry it (an arrangement label defect, D2b's finding), so its cells' ends disagree — and
        // the emitter refuses the class by name rather than reading either end. The guard does
        // not hide the defect; it says no face is built over it. ★ Stated on `emit_unknown`
        // (present cells only), not on `src2_disagree`, which also counts absent cells the
        // emitter rightly ignores — that stays a ledger column (D3 self-check).
        assert!(
            d2b.emit_unknown == 0 || emission.is_err(),
            "a present cell's chamber could not be read, yet the emitter emitted: cyl {k}"
        );
        // ★ The premise at one end: a present cell never meets a line with no circle of its own
        // cylinder (`circle_on_class` marks every ⊥ class within a face's span).
        assert_eq!(
            d2b.nocircle_present, 0,
            "a cell with a face meets a line with no circle of its cylinder: cyl {k}"
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
                // …unless both rims are cut: the emitter sends that interval to the panel rings
                // (the band road's dispatch), and a `Band` there would be refused by `band_loop`.
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
        // is one of the boundary lines (`boundary_lines`). The sub-count across a line a
        // plane face emitted as a `Bound::Circle` is asserted 0 where it is made: a whole circle
        // through a present band with the same chamber on both sides is a T-junction, which the
        // arrangement does not produce.
        let mut z_merge_bandlike_circle = 0usize;
        for i in 1..n_int {
            let t = chart.z_lines[i].t;
            // Two band-shaped neighbours whose chambers *differ* across a non-boundary line: a
            // chamber change with no face edge on the lateral there. ☑ Counted, predicted 0.
            if let (Some(p), Some(q)) = (bandlike[i - 1], bandlike[i])
                && p != q
                && !is_boundary(t)
            {
                d2b.z_flip_nonboundary += 1;
            }
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
        let mut unpaired = 0usize;
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
            // ★ **A wall's rulings need not come in pairs any more (E2-2).** A band's lateral
            // reaches both rulings of every through-axis wall, so the walk below crossed each
            // wall twice; a panel's or a chain's boundary may run along one ruling of a wall
            // while the other is no face's edge at all (a quarter boss on a corner). Around such
            // an interval the walk crosses that wall once and cannot return to its start — not a
            // contradiction, a walk with an open end — so the closure is asserted only where every
            // wall is crossed an even number of times, and the rest are counted (`unpaired`).
            let paired = {
                let mut walls: Vec<usize> = alive.iter().map(|t| t.wall).collect();
                walls.sort_unstable();
                walls.chunk_by(|a, b| a == b).all(|c| c.len() % 2 == 0)
            };
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
            if all_known && !paired {
                unpaired += 1;
            }
            if all_known && paired {
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
                // every interval asserted here carries an even count of rulings per wall (the
                // `paired` guard above; `odd_k` counted 0 until E2-2's panels), so a *per-ruling*
                // systematic flip cancels itself and passes. ☑ Probed both ways:
                // adding one flip per ruling stays green, seeding the accumulator wrong goes red.
                // So this holds the labels against each other; what holds their absolute sense is
                // `ruling_probe::SIDE_CHECK`, which compares against content and is not a walk.
                assert!(
                    acc == [false; 2],
                    "the chart's labels do not close around an interval: cyl {k}, interval {i}"
                );
            }
        }

        probe::push(probe::Row {
            z_lines: chart.z_lines.len(),
            theta: chart.theta.len(),
            rows: mine.len(),
            refused: false,
            cells: cells.len(),
            whole_circle,
            odd_k,
        });
        probe::d2::push(probe::d2::Row {
            rulings,
            wall_flips,
            intervals_with_flip,
            closes,
            does_not_close,
            unpaired,
            grazing_rulings,
        });
        // ── The emitter, seen from the census (D3): its face count per class is predicted by the
        // census's own walk of the same chart — band-shaped intervals minus the merges across
        // non-boundary lines, plus the partial runs. The emitter's `pending` logic and this walk
        // are different code, so the equality is not a tautology; it is the eye on the emitter
        // that the reference road's shadow used to be. ★ Premise: no band-shaped interval has
        // both rims cut (`full_run_both_cut`, `whole_both_cut` — asserted 0 in the reporting
        // test); the day one does, the emitter sends it to the panel road as `m` rings and this
        // needs a term.
        if let Ok(faces) = emission {
            d2b.emitted_faces = faces
                .iter()
                .filter(|f| matches!(f.surf, ClassIx::Cyl(c) if c == k))
                .count();
            // ☑ Measured before promotion (D3): 407 = 254 + 113 − 47 + 87 over the lib suite,
            // 590 over the corpus too, class by class.
            assert_eq!(
                d2b.emitted_faces,
                d2b.whole_emitted + d2b.full_runs - d2b.z_merge_bandlike + d2b.partial_runs,
                "the emitter's face count is not the census's: {d2b:?}"
            );
        }
        probe::d2b::push(d2b);
    }
}

/// The census's ledger — the same shape as `ruling_probe`: filled where the fact is made, read by
/// one test that reports it.
/// **D5, 1b — a lateral's result faces as regions of its chart** (shadow, `cfg(test)`: built
/// beside today's emission and compared by the census; production untouched).
///
/// The plane side's stages, on the chart: cells → labels (already read) → emitted cells →
/// **connected components** → the boundary of each component → the boundary's runs cut into the
/// **neighbouring class's pieces** (a rim's arcs, a wall's ruling pieces) → cycles → `Bound`.
/// No band/panel dispatch, no both-rims-cut rule, no run ladder: a band, a panel, a chain rim, a
/// hole and a notched panel are one thing here.
#[cfg(test)]
pub(crate) mod regions {
    use super::{Cell, CellRead, Chart, Lines};
    use crate::arrangement::Curved;
    use crate::boolean::{Bound, LocalFace, Ring, Wall};
    use crate::combinatorics::NodeId;
    use crate::planes::{ClassIx, WorkingPlane};
    use crate::tolerant::Judge;
    use crate::{BoolError, RejectReason, reject};

    /// What the walk found for one class.
    #[derive(Clone, Debug, Default)]
    pub(crate) struct Report {
        pub(crate) components: usize,
        pub(crate) faces: usize,
        /// A boundary run along a rim whose end is not a rim node (the consistency claim's
        /// violation — predicted 0).
        pub(crate) run_end_not_piece_end: usize,
        /// A boundary run along a ruling whose pieces do not tile it (predicted 0).
        pub(crate) ruling_run_not_tiled: usize,
        /// A run along an uncut rim that is not the whole circle.
        pub(crate) rim_run_not_whole: usize,
        /// Corners where a component's boundary passes twice (a hole touching the outer
        /// boundary at a corner) — walked by the left-turn rule, counted.
        pub(crate) pinch_corners: usize,
        /// Cycles whose pieces did not chain end to end, or a component with no outer cycle.
        pub(crate) unchained: usize,
        /// `classify_cycles` abstained.
        pub(crate) classify_abstain: usize,
        pub(crate) ring_faces: usize,
        pub(crate) band_faces: usize,
    }

    /// A corner of the chart's grid: `(z-line index, global station index)`.
    type Corner = (usize, usize);

    /// One directed boundary edge of the grid, with the region on its left.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Edge {
        /// Along z-line `line`, over unit sector `unit` (from station `unit` to `unit + 1`);
        /// `ccw` when the region lies above the line.
        Arc { line: usize, unit: usize, ccw: bool },
        /// Along station `station`, over interval `interval`; `up` when the region lies on the
        /// station's −θ side (the cell whose sector *ends* at the station).
        Ruling {
            interval: usize,
            station: usize,
            up: bool,
        },
    }

    /// One boundary piece in travel order: a neighbouring class's own edge.
    struct Piece {
        wall: Wall,
        ends: (NodeId, NodeId),
    }

    /// The regions of class `k`, as `LocalFace`s, and the report. `Err` only for what the
    /// emitter would refuse by name (a chart whose stations cannot be ordered).
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(crate) fn emit_regions(
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        kind: crate::BoolKind,
        side: crate::planes::SolidSide,
        chart: &Chart,
        lines: &Lines,
        cells: &[Cell],
        reads: &[CellRead<'_>],
        curved: &Curved,
    ) -> Result<(Vec<LocalFace>, Report), BoolError> {
        let mut rep = Report::default();
        let undecided = || reject(RejectReason::WitnessNotRational);

        // ── 1. Stations: one per (wall, side), in one global θ order. ──
        let mut names: Vec<(usize, i8)> = Vec::new();
        let mut reps: Vec<NodeId> = Vec::new();
        let mut name_of_theta: Vec<usize> = Vec::with_capacity(chart.theta.len());
        for j in 0..chart.theta.len() {
            let n = chart.ruling_name(jd, k, def, j).ok_or_else(undecided)?;
            let ni = match names.iter().position(|x| *x == n) {
                Some(p) => p,
                None => {
                    names.push(n);
                    reps.push(chart.theta[j].end[0]);
                    names.len() - 1
                }
            };
            name_of_theta.push(ni);
        }
        let ns = names.len();
        let mut station_of_name = vec![0usize; ns];
        if ns > 0 {
            let (order, _) =
                crate::arrangement::circular_order(jd, k, def, &reps).map_err(|_| undecided())?;
            for (p, &ni) in order.iter().enumerate() {
                station_of_name[ni] = p;
            }
        }
        let name_at_station = |s: usize| -> (usize, i8) {
            let ni = station_of_name
                .iter()
                .position(|&p| p == s)
                .expect("a station index comes from this table");
            names[ni]
        };
        let station_of_theta = |j: usize| station_of_name[name_of_theta[j]];
        // A cell's sector in global stations: `None` is the whole circle.
        let sector: Vec<Option<(usize, usize)>> = cells
            .iter()
            .map(|c| {
                c.walls
                    .map(|[x, y]| (station_of_theta(x), station_of_theta(y)))
            })
            .collect();
        // Does the cell cover unit sector `u` (from station u to u + 1)?
        let covers = |ci: usize, u: usize| -> bool {
            match sector[ci] {
                None => true,
                Some((a, b)) if a == b => true, // one cut: the circle less a point
                Some((a, b)) => (u + ns - a) % ns < (b + ns - a) % ns,
            }
        };
        let n_int = chart.z_lines.len().saturating_sub(1);
        // Cells by interval, in sector order (as `cells()` pushed them).
        let mut by_int: Vec<Vec<usize>> = vec![Vec::new(); n_int];
        for (ci, c) in cells.iter().enumerate() {
            by_int[c.interval].push(ci);
        }
        let units = ns.max(1);

        // ── 2. Emitted cells and their components (4-adjacency, θ-periodic). ──
        let emitted: Vec<bool> = reads.iter().map(|r| r.emit == Some(true)).collect();
        let neighbours = |ci: usize| -> Vec<usize> {
            let c = &cells[ci];
            let mut out = Vec::new();
            let row = &by_int[c.interval];
            if row.len() > 1 {
                let s = c.sector;
                out.push(row[(s + 1) % row.len()]);
                out.push(row[(s + row.len() - 1) % row.len()]);
            }
            for i2 in [c.interval.wrapping_sub(1), c.interval + 1] {
                if i2 >= n_int {
                    continue;
                }
                for &cj in &by_int[i2] {
                    let overlap = ns == 0 || (0..ns).any(|u| covers(ci, u) && covers(cj, u));
                    if overlap {
                        out.push(cj);
                    }
                }
            }
            out
        };
        let mut comp: Vec<Option<usize>> = vec![None; cells.len()];
        let mut components: Vec<Vec<usize>> = Vec::new();
        for ci in 0..cells.len() {
            if !emitted[ci] || comp[ci].is_some() {
                continue;
            }
            let id = components.len();
            let mut stack = vec![ci];
            let mut members = Vec::new();
            comp[ci] = Some(id);
            while let Some(x) = stack.pop() {
                members.push(x);
                for y in neighbours(x) {
                    if emitted[y] && comp[y].is_none() {
                        comp[y] = Some(id);
                        stack.push(y);
                    }
                }
            }
            members.sort_unstable();
            components.push(members);
        }
        rep.components = components.len();

        // Endpoints (corners) of a directed edge, and its direction of travel
        // (0 = +θ, 1 = +z, 2 = −θ, 3 = −z).
        let ends = |e: Edge| -> (Corner, Corner) {
            match e {
                Edge::Arc { line, unit, ccw } => {
                    let (a, b) = ((line, unit), (line, (unit + 1) % units));
                    if ccw { (a, b) } else { (b, a) }
                }
                Edge::Ruling {
                    interval,
                    station,
                    up,
                } => {
                    let (a, b) = ((interval, station), (interval + 1, station));
                    if up { (a, b) } else { (b, a) }
                }
            }
        };
        let dir = |e: Edge| -> usize {
            match e {
                Edge::Arc { ccw: true, .. } => 0,
                Edge::Arc { ccw: false, .. } => 2,
                Edge::Ruling { up: true, .. } => 1,
                Edge::Ruling { up: false, .. } => 3,
            }
        };
        let same_run = |a: Edge, b: Edge| -> bool {
            match (a, b) {
                (
                    Edge::Arc {
                        line: l1, ccw: c1, ..
                    },
                    Edge::Arc {
                        line: l2, ccw: c2, ..
                    },
                ) => l1 == l2 && c1 == c2,
                (
                    Edge::Ruling {
                        station: s1,
                        up: u1,
                        ..
                    },
                    Edge::Ruling {
                        station: s2,
                        up: u2,
                        ..
                    },
                ) => s1 == s2 && u1 == u2,
                _ => false,
            }
        };
        // The class of a z-line whose circle this cylinder marks — a cut rim first.
        let class_on = |line: usize| -> Option<usize> {
            let t = chart.z_lines[line].t;
            let cs = lines.classes(t);
            cs.iter()
                .copied()
                .find(|&c| curved.cut_rims.contains_key(&(k, c)))
                .or_else(|| {
                    cs.iter()
                        .copied()
                        .find(|&c| curved.disk_labels.contains_key(&(k, c)))
                })
        };
        let station_name = |c: usize, s: usize| -> Option<NodeId> {
            let (wall, sd) = name_at_station(s);
            crate::arrangement::crossing_on_ruling(jd, def, c, wall, k, sd).ok()
        };

        // ── 3–7. Per component: boundary edges → cycles → runs → pieces → rings → bound. ──
        let mut faces: Vec<LocalFace> = Vec::new();
        for (id, members) in components.iter().enumerate() {
            let in_c = |ci: usize| comp[ci] == Some(id);
            let cell_at = |i: usize, u: usize| -> Option<usize> {
                by_int[i].iter().copied().find(|&cj| covers(cj, u))
            };
            let mut edges: Vec<Edge> = Vec::new();
            for &ci in members {
                let c = &cells[ci];
                for u in 0..units {
                    if ns > 0 && !covers(ci, u) {
                        continue;
                    }
                    let below = if c.interval == 0 {
                        None
                    } else {
                        cell_at(c.interval - 1, u)
                    };
                    if !below.is_some_and(in_c) {
                        edges.push(Edge::Arc {
                            line: c.interval,
                            unit: u,
                            ccw: true,
                        });
                    }
                    let above = if c.interval + 1 >= n_int {
                        None
                    } else {
                        cell_at(c.interval + 1, u)
                    };
                    if !above.is_some_and(in_c) {
                        edges.push(Edge::Arc {
                            line: c.interval + 1,
                            unit: u,
                            ccw: false,
                        });
                    }
                }
                if let Some((a, b)) = sector[ci] {
                    if a != b {
                        let row = &by_int[c.interval];
                        let s = c.sector;
                        let right = row[(s + 1) % row.len()];
                        let left = row[(s + row.len() - 1) % row.len()];
                        if !in_c(right) {
                            edges.push(Edge::Ruling {
                                interval: c.interval,
                                station: b,
                                up: true,
                            });
                        }
                        if !in_c(left) {
                            edges.push(Edge::Ruling {
                                interval: c.interval,
                                station: a,
                                up: false,
                            });
                        }
                    }
                }
            }
            edges.sort_by_key(|e| match *e {
                Edge::Arc { line, unit, ccw } => (0, line, unit, usize::from(ccw)),
                Edge::Ruling {
                    interval,
                    station,
                    up,
                } => (1, interval, station, usize::from(up)),
            });
            edges.dedup();

            // ── cycles, by the left-turn rule at a corner with two ways out ──
            let mut used = vec![false; edges.len()];
            let mut cycles: Vec<Vec<Edge>> = Vec::new();
            for start in 0..edges.len() {
                if used[start] {
                    continue;
                }
                let mut cyc = vec![edges[start]];
                used[start] = true;
                let (first, mut at) = ends(edges[start]);
                let mut prev = dir(edges[start]);
                while at != first {
                    let mut cands: Vec<(usize, usize)> = (0..edges.len())
                        .filter(|&j| !used[j] && ends(edges[j]).0 == at)
                        .map(|j| {
                            let rel = (dir(edges[j]) + 4 - prev) % 4;
                            // left turn first, then straight, then right, then back
                            let rank = match rel {
                                1 => 0,
                                0 => 1,
                                3 => 2,
                                _ => 3,
                            };
                            (rank, j)
                        })
                        .collect();
                    if cands.len() > 1 {
                        rep.pinch_corners += 1;
                    }
                    cands.sort_unstable();
                    let Some(&(_, j)) = cands.first() else {
                        rep.unchained += 1;
                        break;
                    };
                    used[j] = true;
                    cyc.push(edges[j]);
                    prev = dir(edges[j]);
                    at = ends(edges[j]).1;
                }
                cycles.push(cyc);
            }

            // ── runs → pieces → rings ──
            let mut rings: Vec<Ring> = Vec::new();
            let mut rims_lo: Vec<usize> = Vec::new();
            let mut rims_hi: Vec<usize> = Vec::new();
            let mut outer: Option<usize> = None;
            let lowest = cells[members[0]].interval;
            let mut ok = true;
            for cyc in &cycles {
                let mut runs: Vec<Vec<Edge>> = Vec::new();
                for &e in cyc {
                    match runs.last_mut() {
                        Some(r) if same_run(r[0], e) => r.push(e),
                        _ => runs.push(vec![e]),
                    }
                }
                if runs.len() > 1 && same_run(runs[0][0], runs[runs.len() - 1][0]) {
                    let last = runs.pop().expect("two runs at least");
                    let first = std::mem::take(&mut runs[0]);
                    runs[0] = last.into_iter().chain(first).collect();
                }
                let mut pieces: Vec<Piece> = Vec::new();
                let mut whole_rim: Option<(usize, bool)> = None;
                for run in &runs {
                    let (ca, _) = ends(run[0]);
                    let (_, cb) = ends(*run.last().expect("a run has an edge"));
                    match run[0] {
                        Edge::Arc { line, ccw, .. } => {
                            let Some(c) = class_on(line) else {
                                rep.run_end_not_piece_end += 1;
                                ok = false;
                                continue;
                            };
                            let full = ns == 0 || run.len() == ns;
                            let Some(rim) = curved.cut_rims.get(&(k, c)) else {
                                if !full {
                                    rep.rim_run_not_whole += 1;
                                    ok = false;
                                    continue;
                                }
                                whole_rim = Some((c, ccw));
                                continue;
                            };
                            let arcs = &curved.arc_labels[&(k, c)];
                            let m = rim.nodes.len();
                            let (pa, pb) = if full {
                                (0, 0)
                            } else {
                                let (Some(na), Some(nb)) =
                                    (station_name(c, ca.1), station_name(c, cb.1))
                                else {
                                    rep.run_end_not_piece_end += 1;
                                    ok = false;
                                    continue;
                                };
                                let (Some(pa), Some(pb)) = (
                                    rim.nodes.iter().position(|&n| n == na),
                                    rim.nodes.iter().position(|&n| n == nb),
                                ) else {
                                    rep.run_end_not_piece_end += 1;
                                    ok = false;
                                    continue;
                                };
                                (pa, pb)
                            };
                            if m < 2 {
                                rep.run_end_not_piece_end += 1;
                                ok = false;
                                continue;
                            }
                            let mut p = pa;
                            loop {
                                let q = if ccw { (p + 1) % m } else { (p + m - 1) % m };
                                let (from, to) = if ccw {
                                    (rim.nodes[p], rim.nodes[q])
                                } else {
                                    (rim.nodes[q], rim.nodes[p])
                                };
                                if !arcs.iter().any(|a| a.ends == [from, to]) {
                                    rep.run_end_not_piece_end += 1;
                                    ok = false;
                                }
                                pieces.push(Piece {
                                    wall: Wall::Arc { cyl: k, ccw },
                                    ends: if ccw { (from, to) } else { (to, from) },
                                });
                                p = q;
                                if p == pb {
                                    break;
                                }
                            }
                        }
                        Edge::Ruling { station, up, .. } => {
                            let (i0, i1) = if up { (ca.0, cb.0) } else { (cb.0, ca.0) };
                            let (t_lo, t_hi) = (chart.z_lines[i0].t, chart.z_lines[i1].t);
                            let (wall, sd) = name_at_station(station);
                            let ni = names
                                .iter()
                                .position(|&n| n == (wall, sd))
                                .expect("a station has a name");
                            let mut segs: Vec<usize> = (0..chart.theta.len())
                                .filter(|&j| {
                                    name_of_theta[j] == ni
                                        && chart.theta[j].z[0] >= t_lo
                                        && chart.theta[j].z[1] <= t_hi
                                })
                                .collect();
                            segs.sort_by_key(|&j| chart.theta[j].z[0]);
                            let tiled = !segs.is_empty()
                                && chart.theta[segs[0]].z[0] == t_lo
                                && chart.theta[*segs.last().expect("non-empty")].z[1] == t_hi
                                && segs.windows(2).all(|w| {
                                    chart.theta[w[0]].z[1] == chart.theta[w[1]].z[0]
                                        && chart.theta[w[0]].end[1] == chart.theta[w[1]].end[0]
                                });
                            if !tiled {
                                rep.ruling_run_not_tiled += 1;
                                ok = false;
                                continue;
                            }
                            let ordered: Vec<usize> = if up {
                                segs
                            } else {
                                segs.into_iter().rev().collect()
                            };
                            for j in ordered {
                                let [e0, e1] = chart.theta[j].end;
                                pieces.push(Piece {
                                    wall: Wall::Ruling {
                                        cyl: k,
                                        side: sd,
                                        up,
                                    },
                                    ends: if up { (e0, e1) } else { (e1, e0) },
                                });
                            }
                        }
                    }
                }
                if let Some((c, ccw)) = whole_rim {
                    if ccw {
                        rims_lo.push(c)
                    } else {
                        rims_hi.push(c)
                    }
                    continue;
                }
                let n = pieces.len();
                if n < 2 {
                    rep.unchained += 1;
                    ok = false;
                    continue;
                }
                let mut nodes = Vec::with_capacity(n);
                let mut walls = Vec::with_capacity(n);
                for (p, piece) in pieces.iter().enumerate() {
                    if piece.ends.1 != pieces[(p + 1) % n].ends.0 {
                        rep.unchained += 1;
                        ok = false;
                    }
                    nodes.push(piece.ends.0);
                    walls.push(piece.wall);
                }
                let holds_lowest_bottom = cyc
                    .iter()
                    .any(|e| matches!(e, Edge::Arc { line, ccw: true, .. } if *line == lowest));
                if holds_lowest_bottom && outer.is_none() {
                    outer = Some(rings.len());
                }
                rings.push(Ring::new(nodes, walls));
            }
            if !ok {
                continue;
            }
            let (own, other) = reads[members[0]]
                .chamber
                .expect("an emitted cell has a chamber");
            let flip = !crate::bands::keep_for(kind, side, own, other);
            if rims_lo.is_empty() && rims_hi.is_empty() {
                let Some(o) = outer else {
                    rep.unchained += 1;
                    continue;
                };
                let mut rings = rings;
                let outer = rings.remove(o);
                faces.push(LocalFace {
                    surf: ClassIx::Cyl(k),
                    outer: Bound::Ring(outer),
                    inner: rings.into_iter().map(Bound::Ring).collect(),
                    flip,
                });
                rep.ring_faces += 1;
            } else {
                let mut stats = crate::boolean::CurvedStats::default();
                match crate::boolean::classify_cycles(
                    k,
                    flip,
                    rings,
                    rims_lo,
                    rims_hi,
                    &curved.cut_rims,
                    &mut stats,
                ) {
                    Ok(f) => {
                        faces.push(f);
                        rep.band_faces += 1;
                    }
                    Err(_) => rep.classify_abstain += 1,
                }
            }
        }
        rep.faces = faces.len();
        Ok((faces, rep))
    }
}

#[cfg(test)]
pub(crate) mod probe {
    use std::sync::Mutex;

    /// One chart's shape, and how today's emitted lateral faces cover its cells.
    #[derive(Clone, Copy, Debug, Default)]
    pub(crate) struct Row {
        pub(crate) z_lines: usize,
        pub(crate) theta: usize,
        /// How many **rows** (lateral faces) the class carries.
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
            /// Intervals where some wall is crossed an odd number of times (a panel's or a
            /// chain's wall with one ruling in the interval): the walk has an open end there and
            /// the closure is not asserted (E2-2).
            pub(crate) unpaired: usize,
            /// Rulings whose every mark is a **graze** — the face stops at the line rather than
            /// crossing it, so whether it reaches the interval is `face_spans`' question and not a
            /// label's. Existence, not membership (cell ㉒'s split, on the vertical axis).
            pub(crate) grazing_rulings: usize,
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

    /// **The fourth rung's ledger (D2b-0, trimmed in D3)**: what the cell reader read, what the
    /// census's own walk of the chart predicts, and how many faces the emitter put on the class.
    /// One row per chart with cells, whether or not the emitter produced a face there.
    ///
    /// ★ Compared **within a row** (`cells` is copied in), never row-by-row against another
    /// ledger: tests run in parallel and the ledgers interleave independently.
    /// **The vertical answer beside the horizontal one** — the shadow the cutover is measured
    /// against (capability D's own pattern: a rung that ships no capability and *measures*).
    ///
    /// The chart's two axes are symmetric in principle, but only the horizontal ones have ever
    /// decided a face, so the sign bridge on the vertical side ([`crate::arrangement::
    /// plus_theta_is_above`]) is unexercised. These count what it would say.
    pub(crate) mod shadow {
        use std::sync::Mutex;

        /// `(agree, disagree, vertical_only, horizontal_only, both_silent, one_ruling)`
        pub(crate) static COUNTS: Mutex<(usize, usize, usize, usize, usize, usize)> =
            Mutex::new((0, 0, 0, 0, 0, 0));

        pub(crate) fn record(
            horizontal: Option<(bool, bool)>,
            vertical: &[(bool, bool)],
            one_ruling: bool,
        ) {
            let mut g = COUNTS
                .lock()
                .expect("the probe's lock is never held across a panic");
            g.5 += usize::from(one_ruling);
            // The vertical readings must agree with each other before they may agree with anyone.
            let v = match vertical {
                [] => None,
                [a] => Some(*a),
                [a, b] if a == b => Some(*a),
                _ => {
                    g.1 += 1; // two walls of one cell disagreeing is a disagreement of its own
                    return;
                }
            };
            match (horizontal, v) {
                (Some(h), Some(x)) if h == x => g.0 += 1,
                (Some(_), Some(_)) => g.1 += 1,
                (None, Some(_)) => g.2 += 1,
                (Some(_), None) => g.3 += 1,
                (None, None) => g.4 += 1,
            }
        }
    }

    /// The counts behind [`super::OtherWhy`], and the whole-circle disagreements in full: which
    /// test, which line, and the per-arc bits that disagreed.
    pub(crate) mod other {
        use super::super::OtherWhy;
        use std::sync::Mutex;

        pub(crate) static COUNTS: Mutex<Vec<(OtherWhy, usize)>> = Mutex::new(Vec::new());

        #[derive(Clone, Debug)]
        pub(crate) struct Whole {
            pub(crate) test: String,
            pub(crate) cyl: usize,
            pub(crate) t: f64,
            pub(crate) end: usize,
            pub(crate) above: bool,
            pub(crate) bits: Vec<(bool, bool)>,
        }

        pub(crate) static WHOLE: Mutex<Vec<Whole>> = Mutex::new(Vec::new());

        pub(crate) fn record(why: OtherWhy) {
            let mut c = COUNTS
                .lock()
                .expect("the probe's lock is never held across a panic");
            match c.iter_mut().find(|(w, _)| *w == why) {
                Some((_, n)) => *n += 1,
                None => c.push((why, 1)),
            }
        }

        pub(crate) fn record_whole(
            cyl: usize,
            t: nacre_scalar::Rat,
            end: usize,
            above: bool,
            bits: Vec<(bool, bool)>,
        ) {
            record(OtherWhy::WholeDisagree);
            WHOLE
                .lock()
                .expect("the probe's lock is never held across a panic")
                .push(Whole {
                    test: std::thread::current().name().unwrap_or("?").to_string(),
                    cyl,
                    t: t.to_f64(),
                    end,
                    above,
                    bits,
                });
        }
    }

    /// **D5, 1b — the shadow beside today's road**, one row per class.
    pub(crate) mod regions {
        use std::sync::Mutex;

        #[derive(Clone, Debug)]
        pub(crate) struct Row {
            pub(crate) test: String,
            pub(crate) cyl: usize,
            pub(crate) emitter_refused: bool,
            pub(crate) shadow_refused: bool,
            pub(crate) shadow_faces: usize,
            pub(crate) today_faces: usize,
            /// The two face sets agree by signature (bound kinds, flip, holes, node sets).
            pub(crate) agree: bool,
            pub(crate) rep: super::super::regions::Report,
        }

        pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

        pub(crate) fn push(r: Row) {
            ROWS.lock()
                .expect("the probe's lock is never held across a panic")
                .push(r);
        }
    }

    /// **D5, stage 0 — the run-end refusal, attributed.** One row per `RulingBoundNotYet`
    /// raised at `emit_lateral`'s run-end node lookup: which z-line, what that line's arc under
    /// the run's first cell says for this solid, whether the wall class has a node there (a piece
    /// of the same ruling ending at `t`), and whether the cut rim on that line carries the
    /// station's canonical name.
    pub(crate) mod rbny {
        use std::sync::Mutex;

        #[derive(Clone, Debug)]
        pub(crate) struct Row {
            /// The test that ran the boolean (libtest names each test's thread).
            pub(crate) test: String,
            pub(crate) cyl: usize,
            pub(crate) t: f64,
            /// `0` = the interval's low line, `1` = its high line.
            pub(crate) end: usize,
            pub(crate) ruling: Option<(usize, i8)>,
            /// What the run's first cell's end on that line is: `disk` / `exact` / `other` /
            /// `nocircle`.
            pub(crate) end_kind: &'static str,
            /// For an `exact` end: every arc of the run carries, for this solid, only
            /// `Transversal` marks (the face passes through — the line is no boundary there).
            pub(crate) transversal_only: Option<bool>,
            /// For an `exact` end: some arc carries a `Graze` for this solid.
            pub(crate) any_graze: Option<bool>,
            /// A piece of the same `(wall, side)` ruling ends at `t` — the wall class's
            /// arrangement has a node there.
            pub(crate) wall_node: bool,
            /// The cut rim on that line lists the station's canonical name; `None` when no rim
            /// of this cylinder on that line is cut.
            pub(crate) rim_node: Option<bool>,
        }

        pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

        #[allow(clippy::too_many_arguments)]
        pub(crate) fn record(
            jd: &crate::tolerant::Judge<'_, crate::planes::WorkingPlane>,
            k: usize,
            def: &nacre_topo::CylinderDef,
            chart: &super::super::Chart,
            lines: &super::super::Lines,
            curved: &crate::arrangement::Curved,
            side: crate::planes::SolidSide,
            read: &super::super::CellRead<'_>,
            r: usize,
            t: nacre_scalar::Rat,
            e: usize,
        ) {
            use super::super::End;
            use crate::arrangement::SegKind;
            let ruling = chart.ruling_name(jd, k, def, r);
            let (end_kind, transversal_only, any_graze) = match &read.ends[e] {
                End::Disk(_) => ("disk", None, None),
                End::Other => ("other", None, None),
                End::NoCircle => ("nocircle", None, None),
                End::Exact(arcs) => {
                    let mine = || {
                        arcs.iter()
                            .flat_map(|a| a.marks.iter())
                            .filter(move |(s, _)| *s == side)
                            .map(|(_, kind)| *kind)
                    };
                    let n = mine().count();
                    (
                        "exact",
                        Some(n > 0 && mine().all(|k| matches!(k, SegKind::Transversal { .. }))),
                        Some(mine().any(|k| matches!(k, SegKind::Graze { .. }))),
                    )
                }
            };
            let wall_node = ruling.is_some_and(|name| {
                (0..chart.theta.len()).any(|j| {
                    chart.ruling_name(jd, k, def, j) == Some(name) && chart.theta[j].z.contains(&t)
                })
            });
            let rim_node = ruling.and_then(|(wall, sd)| {
                let c = lines
                    .classes(t)
                    .iter()
                    .copied()
                    .find(|&c| curved.cut_rims.contains_key(&(k, c)))?;
                let rim = &curved.cut_rims[&(k, c)];
                let name = crate::arrangement::crossing_on_ruling(jd, def, c, wall, k, sd).ok()?;
                Some(rim.nodes.contains(&name))
            });
            ROWS.lock()
                .expect("the probe's lock is never held across a panic")
                .push(Row {
                    test: std::thread::current().name().unwrap_or("?").to_string(),
                    cyl: k,
                    t: t.to_f64(),
                    end: e,
                    ruling,
                    end_kind,
                    transversal_only,
                    any_graze,
                    wall_node,
                    rim_node,
                });
        }
    }

    pub(crate) mod d2b {
        use std::sync::Mutex;

        #[derive(Clone, Debug, Default)]
        pub(crate) struct Row {
            /// The test that ran the boolean (libtest names each test's thread) — D5's stage 0
            /// attributes the ledger's populations to fixtures by this.
            pub(crate) test: String,
            /// What the emitter refused this boolean with, when it did.
            pub(crate) emission_reason: Option<crate::RejectReason>,
            /// D5 stage 0 (P4): every (z-line, station) pair asked for the station's canonical
            /// name on that line (`crossing_on_ruling`), and how many could not be named.
            pub(crate) station_pairs: usize,
            pub(crate) station_name_failures: usize,
            /// Of `end_other`, the ends of a single-cut sector (both walls one ruling).
            pub(crate) end_other_single_cut: usize,
            pub(crate) cells: usize,
            /// The emitter refused this boolean's lateral faces (some class's cells could not be
            /// read) — the census still records every chart of it.
            pub(crate) emitter_refused: bool,
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
            /// Present cells with a `NoCircle` end — the `circle_on_class` premise at one end.
            pub(crate) nocircle_present: usize,
            /// Cells whose two speaking ends disagreed — the two-end refusal
            /// (`CylinderGateUndecided`), on the chart; the emitter refuses the class.
            pub(crate) src2_disagree: usize,
            /// Cells with a face and no speaking end. Asserted 0 at the record; here as a count
            /// so the reporting test can say the assertion was live.
            pub(crate) src0_present: usize,
            pub(crate) exist_disagree: usize,
            /// The reader refused a cell of this chart by name (`face_spans`' refusal, the
            /// two-cut-ends refusal, or an arrangement-table inconsistency) — the whole class is
            /// refused.
            pub(crate) read_refused: usize,
            /// Both-cut cells the trace says the face is not in — the panel road's dropped sectors.
            pub(crate) exist_marks_false: usize,
            /// Cells the reader would emit, and cells it could not decide.
            pub(crate) emit: usize,
            pub(crate) emit_unknown: usize,
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
            /// Emitted whole-circle cells (`Band`s over intervals with no ruling), and runs the
            /// panel road would cut further because an interior ruling has a node on a cut rim.
            pub(crate) whole_emitted: usize,
            pub(crate) run_split_at_node: usize,
            /// Band-shaped intervals with both rims cut (the emitter sends those to the panel
            /// rings); intervals with more than one partial run; and partial-run pairs that would
            /// have to merge across a non-boundary z-line.
            pub(crate) full_run_both_cut: usize,
            pub(crate) intervals_multi_run: usize,
            pub(crate) run_run_nonboundary: usize,
            /// Consecutive band-shaped intervals with one chamber across a line that is not a
            /// boundary (`boundary_lines`) — the merges the emitter makes at the chart.
            pub(crate) z_merge_bandlike: usize,
            /// Band-shaped neighbours whose chambers differ across a non-boundary line.
            pub(crate) z_flip_nonboundary: usize,
            /// Cut ends read, and those carrying no / more than one lateral mark of their own solid
            /// (the one-mark contract).
            pub(crate) arcs_read: usize,
            pub(crate) arcs_no_mark: usize,
            pub(crate) arcs_multi_mark: usize,
            /// The emitter's lateral faces of this class (0 when it refused).
            /// Σ(len − 1) over `End::Exact` runs: the arcs an end reads **beyond its first**.
            /// A sector that spans a rim node the chart has no line for reads the whole run.
            pub(crate) exact_run_arcs: usize,
            pub(crate) emitted_faces: usize,
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
                ..
            } = *r;
            assert!(n >= 1, "a class with no row: {r:?}");
            assert!(z_lines >= 2, "a band needs two boundaries: {r:?}");
            // ★ `theta % 2 == 0` ("a wall cuts two rulings") was asserted here until E2-2: a
            // band's lateral reaches both rulings of a wall, a panel's or a chain's may reach one.
            let _ = theta;
            assert!(!refused, "a chart's theta order could not be formed: {r:?}");
            assert!(cells >= 1, "two z-lines are one interval: {r:?}");
            // A chart has `z_lines - 1` intervals, and an interval carries either no ruling or
            // some number of them — so these two counts are disjoint subsets of that many.
            assert!(
                r.whole_circle + r.odd_k < z_lines,
                "more intervals than the chart has: {r:?}"
            );
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
                r.intervals_with_flip <= r.closes + r.does_not_close + r.unpaired,
                "more flipping intervals than walked ones: {r:?}"
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
    }

    /// **The cells read their chamber off the horizontal lines, and the emitter is predicted by the
    /// census's own walk** (capability D, D2b-0 → D3).
    ///
    /// The absolute claims are asserted in `census`, where the facts are made: a present cell
    /// always has a speaking end (`src0_present == 0`, the reason no wall-side sign is needed),
    /// and the ends' bookkeeping is total. What this holds, universally over every recorded row
    /// so no interleaving can break it, is the rest of the forecast:
    ///
    /// * `emitted_faces == whole_emitted + full_runs − z_merge_bandlike + partial_runs` — the
    ///   headline since D3: the emitter's face count per class is what the census's own walk of
    ///   the chart predicts (vacuous where the emitter refused).
    /// * `emit_unknown == 0 || emitter_refused` — the emitter never puts a face over a present
    ///   cell it could not read, e.g. one whose two ends contradict (the `wal corner-lo` corpus
    ///   family has such cells — `src2_disagree` 8 per row — and is refused by name).
    /// * `exist_disagree == 0`, `read_refused == 0` — the trace and the span tell the same
    ///   existence story wherever both speak.
    /// * `nocircle_present == 0`, `z_flip_nonboundary == 0` — a present cell has a circle at both
    ///   ends, and a band's chamber never flips across a line that is not a boundary (D1b's extra
    ///   lines are harmless).
    /// * `arcs_read == end_exact + exact_run_arcs`, `arcs_multi_mark == 0` — every cut-end read
    ///   carries at most one lateral mark of its own solid (the band road's MARKS contract), and
    ///   an end that spans a run of the rim's arcs reads every one of them. (`arcs_no_mark` is no
    ///   longer 0: a ⊥ cap through a notch reads cut ends inside the lateral's hole.)
    /// * `partial_theta_in_disk_interval == 0` — a disk-ended interval's kept sectors always close
    ///   into a whole circle, so it is emitted as one `Band`.
    ///
    /// And the counters are not vacuous: the fixtures below put every end kind but `Other` on the
    /// ledger, the chained ones drop sectors for existence (**≥ 4**, the band road's own count,
    /// `tests::a_chained_cylinder_bounded_by_the_first_builds`), both merge directions are seen,
    /// arcs are read and faces are emitted. `end_other` is **unmeasured** if it stays 0 — the
    /// θ-placement it would need is not built, and a zero here is a population claim only once a
    /// fixture reaches it.
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
        // ★ A **through-axis wall** on the wall boss: its plane holds the axis, so it cuts the
        // lateral along two rulings — and the chart has a θ line for only the one the face
        // reaches. The sector between them spans a rim node, which is the only shape that makes
        // an end read a *run* of arcs. The **Cut** is the kind that builds through it (the wall
        // boss's Fuse waits on the emitter's rim-station ladder — `RulingBoundNotYet`). Without this the run equality below is vacuous (measured:
        // the population above produced `exact_run_arcs == 0`).
        {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(Point3::from_array([2.0, 0.0, -1.0]), up, 0.5, 4.0);
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Cut, plate, boss).expect("the notch builds");
            m.rebuild_adjacency();
            let wall = m.add_cuboid(
                Point3::from_array([2.0, -3.0, -2.0]),
                Point3::from_array([5.0, 3.0, 6.0]),
            );
            m.rebuild_adjacency();
            crate::boolean(&mut m, BoolKind::Cut, out[0], wall)
                .expect("the through-axis wall cuts");
        }
        let wall_boss = ([2.0, 0.0, -1.0], 4.0, BoolKind::Fuse);
        for second in [
            ([6.0, 2.0, -1.0], 4.0, BoolKind::Cut),
            ([6.0, 2.0, 2.0], 1.0, BoolKind::Fuse),
            // A short boss whose cap (z = 2.5) is a ⊥ line strictly inside the wall boss's span:
            // the wall boss's circle is traced there but bounds no face, so the line is not a
            // boundary — the band-shaped merge population.
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
            assert!(
                r.emit_unknown == 0 || r.emitter_refused,
                "a cell's chamber could not be read, yet the emitter emitted: {r:?}"
            );
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
            // D3's successors of the reference-road checks, held as counts (measured 0 first).
            assert_eq!(
                r.nocircle_present, 0,
                "a cell with a face meets a line with no circle: {r:?}"
            );
            assert_eq!(
                r.z_flip_nonboundary, 0,
                "the chamber changes across a line that is no boundary: {r:?}"
            );
            // ★ `arcs_no_mark` is no longer held at 0 here: a ⊥ cap through a wall boss's notch
            // reads cut ends inside the lateral's hole, and those carry no mark by right (E3).
            // The record site asserts the true proposition — none only where the cell is absent.
            assert_eq!(
                r.arcs_multi_mark, 0,
                "a cut end read carries two lateral marks of one solid: {r:?}"
            );
            // ★★★★★ **A run expired the old equality.** This read `arcs_read == end_exact` until
            // a sector was allowed to span several of the rim's arcs; the ledger then measured
            // 4,370 against 4,312 over the whole suite while this assertion stayed green, because
            // it sees only its own fixtures (the recorded trap: a lock in a `#[test]` watches what
            // ran before it). The proposition that is still true — and still crosses the two loops
            // that count these, which is what it is for — carries the run's extra arcs by name.
            assert_eq!(
                r.arcs_read,
                r.end_exact + r.exact_run_arcs,
                "every exact end is a cut end read, and only those: {r:?}"
            );
            // ★ The emitter, predicted by the census's own walk (the successor of the shadow's
            // face-by-face comparison): band-shaped intervals minus the merges across
            // non-boundary lines, plus the partial runs. Vacuous where the emitter refused.
            if !r.emitter_refused {
                assert_eq!(
                    r.emitted_faces,
                    r.whole_emitted + r.full_runs - r.z_merge_bandlike + r.partial_runs,
                    "the emitter's face count is not the census's: {r:?}"
                );
            }
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
        }
        // ── D2b-1: what the cutover would emit, held as counts (their production twins are
        // honest rejects, so they are not asserted where the fact is made).
        for r in &rows {
            // ★ Restated (D5, 1a): the bare 0 was true only while those run ends read `Other`
            // — a station placed under another line's name that the θ order refused — so the
            // census's own walk never reached the boundary ruling. Placed by name the ends read
            // `Exact`, the walk sees a boundary ruling with no rim node (the crossing census's
            // whole `RulingBoundNotYet` column), and what is true of it is that the emitter
            // refuses such a run rather than building it. Asserted at the fact too (`census`).
            assert!(
                r.run_boundary_no_node == 0 || r.emitter_refused,
                "a run's boundary ruling has no rim node and the emitter built it: {r:?}"
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
        assert!(sum(&rows, |r| r.arcs_read) > 0, "no cut end was ever read");
        assert!(
            sum(&rows, |r| r.exact_run_arcs) > 0,
            "no end ever read a run of the rim's arcs"
        );
        assert!(
            sum(&rows, |r| r.emitted_faces) > 0,
            "the emitter never emitted a face"
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
