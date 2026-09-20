//! **A cylinder's own chart, as an arrangement's line set** (capability D, first rung).
//!
//! The plane side long ago stopped walking faces by hand: each plane class gets a **cell complex**
//! (`arrangement`'s `walk_cells → nest_cells → label_cells → emit_faces`) and the case-work went
//! with it. The cylinder side follows: this module emits the lateral's faces, and the merge
//! reads one rule (a cycle's winding, `seam_step`); of the hand-written walks only `boolean`'
//! `band_loop` slit-leg remains, generalized to a chain rim.
//!
//! `docs/design.md` names the way out: the lateral has an **isometric chart** `(z, r·θ)`, so the
//! same engine can run there — the seam is only the chart's cut line, and notches, holes and
//! θ-panels all become cells. This module is that road's **first** piece: the arrangement's two
//! axes, and nothing else.
//!
//! ## Why the lines are orthogonal — a derivation, not a measurement
//!
//! The population gate (`planes`) sorts every `(plane class, cylinder class)` pair into
//! five, and only two leave a mark on the lateral:
//!
//! | plane vs axis | distance | on the lateral | ☑ suite |
//! |---|---|---|---|
//! | `n ∥ m` | — | a **circle** — the chart's *horizontal* | 1128 |
//! | `n ⊥ m` | `0 ≤ d < r` | **two rulings** — the chart's *vertical* | 75 |
//! | `n ⊥ m` | `d = r` | **nothing** from another solid's plane (a tangency grazes and divides nothing); the face's **own** tangent wall is a **station with no crossing** — a ruling of `side = 0`, where the face ends | — |
//! | `n ⊥ m` | `d > r` | nothing | 1494 |
//! | oblique | — | refused, `ObliqueCylinderCut` | 2 |
//!
//! So within this milestone the chart carries **axis-parallel segments only**. Not a grid, though:
//! a ruling exists over a finite axis interval and a cut circle is arcs, so the cells are what
//! rectilinear *segments* cut out — which is why building them is its own rung.
//!
//! ★ The gate is the reason, so when the gate changes this does too; the counts above predate
//! the offset and tangent walls. The **offset** wall (`0 < d < r`) adds
//! vertical lines at *irrational* θ; the chart never reads a θ **value**, only an order, so
//! nothing breaks. The **tangent** wall (`d = r`)
//! adds nothing at all: a grazing line stations no sector, so the gate records no
//! crossing and this table's third row is «nothing», not «one vertical». ★ The row's
//! other half: a fillet's or a slot's wall is tangent to its **own** cylinder, and that ruling is
//! not a crossing but the face's edge — a station of `side = 0` carrying no label (the seated
//! face states the piece itself, `SegKind::Tangent`), which the cell reader takes as «the face
//! ends here» and the vertical-answer assertion exempts. ★ And a third reading: when the
//! *other* solid's plane runs through the axis, one of its two verticals **is** that tangent
//! ruling — a line two planes and the cylinder share. Its ends have one name (the alias table's,
//! seeded from the operands — `Curved::aliases` is what this chart's stations are canonicalized
//! through), and the line is stated in the plane vocabulary by both faces that touch it — the
//! tangent wall's run and the lateral's own graze — so the ruling vocabulary states nothing there.
//!
//! ## What is here, and what is not
//!
//! The two axes, the **cells** they cut, the vertical lines' answers, the
//! **cell reader** ([`Chart::read_cell`]), which reads a cell's chamber off the lines at
//! its ends (`DiskLabels`/`ArcLabels` through `bands::read_bits`, the rulings' labels where the
//! rims are silent) and its existence off the trace (`face_spans`), and the **emitter**
//! ([`emit_lateral`]) — the **region walk** ([`regions`]): emitted cells → connected
//! components → each component's boundary on the chart's grid → runs cut into the neighbouring
//! class's own pieces → cycles → a `Bound`. The census holds the chart against its own rules
//! and the emitter's faces against the reads they came from ([`census`]).

use crate::arrangement::{ArcLabel, Curved, Label, RulingExtent};
use crate::boolean::{Bound, LocalFace};
use crate::planes::{ClassIx, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, combinatorics, reject};
use nacre_scalar::Rat;

#[cfg(test)]
mod census;
mod emit;
#[cfg(test)]
pub(crate) mod probe;
mod read_cell;
pub(crate) mod regions;

#[cfg(test)]
pub(crate) use census::*;
pub(crate) use emit::*;
pub(crate) use read_cell::*;

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
    /// The ruling's side on its wall — `0` a **tangent** station: a line the face ends
    /// on, with no chamber to read behind it, so no `label`.
    pub(crate) side: i8,
    /// The piece's own two pierce nodes: the θ **order** is asked of these
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

/// **Which of a wall's two rulings a pierce node lies on** — `arrangement::node_ruling_side`, the
/// one spelling (`ruling_side` asked of a name) the assembly's `Wall::Ruling { side }` reads.
/// ★ This used to be a second copy of that function's body: the same `pierce_meet` +
/// `ruling_side`, spelled twice one module apart.
///
/// ★★★★★ **Not the node's `QuadRoot`.** A `Pierce` name's root is `Lo`/`Hi` along
/// `ℓ = n₁ × n₂` of *its own* plane pair, so the same physical ruling reads `Hi` where its node
/// pairs the wall with the plate's top and `Lo` where it pairs it with the boss's cap (the pair
/// order and the ⊥ normal's sense both reverse `ℓ` — `QuadRoot::canonical`'s doc). ☑ Measured on
/// the corner boss: two different rulings of one wall carried one name and
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
    let (_, cyl, _) = combinatorics::pierce_name(n)?;
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
/// deriving it is the reader's job — building it here would be a field with no consumer.
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
    /// that was the key the band road's `panel_faces` joined on — but then the ruling's *label* could not be
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
    /// node there (`NodeId::pierce` folds the pair and the root), so «does this rim carry this
    /// station» is **name equality**, and a station the rim does not carry joins the θ order under
    /// the name a piece ending there would have worn.
    ///
    /// ★★★★★ **Not a node of another line**. Placing such a station by
    /// its piece's `end[0]` — a node on a *different* z-line at the same θ — would, where the rim
    /// does carry the station (a tool wall crossing the circle there), hand `circular_order` one
    /// point under two names, which it refuses as the coincidence it checks for. Every one of
    /// those ends would read `End::Other` (one point, two names; 321 ends over the suite).
    fn station_on(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
        c: usize,
        i: usize,
        aliases: &crate::arrangement::Aliases,
    ) -> Option<combinatorics::NodeId> {
        let (wall, side) = self.ruling_name(jd, k, def, i)?;
        // As the table knows it: a station on a corner is the corner.
        crate::arrangement::crossing_on_ruling(jd, def, c, wall, k, side)
            .ok()
            .map(|n| aliases.canon_point(n))
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
/// sources, as axis parameters: the rim classes ([`rim_classes`]) and
/// every row's own span ends. A band-shaped emission may continue across any other line, so
/// [`emit_lateral`] merges there and nowhere else, and the census asserts the chart's z-lines ⊇
/// this list where it is made.
#[cfg(test)]
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
/// ★★★★★ **The sources are the band road's, and one thing that road did is
/// deliberately *not* done here.** It ended with `keyed.retain(|(t, _)| t ∈ row.span)`:
/// it clipped to the **face** it was asked about, because its row was a face and not a class. A
/// class has no such span — one lateral
/// surface can carry several faces with gaps between them, which is exactly why `CylRow` became
/// per-face. So the clip is dropped, and that is a *repair* rather than a loss: `CylRow`'s own doc
/// warns that merging the spans into one `min..max` would **invent a band where the solid has no
/// face at all**, and on a chart that gap is simply a cell that keeps nothing — the argument the
/// plane side already makes for its unbounded cells.
///
/// ★ The cylinder class itself is one solid's surface, not both operands': `planes` interns
/// cylinders by `Handle<Surface>` (planes by geometry), and two solids share no handles. That is
/// sound here only because a **coincident pair is refused** — one handle carrying rows of both
/// solids, or one surface stated under two handles (`same_surface`), both `CylinderPairContact` at
/// the gate. What does see both
/// operands is what the chart is made *of*: the plane classes, which are interned geometrically.
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
    // so sorting `z` alone would silently break that pairing — and the cell reader asks
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
                        side: r.side,
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

#[cfg(test)]
#[path = "../tests/cyl_chart.rs"]
mod tests;
