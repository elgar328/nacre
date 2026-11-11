//! **The arrangement's second arm: a cylinder's own chart.**
//!
//! The plane side does not walk faces by hand: each plane class gets a **cell complex** (`super`'s
//! `walk_cells → nest_cells → label_cells → emit_faces`) in place of per-face case-work. The
//! cylinder side does the same on its own grid: the lateral has an **isometric
//! chart** `(z, r·θ)` where the seam is only the chart's cut line, so notches, holes and
//! θ-panels all become cells. This module runs the whole arm — the two axes, the cells
//! (`read_cell`), the regions walk (`regions`) and the emit (`emit`) — and a region's bounds are
//! its cycles as they stand, the outer one picked by its winding (its arc units), with no seam
//! edge written anywhere.
//!
//! **Same stages, not same code.** The chart is an annulus where the plane side is an unbounded
//! DCEL, so it runs the stages on its own cells rather than calling the plane side's — `regions`
//! states what that costs and why.
//!
//! ## What the lines are — a derivation, not a measurement
//!
//! A curve on the lateral is a ruling (θ fixed) or a graph over θ: a plane not along the axis
//! meets the cylinder in `z(θ) = a + b·cosθ + c·sinθ`, a circle at one height when it is across
//! the axis. The chart's **lines** are such sections, and all it asks of two of them is which is
//! higher on a generator ([`combinatorics::line_order`] — one comparator); a slab's order and a
//! generator's heights hold the answers ([`Chart::slabs`], [`Event::height`]). The population gate
//! (`planes`) sorts every `(plane class, cylinder class)` pair into five, and three leave a mark
//! on the lateral:
//!
//! | plane vs axis | distance | on the lateral | suite |
//! |---|---|---|---|
//! | `n ∥ m` | — | a **circle** — the chart's *horizontal* | 22,220 |
//! | `n ⊥ m` | `0 ≤ d < r` | **two rulings** — the chart's *vertical* | 4,616 |
//! | `n ⊥ m` | `d = r` | **nothing** from another solid's plane (a tangency grazes and divides nothing); the face's **own** tangent wall is a **station with no crossing** — a ruling of `side = 0`, where the face ends | 2,062 |
//! | `n ⊥ m` | `d > r`, or no face reaching the lateral | nothing | 24,661 |
//! | oblique | — | a **section** — a line of the chart, crossing others, which the class draws on its own arrangement (`CylRecord::Sections`) wherever the plane touches none of the lateral's circles; a cell it bounds on one side only is a face with an elliptical edge | 298 |
//!
//! So the chart's lines are the sections the arrangement traced on this cylinder — every class with
//! a label or a cut rim of it — and its events the rulings' stations,
//! the rims' cuts and the points two lines cross or touch at. Not a grid, though: a ruling exists
//! between two lines and a cut circle is arcs, so the cells are what those leave between the
//! lines — a slab crossed with a gap.
//!
//! The gate is the reason, so when the gate changes this does too. The counts above are the gate's
//! verdicts over the nacre-ops unit suite, one per (named plane class, cylinder) pair each
//! boolean's gate sorts — the oblique row's are sections, beside 115 oblique classes with no face
//! reaching the lateral and 16 the gate refuses as sections no road draws. The **offset** wall (`0
//! < d < r`) adds vertical lines at *irrational* θ; the chart never reads a θ **value**, only an
//! order, so nothing breaks. The **tangent** wall (`d = r`) adds nothing at all: a grazing line
//! stands at no event, so the gate records no crossing and this table's third row is «nothing», not
//! «one vertical». The row's other half: a fillet's or a slot's wall is tangent to its **own**
//! cylinder, and that ruling is not a crossing but the face's edge — a station of `side = 0`
//! carrying no label (the seated face states the piece itself, `SegKind::Tangent`), which the cell
//! reader takes as «the face ends here» and the vertical-answer assertion exempts. And a third
//! reading: when the *other* solid's plane runs through the axis, one of its two verticals **is**
//! that tangent ruling — a line two planes and the cylinder share. Its ends have one name (the
//! alias table's, seeded from the operands, which every name the chart reads already wears), and
//! the line is stated in the plane vocabulary by both faces that touch it — the tangent wall's run
//! and the lateral's own graze — so the ruling vocabulary states nothing there. The same holds on
//! the *other* solid's tangent wall where the face's own edge lies on its tangent ruling (a half
//! cylinder's flat edge): the flat and the lateral each state their side of the plane-pair line on
//! that wall's class (`rulings::lateral_half_of_plane_pair_lines`), and the chart takes nothing
//! from it.
//!
//! ## What is here, and what is not
//!
//! The two axes, the **cells** they cut, the vertical lines' answers, the
//! **cell reader** ([`Chart::read_cell`]), which reads a cell's chamber off the lines at
//! its ends and walls (`DiskLabels`/`ArcLabels` and the rulings' labels, all through
//! `read_cell::read_bits`, which must agree) and its existence off the trace (`face_spans`), and
//! the **emitter**
//! ([`emit_lateral`]) — the **region walk** ([`regions`]): emitted cells → connected
//! components → each component's boundary on the chart's grid → runs cut into the neighbouring
//! class's own pieces → cycles → a `Bound`. The census holds the chart against its own rules
//! and the emitter's faces against the reads they came from (`census`).

use crate::arrangement::{ArcLabel, Curved, DiskLabel, Label, SegKind};
use crate::draft::LocalFace;
use crate::planes::{SolidSide, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, combinatorics, reject};
use core::cmp::Ordering;
use nacre_exact::Rat;
use nacre_exact::quad::Toward;

#[cfg(test)]
#[path = "../../tests/probes/chart_census.rs"]
mod census;
mod emit;
#[cfg(test)]
#[path = "../../tests/probes/chart_probe.rs"]
pub(crate) mod probe;
mod read_cell;
pub(crate) mod regions;

#[cfg(test)]
pub(crate) use census::*;
pub(crate) use emit::*;
pub(crate) use read_cell::*;

/// One **vertical** segment: a ruling, and the axis interval it spans.
///
/// Carried straight from [`RulingExtent`](crate::arrangement::RulingExtent) — the plane pass decides this and the chart must not
/// decide it again.
#[derive(Clone, Debug)]
pub(crate) struct ThetaSeg {
    pub(crate) wall: usize,
    /// The ruling's side on its wall — `0` a **tangent** station: a line the face ends
    /// on, with no chamber to read behind it, so no `label`.
    pub(crate) side: i8,
    /// The piece's own two pierce nodes: the θ **order** is asked of these
    /// (`circular_order_about_seam` via `combinatorics::theta_groups`), never of a coordinate.
    pub(crate) end: [combinatorics::NodeId; 2],
    /// The chart lines its two ends lie on — indices into [`Chart::lines`], the lower first. Read
    /// off the ends' names: a ruling ends where a plane of its end's name meets it, and that plane
    /// is a line of the chart.
    pub(crate) on: [usize; 2],
    /// **The chart's vertical answer** — the label of the cell this ruling borders inside the
    /// cylinder, carried straight from [`RulingExtent::label`](crate::arrangement::RulingExtent::label). The horizontal lines' answers are
    /// `DiskLabels`/`ArcLabels`; this is the half they do not carry, and [`Chart::read_cell`] reads
    /// it beside them — the cell's chamber stands only where the two agree.
    pub(crate) label: Option<Label>,
    /// The label of the cell beside this ruling **outside** the cylinder — [`RulingExtent::outer`](crate::arrangement::RulingExtent::outer).
    /// A cell's chamber never reads it; [`Chart::slit_at`] does, where the lateral keeps the
    /// outside.
    pub(crate) outer: Option<Label>,
    /// Who traced this line — [`RulingExtent::marks`](crate::arrangement::RulingExtent::marks). Membership is `label`'s question; whether
    /// this lateral face is even here is this one's.
    #[cfg(test)]
    pub(crate) marks: Vec<(SolidSide, SegKind)>,
    /// **Every secant wall through this very ruling, in chord order** — the line two or more
    /// secant classes share on the cylinder ([`crate::planes::SharedRuling`]) is one station, and
    /// the chamber beside the lateral on each side of it is closed by a different one of their
    /// planes ([`Station::chords`]). `None` for a station one wall alone runs through.
    pub(crate) station: Option<Station>,
}

/// **The secant walls of one station** — see [`ThetaSeg::station`].
#[derive(Clone, Debug)]
pub(crate) struct Station {
    /// **The walls in chord order.** Inside the cylinder every secant plane through the line runs
    /// from it as a chord, and all those chords enter the disk: turning from the lateral's `+θ`
    /// tangent `τ = m × r` inward (counter-clockwise about `m`, towards `−r`), they come one after
    /// another within a half turn, so the order is total ([`chord_before`]). The first chord closes
    /// the chamber beside the `+θ` sector, the last the one beside `−θ`, and each two neighbours
    /// the chamber between them — OCCT's `GetFaceOff` picks the next face about an edge by the same
    /// angle, asked here as an exact sign. The station's own wall is one of them.
    ///
    /// Outside, the line cuts the cross-section into as many chambers again: the one between two
    /// neighbouring walls' outer halves is the vertically opposite angle of the one between their
    /// chords. Every chamber between chords `k` and `k + 1` lies on chord `k`'s far side from `τ`
    /// and on chord `k + 1`'s near side; its outer opposite the other way round.
    pub(crate) chords: Vec<Chord>,
}

/// One secant wall of a [`Station`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct Chord {
    pub(crate) wall: usize,
    /// Which of the wall's two rulings the line is.
    pub(crate) side: i8,
    pub(crate) says: ChordSays,
}

/// **What a secant wall states on the station's line over the station's axis interval.**
#[derive(Clone, Copy, Debug)]
pub(crate) enum ChordSays {
    /// A ruling piece over the same interval, end for end: its labels inside and outside the
    /// cylinder — [`ThetaSeg::label`] and [`ThetaSeg::outer`] of that piece.
    Piece {
        label: Option<Label>,
        outer: Option<Label>,
    },
    /// **No piece over this interval.** The wall's class traces a lateral's ruling only where the
    /// lateral's boundary does not run along a plane pair it holds with another class; where it
    /// does, the line is the lateral's own edge at that height and the class states it in the
    /// plane vocabulary alone (`rulings::cycle_on_wall`'s `quiet` runs). A record is written only for a line
    /// on a lateral *face* (`ruling_on_a_face`), and a tangent class is no secant chord, so this is
    /// the one way a secant wall of a recorded line draws nothing there (measured: every one in the
    /// suite, the census and the ignored sweep). The chord still closes its chamber — the two
    /// planes and the axis say so — but no label of it can be read: the side it bounds is silent,
    /// and the lateral ends at the line there.
    NoPiece,
}

/// One cylinder class's chart: its lines, its rulings and events, and the slabs' orders.
#[derive(Clone, Debug, Default)]
pub(crate) struct Chart {
    /// **The lines**: the plane classes the chart carries, each a graph `z(θ)` on the lateral — a
    /// class across the axis is a circle at one height. Ordered nowhere: their order is a slab's
    /// ([`Chart::slabs`]).
    pub(crate) lines: Vec<usize>,
    pub(crate) theta: Vec<ThetaSeg>,
    /// **Every θ something happens at, in θ order** ([`combinatorics::theta_groups`], the seam
    /// first): a ruling's station, a node a rim is cut at, a point two lines cross or touch at.
    /// Between two consecutive events (a *slab*) no line is cut, no two lines meet and no ruling
    /// stands, so a slab of a gap is one cell, and a rim's arc covering it is one arc.
    pub(crate) events: Vec<Event>,
    /// For each class whose rim the split cut on this cylinder, the event of each of its
    /// [`crate::draft::CutRim`] nodes, in the rim's own order.
    pub(crate) rim_events: std::collections::HashMap<usize, Vec<usize>>,
    /// **The lines of each slab, lowest first** — `slabs[u]` the order just after event `u`
    /// ([`combinatorics::line_order`], the one comparator), which holds until the next event; one
    /// slab, asked at no generator, when there is no event.
    pub(crate) slabs: Vec<Vec<usize>>,
}

/// **One event of the chart** — a generator of the cylinder where a ruling stands or a rim is cut.
#[derive(Clone, Debug)]
pub(crate) struct Event {
    /// A point on the generator — where the lines are compared.
    pub(crate) at: combinatorics::ExactPoint,
    /// The rulings standing on it — indices into [`Chart::theta`], one name among them
    /// ([`Chart::ruling_name`]), between several pairs of lines.
    pub(crate) rulings: Vec<usize>,
    /// Each line's height on this generator, as a rank: lines that meet here share one.
    pub(crate) height: Vec<usize>,
}

/// A ruling's **name**: the wall class that cut it, and which of the two roots it is.
///
/// One name, not two. `RulingExtent` also carries a `side`, but that is *derived* from the same
/// pair (`combinatorics::ruling_side_signed`); a second spelling of one identity is where two copies
/// drift apart.
type RulingName = (usize, i8);

/// **Which of a wall's two rulings a pierce node lies on** — `arrangement::node_ruling_side`, the
/// one spelling (`ruling_side` asked of a name) the chart's `Wall::Ruling { side, plane }` states,
/// `plane` the wall it was asked against. A copy of that function's body here (the same
/// `pierce_meet` + `ruling_side`) would be a second spelling one module apart.
///
/// **Not the node's `QuadRoot`.** A `Pierce` name's root is `Lo`/`Hi` along
/// `ℓ = n₁ × n₂` of *its own* plane pair, so the same physical ruling reads `Hi` where its node
/// pairs the wall with the plate's top and `Lo` where it pairs it with the boss's cap (the pair
/// order and the ⊥ normal's sense both reverse `ℓ` — `QuadRoot::canonical`'s doc). Measured on
/// the corner boss: two different rulings of one wall carried one name and
/// the panel join claimed the wrong sectors. The side is a fact about the *point*, so every piece
/// of one ruling answers alike, whatever its end nodes pair it with.
fn ruling_side_of(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    wall: usize,
    n: combinatorics::NodeId,
) -> Result<i8, RejectReason> {
    // A ruling piece's node is a pierce name of this chart's own cylinder, and every class has a
    // narrow world statement past the cylinder gate — so a node of another kind or another
    // cylinder (M6b's pair, which no piece carries), or a class with none, is the split and the
    // chart disagreeing. What is left is `pierce_meet`'s and the side's `None`: an overflow.
    let stages = RejectReason::CylinderStagesDisagree;
    let (_, cyl, _) = combinatorics::pierce_name(n).ok_or(stages)?;
    if cyl != k {
        return Err(stages);
    }
    jd.world_rat(wall).ok_or(stages)?;
    crate::arrangement::node_ruling_side(jd, def, wall, n).ok_or(RejectReason::WitnessNotRational)
}

/// **One cell of the chart** — an axis interval crossed with a slab, the θ between two
/// consecutive [`Chart::events`].
///
/// Addressed by `(slab, gap)` rather than kept as a flat list, because the region walk
/// ([`regions`]) asks for **neighbours**: crossing a line moves to the gap above or below,
/// crossing an event to the slab beside. That adjacency is *derivable* from this address,
/// and deriving it is the reader's job — building it here would be a field with no consumer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Cell {
    /// The slab from event `slab` to the next one counter-clockwise; with no event at all, `0`
    /// and the whole circle. With one event the slab runs from it round to it again.
    pub(crate) slab: usize,
    /// The gap between the slab's lines `gap` and `gap + 1`, lowest first.
    pub(crate) gap: usize,
    /// The two lines bounding it, below and above — indices into [`Chart::lines`].
    pub(crate) lines: [usize; 2],
    /// The rulings standing at the slab's two events across this gap **as indices into
    /// [`Chart::theta`]** — `None` at an event where no ruling spans it (a rim cut there, or a
    /// ruling between other lines). With one event both sides are the same ruling.
    ///
    /// **Indices, not names.** Kept as `(wall, side)` names, the ruling's *label* could not be
    /// reached from a cell without searching for it. The name is one call away
    /// ([`Chart::ruling_name`]).
    pub(crate) walls: [Option<usize>; 2],
}

impl Chart {
    /// **The cells**: every slab crossed with every gap between its lines, gap by gap (lowest
    /// first), slab by slab within a gap.
    ///
    /// **Every ruling ends on a line** — a ruling's node is named `[wall, other plane]`, and the
    /// other plane meets the lateral there, so it is a line ([`chart_of`] reads which). So a ruling
    /// standing at an event spans a gap whole or not at all.
    pub(crate) fn cells(&self) -> Vec<Cell> {
        let n = self.events.len();
        let mut out = Vec::new();
        for g in 0..self.lines.len().saturating_sub(1) {
            for u in 0..n.max(1) {
                let lines = [self.slabs[u][g], self.slabs[u][g + 1]];
                let walls = if n == 0 {
                    [None, None]
                } else {
                    [self.wall_at(u, lines), self.wall_at((u + 1) % n, lines)]
                };
                out.push(Cell {
                    slab: u,
                    gap: g,
                    lines,
                    walls,
                });
            }
        }
        out
    }

    /// The ruling standing at event `e` across the gap between `lines` — one whose ends lie at or
    /// beyond both, on that generator; `None` where the two lines meet there (the gap pinches to
    /// a point) or no ruling spans them.
    pub(crate) fn wall_at(&self, e: usize, [lo, hi]: [usize; 2]) -> Option<usize> {
        let h = &self.events[e].height;
        if h[lo] >= h[hi] {
            return None;
        }
        self.events[e].rulings.iter().copied().find(|&j| {
            let [a, b] = self.theta[j].on;
            h[a] <= h[lo] && h[hi] <= h[b]
        })
    }

    /// Where line `l` stands in slab `u`'s order.
    pub(crate) fn rank(&self, u: usize, l: usize) -> usize {
        self.slabs[u]
            .iter()
            .position(|&x| x == l)
            .expect("every slab orders every line")
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
    ) -> Result<RulingName, RejectReason> {
        let t = &self.theta[i];
        Ok((t.wall, ruling_side_of(jd, k, def, t.wall, t.end[0])?))
    }
}

/// **Every section the arrangement traced on this cylinder, and every ruling on it** — per cylinder
/// *class*.
///
/// **No clip to a face.** Clipping the lines to a face's ends would answer about a face; a class
/// has no such ends — one lateral surface can carry several faces with gaps between them, which is
/// why `CylRow` is per-face. So there is no clip: `CylRow`'s own doc warns that merging the spans
/// into one `min..max` would **invent a band where the solid has no face at all**, and on a chart
/// that gap is simply a cell that keeps nothing — the argument the plane side already makes for its
/// unbounded cells.
///
/// The cylinder class itself is one solid's surface, not both operands': `planes` interns
/// cylinders by `Handle<Surface>` (planes by geometry), and two solids share no handles. That is
/// sound here only because a **coincident pair is refused** — one handle carrying rows of both
/// solids, or one surface stated under two handles (`same_surface`), both `CylinderPairContact` at
/// the gate. What does see both
/// operands is what the chart is made *of*: the plane classes, which are interned geometrically.
pub(crate) fn chart_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    k: usize,
    curved: &Curved,
) -> Result<Chart, BoolError> {
    let def = &cyls[k].def;
    // **The lines are the curves the arrangement traced on this cylinder** — every class with a
    // label or a cut rim of it. Not every class not along the axis: the chart reads a cell off the
    // traces at its ends («a face that is here has a circle at both its ends»), and a section no
    // class drew carries none. A class across the axis traces its circle within every face's span
    // (`circle_on_class`), so those are all there, a rim no face lands on included; an oblique class
    // traces its section where the gate recorded it (`CylRecord::Sections`). Measured: taking
    // every oblique class instead, a section the gate passed as clear of the lateral ran through a
    // fillet's span outside its quarter, and a cell between two such lines read present with no
    // side to say what it held.
    let mut lines: Vec<usize> = curved
        .disk_labels
        .keys()
        .chain(curved.arc_labels.keys())
        .chain(curved.split_rims.keys())
        .filter(|&&(kk, _)| kk == k)
        .map(|&(_, c)| c)
        .collect();
    lines.sort_unstable();
    lines.dedup();
    // A traced class with no world statement is width; one along the axis leaves rulings, never a
    // closed trace, so a label filed under it is the stages disagreeing.
    for &c in &lines {
        let Some(w) = jd.world_rat(c) else {
            return Err(reject(RejectReason::WitnessNotRational));
        };
        if nacre_exact::axis_relation(&[w[0], w[1], w[2]], &def.dir())
            == nacre_exact::AxisRelation::Along
        {
            debug_assert!(
                false,
                "a label of cylinder {k} filed under class {c}, which runs along its axis"
            );
            return Err(reject(RejectReason::CylinderStagesDisagree));
        }
    }
    let line_of = |c: usize| lines.iter().position(|&x| x == c);
    let height = |a: usize, b: usize, at: &combinatorics::ExactPoint| {
        combinatorics::line_order(jd, def, lines[a], lines[b], Some(at), Toward::At).map_err(reject)
    };

    // **A ruling's ends are on lines, read off their names**: the plane of an end's name that is a
    // line is where the ruling stops. A ruling piece states its ends in ascending axis order, so
    // ends that are not — or are level — are the stages disagreeing.
    let mut theta: Vec<ThetaSeg> = Vec::new();
    for r in curved.rulings.get(&k).map_or(&[][..], Vec::as_slice) {
        let on = [
            line_at(jd, &line_of, r.end[0])?,
            line_at(jd, &line_of, r.end[1])?,
        ];
        let at = combinatorics::exact_point(jd, k, def, r.end[0]).map_err(reject)?;
        if height(on[0], on[1], &at)? != Ordering::Less {
            return Err(reject(RejectReason::CylinderStagesDisagree));
        }
        theta.push(ThetaSeg {
            wall: r.wall,
            side: r.side,
            end: r.end,
            on,
            label: r.label,
            outer: r.outer,
            #[cfg(test)]
            marks: r.marks.clone(),
            station: None,
        });
    }
    let theta = stations_of_lines(jd, cyls, k, theta, &height)?;
    let (events, rim_events) = events_of(jd, k, def, &theta, curved, &lines, &height)?;
    // Each slab's order, asked just after its event — or, with no event, of the whole circle.
    let order = |at: Option<&combinatorics::ExactPoint>| -> Result<Vec<usize>, BoolError> {
        let mut v: Vec<usize> = (0..lines.len()).collect();
        combinatorics::try_sort_by(&mut v, |&a, &b| {
            match combinatorics::line_order(jd, def, lines[a], lines[b], at, Toward::After)? {
                // Two lines level all along a slab are one plane, which the classes are not.
                Ordering::Equal if a != b => Err(RejectReason::CylinderStagesDisagree),
                o => Ok(o),
            }
        })
        .map_err(reject)?;
        Ok(v)
    };
    let slabs = if events.is_empty() {
        vec![order(None)?]
    } else {
        events
            .iter()
            .map(|e| order(Some(&e.at)))
            .collect::<Result<_, _>>()?
    };
    Ok(Chart {
        lines,
        theta,
        events,
        rim_events,
        slabs,
    })
}

/// **The line a ruling's end stops on** — the plane of its name that is a line of the chart. The
/// name need not carry the wall (a corner the operand named, a fillet's tangent corner, lies on
/// the ruling too), and it carries at most one line unless two lines meet there, where either is
/// the place. A name with no line is width where a plane has no world statement and the stages
/// disagreeing otherwise.
fn line_at(
    jd: &Judge<'_, WorkingPlane>,
    line_of: &impl Fn(usize) -> Option<usize>,
    n: combinatorics::NodeId,
) -> Result<usize, BoolError> {
    let (planes, _, _) = combinatorics::pierce_name(n)
        .ok_or_else(|| reject(RejectReason::CylinderStagesDisagree))?;
    if let Some(l) = planes.iter().find_map(|&c| line_of(c)) {
        return Ok(l);
    }
    Err(reject(
        if planes.iter().any(|&c| jd.world_rat(c).is_none()) {
            RejectReason::WitnessNotRational
        } else {
            RejectReason::CylinderStagesDisagree
        },
    ))
}

/// **The chart's events** — every ruling's station, every node a rim is cut at and every point two
/// lines cross or touch at, grouped by θ ([`combinatorics::theta_groups`]), each with the lines'
/// heights on it. A generator carries one
/// ruling name: two would be two stations on one line, refused as the one point under two names it
/// is (`CoincidentNodes`). A rim's nodes stand at distinct events in the rim's own order — two at
/// one event are that point again.
#[allow(clippy::type_complexity)]
fn events_of(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    theta: &[ThetaSeg],
    curved: &Curved,
    lines: &[usize],
    height: &impl Fn(usize, usize, &combinatorics::ExactPoint) -> Result<Ordering, BoolError>,
) -> Result<(Vec<Event>, std::collections::HashMap<usize, Vec<usize>>), BoolError> {
    // The points asked, each once: a ruling's lower end (a ruling is vertical, so one end states
    // its θ), and every cut rim's nodes.
    let mut nodes: Vec<combinatorics::NodeId> = Vec::new();
    let mut index = |n: combinatorics::NodeId| match nodes.iter().position(|&x| x == n) {
        Some(p) => p,
        None => {
            nodes.push(n);
            nodes.len() - 1
        }
    };
    let ruling_at: Vec<usize> = theta.iter().map(|t| index(t.end[0])).collect();
    let mut rims: Vec<(usize, Vec<usize>)> = curved
        .split_rims
        .iter()
        .filter(|((kk, _), _)| *kk == k)
        .map(|(&(_, c), rim)| (c, rim.nodes.iter().map(|&n| index(n)).collect()))
        .collect();
    rims.sort_unstable_by_key(|r| r.0); // a HashMap's order decides nothing; stated once
    // **Where two lines cross or touch** — the line of their two planes pierces the cylinder there,
    // twice or once ([`nacre_exact::quad::section_meet`], the comparator's own reading).
    for a in 0..lines.len() {
        for b in a + 1..lines.len() {
            use nacre_exact::quad::{QuadRoot, SectionMeet};
            let w = |c: usize| jd.world_rat(c).expect("a line has a world statement");
            let at = |root| combinatorics::NodeId::pierce(lines[a], lines[b], k, root);
            match nacre_exact::quad::section_meet(&w(lines[a]), &w(lines[b]), def) {
                Some(SectionMeet::Cross) => {
                    index(at(QuadRoot::Lo));
                    index(at(QuadRoot::Hi));
                }
                Some(SectionMeet::Touch) => {
                    index(at(QuadRoot::Double));
                }
                Some(SectionMeet::Apart) => {}
                // Two classes are two planes, and a line runs not along the axis.
                Some(SectionMeet::Same) | None => {
                    return Err(reject(RejectReason::CylinderStagesDisagree));
                }
            }
        }
    }
    let groups = combinatorics::theta_groups(jd, k, def, &nodes).map_err(|e| reject(e.reason()))?;
    let mut event_of = vec![0usize; nodes.len()];
    for (e, g) in groups.iter().enumerate() {
        for &i in g {
            event_of[i] = e;
        }
    }
    let mut events: Vec<Event> = groups
        .iter()
        .map(|g| {
            Ok(Event {
                at: combinatorics::exact_point(jd, k, def, nodes[g[0]]).map_err(reject)?,
                rulings: Vec::new(),
                height: Vec::new(),
            })
        })
        .collect::<Result<_, BoolError>>()?;
    for (j, &n) in ruling_at.iter().enumerate() {
        events[event_of[n]].rulings.push(j);
    }
    for ev in &mut events {
        // The heights on this generator: the lines sorted by the comparator asked at it, lines
        // that meet here sharing a rank.
        let mut by_height: Vec<usize> = (0..lines.len()).collect();
        let at = &ev.at;
        combinatorics::try_sort_by(&mut by_height, |&a, &b| height(a, b, at))?;
        ev.height = vec![0; lines.len()];
        for w in 1..by_height.len() {
            let (a, b) = (by_height[w - 1], by_height[w]);
            ev.height[b] = ev.height[a] + usize::from(height(a, b, at)? == Ordering::Less);
        }
        // **One stretch of a generator, one name.** Rulings of several walls may stand on one
        // generator between different lines; two over one stretch are one line under two names
        // (`CoincidentNodes`) — or, one wall's, a piece the split cut twice.
        for (x, &i) in ev.rulings.iter().enumerate() {
            for &j in &ev.rulings[x + 1..] {
                let ([a0, a1], [b0, b1]) = (theta[i].on, theta[j].on);
                let h = &ev.height;
                if h[a0] < h[b1] && h[b0] < h[a1] {
                    let name = |j: usize| -> Result<RulingName, BoolError> {
                        let t = &theta[j];
                        Ok((
                            t.wall,
                            ruling_side_of(jd, k, def, t.wall, t.end[0]).map_err(reject)?,
                        ))
                    };
                    return Err(reject(if name(i)? == name(j)? {
                        RejectReason::CylinderStagesDisagree
                    } else {
                        RejectReason::CoincidentNodes
                    }));
                }
            }
        }
    }
    let mut rim_events = std::collections::HashMap::new();
    for (c, idx) in rims {
        let evs: Vec<usize> = idx.iter().map(|&i| event_of[i]).collect();
        // The split ordered the rim's nodes by the same θ order, the seam first.
        if evs.windows(2).any(|w| w[0] == w[1]) {
            return Err(reject(RejectReason::CoincidentNodes));
        }
        if evs.windows(2).any(|w| w[0] > w[1]) {
            return Err(reject(RejectReason::CylinderStagesDisagree));
        }
        rim_events.insert(c, evs);
    }
    Ok((events, rim_events))
}

/// **One line, one station** — read off the gate's records ([`crate::planes::SharedRuling`] with
/// two stating classes, one per pair of secant walls through the line): a piece of the chart's
/// rulings gathers every partner wall on its line, on the same side of its own class
/// ([`combinatorics::shared_line`]). A partner piece over the same axis interval, end for end,
/// folds into it ([`ChordSays::Piece`]) — the alias seed has given those points one name each, so
/// the pair is found by **name**, and left as two, the events would hold one stretch of a generator
/// under two names, which they refuse (`events_of`). A partner with no piece overlapping the interval states the line in
/// the plane vocabulary there ([`ChordSays::NoPiece`]); one overlapping it without matching cuts
/// the line elsewhere than this class does, which is the stages disagreeing. The walls are then put
/// in chord order ([`Station::chords`]).
///
/// The piece kept for the station is the one of the **smallest** class among those folded, so the
/// station's name ([`Chart::ruling_name`]) does not depend on the order the pieces arrive in.
fn stations_of_lines(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    k: usize,
    theta: Vec<ThetaSeg>,
    height: &impl Fn(usize, usize, &combinatorics::ExactPoint) -> Result<Ordering, BoolError>,
) -> Result<Vec<ThetaSeg>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let cyl = &cyls[k];
    let records: Vec<&crate::planes::SharedRuling> = cyl
        .shared
        .iter()
        .filter(|sr| sr.stated_by.len() == 2)
        .collect();
    // Each record's two sides, read the first time a piece asks — only records a piece of the chart
    // reaches are read, as the ruling split reads only its own class's.
    let mut sides: Vec<Option<(i8, i8)>> = vec![None; records.len()];
    let mut side_of = |r: usize, t: usize| -> Result<i8, BoolError> {
        let sr = records[r];
        let (a, b) = match sides[r] {
            Some(s) => s,
            None => {
                let read = |t: usize| {
                    combinatorics::shared_line(jd, cyls, k, sr, t)
                        .map(|l| l.side)
                        .ok_or_else(|| {
                            reject(combinatorics::undecided_reason(
                                cyls,
                                &[sr.line.0, sr.line.1],
                            ))
                        })
                };
                let s = (read(sr.line.0)?, read(sr.line.1)?);
                sides[r] = Some(s);
                s
            }
        };
        Ok(if t == sr.line.0 { a } else { b })
    };
    let mut gone = vec![false; theta.len()];
    let mut stations: Vec<Option<Station>> = vec![None; theta.len()];
    for i in 0..theta.len() {
        if gone[i] || stations[i].is_some() {
            continue;
        }
        let (wall, side) = (theta[i].wall, theta[i].side);
        // Every partner of this piece's line, with its own side of it and a record to read the
        // line's direction off.
        let mut partners: Vec<(usize, i8, &crate::planes::SharedRuling)> = Vec::new();
        for (r, &sr) in records.iter().enumerate() {
            let other = match sr.line {
                (a, b) if a == wall => b,
                (a, b) if b == wall => a,
                _ => continue,
            };
            if side_of(r, wall)? == side && !partners.iter().any(|p| p.0 == other) {
                partners.push((other, side_of(r, other)?, sr));
            }
        }
        let Some(&(_, _, line)) = partners.first() else {
            continue;
        };
        let mut members = vec![i];
        let mut chords = vec![(
            i,
            Chord {
                wall,
                side,
                says: ChordSays::Piece {
                    label: theta[i].label,
                    outer: theta[i].outer,
                },
            },
        )];
        for &(q, qs, _) in &partners {
            let on_line =
                |j: usize| !gone[j] && j != i && theta[j].wall == q && theta[j].side == qs;
            // Two pieces on one line overlap where each starts below the other's top, on that
            // generator — asked only where no piece matches end for end.
            let overlapping = || -> Result<bool, BoolError> {
                let at =
                    combinatorics::exact_point(jd, k, &cyl.def, theta[i].end[0]).map_err(reject)?;
                let below = |a: usize, b: usize| -> Result<bool, BoolError> {
                    Ok(height(a, b, &at)? == Ordering::Less)
                };
                for j in (0..theta.len()).filter(|&j| on_line(j)) {
                    if below(theta[j].on[0], theta[i].on[1])?
                        && below(theta[i].on[0], theta[j].on[1])?
                    {
                        return Ok(true);
                    }
                }
                Ok(false)
            };
            let says = match (0..theta.len()).find(|&j| on_line(j) && theta[j].end == theta[i].end)
            {
                Some(j) => {
                    members.push(j);
                    chords.push((
                        j,
                        Chord {
                            wall: q,
                            side: qs,
                            says: ChordSays::Piece {
                                label: theta[j].label,
                                outer: theta[j].outer,
                            },
                        },
                    ));
                    continue;
                }
                None if overlapping()? => {
                    return Err(reject(RejectReason::CylinderStagesDisagree));
                }
                None => ChordSays::NoPiece,
            };
            chords.push((
                usize::MAX,
                Chord {
                    wall: q,
                    side: qs,
                    says,
                },
            ));
        }
        // Chord order, by insertion: the comparison can fail, and a failing comparison handed to a
        // library sort is not a total order.
        let r = line.radial(&cyl.def).ok_or_else(undecided)?;
        let mut ordered: Vec<(usize, Chord)> = Vec::with_capacity(chords.len());
        for c in chords {
            let mut at = ordered.len();
            for (p, o) in ordered.iter().enumerate() {
                if chord_before(jd, &cyl.def, &r, c.1.wall, o.1.wall).map_err(reject)? {
                    at = p;
                    break;
                }
            }
            ordered.insert(at, c);
        }
        let keep = *members
            .iter()
            .min_by_key(|&&j| theta[j].wall)
            .expect("the piece itself is a member");
        for &j in &members {
            gone[j] = j != keep;
        }
        stations[keep] = Some(Station {
            chords: ordered.into_iter().map(|(_, c)| c).collect(),
        });
    }
    Ok(theta
        .into_iter()
        .zip(stations)
        .zip(gone)
        .filter(|(_, g)| !g)
        .map(|((mut seg, st), _)| {
            seg.station = st;
            seg
        })
        .collect())
}

/// **Whether plane class `a`'s chord is met before `b`'s** turning from the `+θ` tangent inward at
/// the line whose direction across the cylinder is `r` ([`crate::planes::SharedRuling::radial`]) —
/// [`Station::chords`]'s order. Each chord's direction is its class's normal turned about the axis,
/// pointed into the cylinder (against `r`); `a` comes first iff `(d_a × d_b)·m > 0`. Refused as
/// `WitnessNotRational` where the arithmetic cannot answer, and as `CylinderStagesDisagree` for two
/// chords along one direction — two secant classes through one line are two planes, so that is a
/// broken premise rather than a tie.
fn chord_before(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    r: &[Rat; 3],
    a: usize,
    b: usize,
) -> Result<bool, RejectReason> {
    use nacre_exact::{cross3_rat, dot3_rat};
    let m = def.dir();
    let inward = |c: usize| -> Option<[Rat; 3]> {
        let w = jd.world_rat(c)?;
        let d = cross3_rat(&m, &[w[0], w[1], w[2]])?;
        if dot3_rat(&d, r)? > Rat::from_int(0) {
            Some([
                Rat::from_int(0).checked_sub(d[0])?,
                Rat::from_int(0).checked_sub(d[1])?,
                Rat::from_int(0).checked_sub(d[2])?,
            ])
        } else {
            Some(d)
        }
    };
    let turn = (|| dot3_rat(&cross3_rat(&inward(a)?, &inward(b)?)?, &m))()
        .ok_or(RejectReason::WitnessNotRational)?;
    match turn.cmp(&Rat::from_int(0)) {
        core::cmp::Ordering::Greater => Ok(true),
        core::cmp::Ordering::Less => Ok(false),
        core::cmp::Ordering::Equal => Err(RejectReason::CylinderStagesDisagree),
    }
}

#[cfg(test)]
#[path = "../../tests/cyl_chart.rs"]
mod tests;
