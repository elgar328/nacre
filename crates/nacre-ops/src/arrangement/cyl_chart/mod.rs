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
//! **Same stages, not same code.** The chart is an annulus and rectilinear where the plane
//! side is an unbounded DCEL ordered by angle, so it runs the stages on its own grid rather
//! than calling the plane side's — `regions` states what that costs and why.
//!
//! ## Why the lines are orthogonal — a derivation, not a measurement
//!
//! The population gate (`planes`) sorts every `(plane class, cylinder class)` pair into
//! five, and only two leave a mark on the lateral:
//!
//! | plane vs axis | distance | on the lateral | suite |
//! |---|---|---|---|
//! | `n ∥ m` | — | a **circle** — the chart's *horizontal* | 1128 |
//! | `n ⊥ m` | `0 ≤ d < r` | **two rulings** — the chart's *vertical* | 75 |
//! | `n ⊥ m` | `d = r` | **nothing** from another solid's plane (a tangency grazes and divides nothing); the face's **own** tangent wall is a **station with no crossing** — a ruling of `side = 0`, where the face ends | — |
//! | `n ⊥ m` | `d > r` | nothing | 1494 |
//! | oblique | — | **nothing**: the class's faces stand apart from every lateral face and its section — a whole ellipse or arcs of one — is drawn on its own arrangement (`CylRecord::Sections`); otherwise refused, `ObliqueCylinderCut` | 2 |
//!
//! So the chart carries **axis-parallel segments only**. Not a grid, though: a ruling exists over a
//! finite axis interval and a cut circle is arcs, so the cells are what rectilinear *segments* cut
//! out.
//!
//! The gate is the reason, so when the gate changes this does too; the counts above were
//! measured without the offset and tangent walls. The **offset** wall (`0 < d < r`) adds
//! vertical lines at *irrational* θ; the chart never reads a θ **value**, only an order, so
//! nothing breaks. The **tangent** wall (`d = r`)
//! adds nothing at all: a grazing line stations no sector, so the gate records no
//! crossing and this table's third row is «nothing», not «one vertical». The row's
//! other half: a fillet's or a slot's wall is tangent to its **own** cylinder, and that ruling is
//! not a crossing but the face's edge — a station of `side = 0` carrying no label (the seated
//! face states the piece itself, `SegKind::Tangent`), which the cell reader takes as «the face
//! ends here» and the vertical-answer assertion exempts. And a third reading: when the
//! *other* solid's plane runs through the axis, one of its two verticals **is** that tangent
//! ruling — a line two planes and the cylinder share. Its ends have one name (the alias table's,
//! seeded from the operands — `Curved::aliases` is what this chart's stations are canonicalized
//! through), and the line is stated in the plane vocabulary by both faces that touch it — the
//! tangent wall's run and the lateral's own graze — so the ruling vocabulary states nothing there.
//! The same holds on the *other* solid's tangent wall where the face's own edge lies on its
//! tangent ruling (a half cylinder's flat edge): the flat and the lateral each state their side of
//! the plane-pair line on that wall's class (`rulings::lateral_half_of_plane_pair_lines`), and the
//! chart takes nothing from it.
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

use crate::arrangement::{ArcLabel, Curved, Label, RulingExtent, SegKind};
use crate::draft::{Bound, LocalFace};
use crate::planes::{ClassIx, SolidSide, WorkingCyl, WorkingPlane, class_carries_circle};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, combinatorics, reject};
use nacre_exact::Rat;

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

/// One **horizontal** line: a ⊥ plane class, and where it crosses the axis.
///
/// `t` is exact (`Rat`) — a ⊥ class's axis parameter is a rational, and `bands::axis_param` is the
/// one spelling of it. A class whose parameter is too wide for `Rat` refuses rather than being
/// skipped: a missing band boundary merges two regions whose membership differs, which is a closed
/// and wrong answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ZLine {
    pub(crate) class: usize,
    pub(crate) t: Rat,
}

/// One **vertical** segment: a ruling, and the axis interval it spans.
///
/// Carried straight from [`RulingExtent`] — the plane pass decides this and the chart must not
/// decide it again.
#[derive(Clone, Debug)]
pub(crate) struct ThetaSeg {
    pub(crate) wall: usize,
    /// The ruling's side on its wall — `0` a **tangent** station: a line the face ends
    /// on, with no chamber to read behind it, so no `label`.
    pub(crate) side: i8,
    /// The piece's own two pierce nodes: the θ **order** is asked of these
    /// (`circular_order_about_seam` via `combinatorics::circular_order`), never of a coordinate.
    pub(crate) end: [combinatorics::NodeId; 2],
    /// `[lower, upper]` on the axis.
    pub(crate) z: [Rat; 2],
    /// **The chart's vertical answer** — the label of the cell this ruling borders inside the
    /// cylinder, carried straight from [`RulingExtent::label`]. The horizontal lines' answers are
    /// `DiskLabels`/`ArcLabels`; this is the half they do not carry, and [`Chart::read_cell`] reads
    /// it beside them — the cell's chamber stands only where the two agree.
    pub(crate) label: Option<Label>,
    /// The label of the cell beside this ruling **outside** the cylinder — [`RulingExtent::outer`].
    /// A cell's chamber never reads it; [`Chart::slit_at`] does, where the lateral keeps the
    /// outside.
    pub(crate) outer: Option<Label>,
    /// Who traced this line — [`RulingExtent::marks`]. Membership is `label`'s question; whether
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
    /// plane vocabulary alone (`ruling_sweep`'s `quiet_nodes`). A record is written only for a line
    /// on a lateral *face* (`ruling_on_a_face`), and a tangent class is no secant chord, so this is
    /// the one way a secant wall of a recorded line draws nothing there (measured: every one in the
    /// suite, the census and the ignored sweep). The chord still closes its chamber — the two
    /// planes and the axis say so — but no label of it can be read: the side it bounds is silent,
    /// and the lateral ends at the line there.
    NoPiece,
}

/// One cylinder class's chart: the two axes, and the cells they cut.
#[derive(Clone, Debug, Default)]
pub(crate) struct Chart {
    pub(crate) z_lines: Vec<ZLine>,
    pub(crate) theta: Vec<ThetaSeg>,
    /// How many rulings arrived with `z` descending, so that [`chart_of`] had to swap `end`
    /// alongside `z` to keep `end[e] ↔ z[e]` true. Always 0 — a ruling `Piece`'s `end` is stated in
    /// ascending axis order — and the census asserts it, which is the count's one reader.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) end_swapped: usize,
}

/// A ruling's **name**: the wall class that cut it, and which of the two roots it is.
///
/// One name, not two. `RulingExtent` also carries a `side`, but that is *derived* from the same
/// pair (`combinatorics::ruling_side`); a second spelling of one identity is where two copies
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

/// **One cell of the chart** — an axis interval crossed with a θ-sector.
///
/// Addressed by `(interval, sector)` rather than kept as a flat list, because the region walk
/// ([`regions`]) asks for **neighbours**: crossing a horizontal line moves to the interval above or
/// below, crossing a ruling to the sector beside. That adjacency is *derivable* from this address,
/// and deriving it is the reader's job — building it here would be a field with no consumer.
///
/// The adjacency is not 1:1 and the address already says so: an interval carrying one cell sits
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
    /// **Indices, not names.** Kept as `(wall, side)` names, the ruling's *label* could not be
    /// reached from a cell without searching for it. The name is one call away
    /// ([`Chart::ruling_name`]).
    pub(crate) walls: Option<[usize; 2]>,
}

impl Chart {
    /// **The cells the two axes cut** — `Err` names why a θ order could not be formed.
    ///
    /// **Every ruling's axis endpoints are themselves z-lines** — measured across the suite
    /// (263 charts, 0 exceptions) and structural besides: a ruling's node is named
    /// `[wall, other plane]`, and a plane that meets the lateral while the wall holds the axis can
    /// only be ⊥ to it, which is exactly what `chart_of` collects. So the arrangement **splits at
    /// every axis interval**, and inside one interval a cell is just a θ-sector: no face walk, no
    /// angular order at a vertex.
    ///
    /// **A refusal loses the whole chart, not one interval.** Skipping the interval that
    /// refused would leave a *partial* chart that the census then compares as if it were whole —
    /// the exact shape in which a counting claim goes blind. Measured 0 refusals; the shape is
    /// not left to that.
    pub(crate) fn cells(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        k: usize,
        def: &nacre_topo::CylinderDef,
    ) -> Result<Vec<Cell>, RejectReason> {
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
            // **One node per ruling.** A ruling is vertical, so its two ends share a θ —
            // handing both to `circular_order` asks it to rank one point twice, and it honestly
            // refuses that as the two-names-for-one-point it checks for. The order is a property
            // of the *segment*, and one endpoint states it.
            let nodes: Vec<combinatorics::NodeId> =
                alive.iter().map(|&j| self.theta[j].end[0]).collect();
            let order = crate::combinatorics::circular_order(jd, k, def, &nodes)
                .map_err(crate::combinatorics::CircleOrderFail::reason)?;
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
        Ok(out)
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

    /// **A station's canonical name on a z-line** — `crossing_on_ruling(c, wall, k, side)`: the
    /// very name `ruling_sweep` gives a piece ending on that line and `split_circles` gives a rim
    /// node there (`NodeId::pierce` folds the pair and the root), so «does this rim carry this
    /// station» is **name equality**, and a station the rim does not carry joins the θ order under
    /// the name a piece ending there would have worn.
    ///
    /// **Not a node of another line**. Placing such a station by
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
    ) -> Result<combinatorics::NodeId, RejectReason> {
        let (wall, side) = self.ruling_name(jd, k, def, i)?;
        // As the table knows it: a station on a corner is the corner. The wall's ruling reaches
        // this line (the cell's wall spans the interval), so no crossing is the stages
        // disagreeing ([`crate::combinatorics::RulingNameFail::reason`]).
        crate::combinatorics::crossing_on_ruling(jd, def, c, wall, k, side)
            .map(|n| aliases.canon_point(n))
            .map_err(crate::combinatorics::RulingNameFail::reason)
    }
}

/// **The plane classes where this cylinder's rim is a boundary the arrangement already made**:
/// (i) the circles a plane face emitted as a `Bound::Circle`, and (ii) the circles the arrangement
/// split, which emit arcs instead. The chart's line set, its boundary rule and its census all read
/// this one list.
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
        .split_rims
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

/// **The lines a lateral face's band may not run across** — the rim classes (both kinds,
/// [`rim_classes`]) and every row's own span ends, as axis parameters. A band-shaped emission may
/// continue across any other line, so [`emit_lateral`] merges there and nowhere else, and the
/// census asserts the chart's z-lines ⊇ this list where it is made.
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
/// **No clip to a face.** Clipping the lines to a face's span would answer about a face; a class
/// has no such span — one lateral surface can carry several faces with gaps between them, which is
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
    plane_faces: &[LocalFace],
    curved: &Curved,
) -> Result<Chart, BoolError> {
    let def = &cyls[k].def;
    // The rim classes are read without asking whether they are horizontals: one that is not (a cut
    // rim of an oblique section, which the split refuses before it exists) would be filed at its
    // axis meet here, and `Lines::of`, which every emit reads after this chart, refuses it by name.
    let mut classes: Vec<usize> = rim_classes(k, plane_faces, curved);
    let push = |c: usize, v: &mut Vec<usize>| {
        if !v.contains(&c) {
            v.push(c);
        }
    };
    // (iii) every ⊥ class, full stop — not only those landing on a face's rim, which would answer
    // about a face; a chart covers the whole axis.
    for c in 0..jd.planes.len() {
        let Some(coeffs) = jd.world_rat(c) else {
            continue;
        };
        if class_carries_circle(&[coeffs[0], coeffs[1], coeffs[2]], &def.dir()) {
            push(c, &mut classes);
        }
    }
    let mut z_lines: Vec<ZLine> = classes
        .into_iter()
        .map(|class| crate::bands::axis_param(jd, class, def).map(|t| ZLine { class, t }))
        .collect::<Result<_, BoolError>>()?;
    z_lines.sort_by_key(|l| l.t);
    // Deduped by **parameter**, and that is exact rather than lucky: two classes ⊥ to this axis
    // at one `t` are the *same plane*, which is why span ends match against `t` without a
    // tolerance (`boundary_lines`). So one of the two indices survives arbitrarily — which is
    // why everything downstream compares `t`, never the class index.
    z_lines.dedup_by_key(|l| l.t);

    // **`end` travels with `z`.** `per_class` states `z[e]` as the axis parameter of `end[e]`,
    // so sorting `z` alone would silently break that pairing — and the cell reader asks
    // "this ruling's node *on this z-line*" of exactly that pairing. Swapping both keeps it a
    // carried fact rather than a rule spelled a second time (`pierce_axis_param` is not re-asked).
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
                        outer: r.outer,
                        #[cfg(test)]
                        marks: r.marks.clone(),
                        station: None,
                    }
                })
                .collect()
        })
        .unwrap_or_default();
    let theta = stations_of_lines(jd, cyls, k, theta)?;
    Ok(Chart {
        z_lines,
        theta,
        end_swapped,
    })
}

/// **One line, one station** — read off the gate's records ([`crate::planes::SharedRuling`] with
/// two stating classes, one per pair of secant walls through the line): a piece of the chart's
/// rulings gathers every partner wall on its line, on the same side of its own class
/// ([`combinatorics::shared_line`]). A partner piece over the same axis interval, end for end,
/// folds into it ([`ChordSays::Piece`]) — the alias seed has given those points one name each, so
/// the pair is found by **name**, and left as two, the cell builder would hand `circular_order` one
/// point under two stations. A partner with no piece overlapping the interval states the line in
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
            let says = match (0..theta.len())
                .find(|&j| on_line(j) && theta[j].end == theta[i].end && theta[j].z == theta[i].z)
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
                None if (0..theta.len()).any(|j| {
                    on_line(j)
                        && theta[j].z[0].max(theta[i].z[0]) < theta[j].z[1].min(theta[i].z[1])
                }) =>
                {
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

/// **Is class `c` a horizontal of cylinder `def`'s chart** — a class across the axis, whose trace on
/// the lateral is a circle at one station of the axis ([`class_carries_circle`])? A class oblique to
/// the axis is not: its section is no line of the chart. A class with no world statement reads as
/// one; the chart asks its station by name and refuses it there (`bands::axis_param`). The chart's
/// backstop for a label filed off its horizontals ([`Lines::of`]); the per-class stage asserts it
/// agrees with the gate's record.
pub(crate) fn is_horizontal(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    def: &nacre_topo::CylinderDef,
) -> bool {
    jd.world_rat(c)
        .is_none_or(|w| class_carries_circle(&[w[0], w[1], w[2]], &def.dir()))
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
            .chain(curved.split_rims.keys());
        for &(kk, c) in keys {
            if kk != k {
                continue;
            }
            // Not swallowed: `axis_param` refuses through `reject()`, which a test build's
            // reject census records at the raise whether or not the caller looks. `chart_of`
            // (iii) already asked this of every ⊥ class, so a refusal here cannot be reached after
            // a chart was built — the `?` removes a silent arm rather than adding a population.
            // An oblique class's axis meet is no station of the chart, and filed here it would
            // relabel the class whose station it equals. The per-class stage keeps a recorded
            // section's labels out by the gate's record; this backstop asks the geometry
            // ([`is_horizontal`]) of every key, the cut-rim keys too, which no stage filters — so
            // the two readings parting is a refusal by name, never a label filed at a false
            // station.
            if !is_horizontal(jd, c, def) {
                debug_assert!(
                    false,
                    "a label of cylinder {k} filed under class {c}, which is not across its axis"
                );
                return Err(reject(RejectReason::ObliqueCircleClass));
            }
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
#[path = "../../tests/cyl_chart.rs"]
mod tests;
