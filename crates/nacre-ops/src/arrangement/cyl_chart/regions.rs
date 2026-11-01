//! **A lateral's result faces are the regions of its chart.**
//!
//! The plane side's stages, on the chart: cells → labels (read off the neighbouring classes)
//! → emitted cells → **connected components** → the boundary of each component → the boundary's
//! runs cut into the **neighbouring class's own pieces** (a rim's arcs between its cut nodes, a
//! wall's ruling pieces) → cycles → a `Bound`. A band, a panel, a chain rim, a hole and a
//! notched panel are one thing here: a component and its boundary cycles. There is no band /
//! panel dispatch, no both-rims-cut rule and no run ladder — a per-interval vocabulary
//! cannot spell a region that crosses a line transversally in one slab and ends on it in another:
//! the corner it asks for is a node no class has.
//!
//! **«Same engine» is the same stages, not the same code.** The planar DCEL's `walk_cells`
//! orders half-edges by angle at a vertex, `nest_cells` asks `point_in_ring` on the plane and
//! `label_cells` propagates from an unbounded root; the chart is an annulus (no unbounded
//! region, wrapping cycles) whose cells are known by construction and whose ways out of a corner
//! are ordered by the slabs' height orders it already holds (`departure_key`), and its labels are
//! read rather than propagated — so it runs the stages on its own cells.
//!
//! **Vertices are where the neighbouring class has them.** A boundary run along a rim is cut
//! at the rim's nodes as the split placed them (`Curved::split_rims`) and a run along a ruling at
//! the wall class's piece ends (`Curved.rulings`) — the nodes the neighbouring faces were cut at
//! too; a station the boundary passes straight through gets no vertex unless that class split its
//! edge there. The result's one cleaning pass then drops a node no ring turns at from every ring
//! at once (`arrangement::whole_result_faces`), the lateral's with its neighbour's. That is what
//! keeps the lateral's edges welded to the faces beside them, and what the refused corner
//! violated.

use super::{Cell, CellRead, Chart};
use crate::arrangement::Curved;
use crate::combinatorics::NodeId;
use crate::combinatorics::Wall;
use crate::draft::{Bound, LocalFace, Ring};
use crate::planes::{ClassIx, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};

/// A corner of the chart: `(event, height)` — a line's height on the event's generator
/// ([`super::Event::height`]), so lines that meet there share it.
type Corner = (usize, usize);

/// One directed boundary edge of the chart, with the region on its left.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Edge {
    /// Along line `line`, over slab `unit` (from event `unit` to `unit + 1`); `ccw` when the
    /// region lies above the line.
    Arc { line: usize, unit: usize, ccw: bool },
    /// Along the ruling at event `station`, between lines `lo` and `hi`; `up` when the region
    /// lies on the station's −θ side (the cell whose slab *ends* at the station).
    Ruling {
        station: usize,
        lo: usize,
        hi: usize,
        up: bool,
    },
}

/// One boundary piece in travel order: a neighbouring class's own edge.
struct Piece {
    wall: Wall,
    ends: (NodeId, NodeId),
}

/// The faces of class `k`, each with the cells it covers, and every cell's neighbours (the
/// adjacency the components were joined on — read back by the census).
pub(crate) struct Walked {
    pub(crate) faces: Vec<LocalFace>,
    /// Read by the census only — an instrument's fields, stated per build rather than
    /// blanket-allowed (the `RulingExtent` precedent).
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) cells_of: Vec<Vec<usize>>,
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) neighbours: Vec<Vec<usize>>,
}

/// **Which way out of a corner the cycle takes first** — the way out nearest clockwise from the way
/// back along the edge it arrived by, which keeps the region on the left and each cycle simple
/// where the region touches itself at a corner. The ways out of a corner, counter-clockwise from
/// `−z`: down the ruling, the lines leaving `+θ` from the lowest to the highest in the slab they
/// enter, up the ruling, the lines leaving `−θ` from the highest to the lowest in theirs — the
/// slabs' height orders, which the chart already holds, are the angles. Going back along the edge
/// arrived by comes last. With one line each way this is left, straight, right, back.
fn departure_key(chart: &Chart, arrived: Edge, out: Edge) -> usize {
    let n = chart.lines.len();
    let total = 2 * n + 2;
    // A way out's place counter-clockwise from `−z`: `0` down, `1..=n` east, `n + 1` up,
    // `n + 2..` west.
    let east = |unit: usize, line: usize| 1 + chart.rank(unit, line);
    let west = |unit: usize, line: usize| n + 2 + (n - 1 - chart.rank(unit, line));
    let place = |e: Edge| match e {
        Edge::Arc {
            line,
            unit,
            ccw: true,
        } => east(unit, line),
        Edge::Arc {
            line,
            unit,
            ccw: false,
        } => west(unit, line),
        Edge::Ruling { up: true, .. } => n + 1,
        Edge::Ruling { up: false, .. } => 0,
    };
    // The way back along the edge arrived by: `+θ` along a line over slab `unit` came in from
    // that slab, so back is `−θ` into it; and the other way round.
    let back = match arrived {
        Edge::Arc {
            line,
            unit,
            ccw: true,
        } => west(unit, line),
        Edge::Arc {
            line,
            unit,
            ccw: false,
        } => east(unit, line),
        Edge::Ruling { up: true, .. } => 0,
        Edge::Ruling { up: false, .. } => n + 1,
    };
    match (back + total - place(out)) % total {
        0 => total,
        k => k,
    }
}

/// The walk's refusal where the split's own record is broken — the trace and the split disagree
/// about where a boundary runs: a run along a rim whose end is not one of the rim's nodes, a run
/// along a ruling the wall's pieces do not tile, a cycle whose pieces do not chain, a component
/// with no outer cycle, a run along an uncut rim — or a rim cut at one node — that is not the whole
/// circle, a cycle whose arc units are not a whole number of turns, the lateral kept on both sides
/// of a station it ends at (`ChordSays::NoPiece`). What the walk cannot spell is
/// named where it fires: wrapping cycles other than one lower and one upper (`ArcBoundNotYet`).
fn split_disagrees() -> BoolError {
    reject(RejectReason::CylinderStagesDisagree)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn walk(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    def: &nacre_topo::CylinderDef,
    kind: crate::BoolKind,
    side: crate::planes::SolidSide,
    chart: &Chart,
    cells: &[Cell],
    reads: &[CellRead<'_>],
    curved: &Curved,
) -> Result<Walked, BoolError> {
    // ── 1. Events: the chart's, in θ order; a cell is a slab crossed with a gap. ──
    let ns = chart.events.len();
    let n_gap = chart.lines.len().saturating_sub(1);
    // Cells by gap, in slab order (as `cells()` pushed them).
    let mut by_gap: Vec<Vec<usize>> = vec![Vec::new(); n_gap];
    for (ci, c) in cells.iter().enumerate() {
        by_gap[c.gap].push(ci);
    }
    let units = ns.max(1);
    // A line's height on event `e`'s generator; with no event, its place in the one slab.
    let height = |e: usize, l: usize| -> usize {
        if ns == 0 {
            chart.rank(0, l)
        } else {
            chart.events[e].height[l]
        }
    };
    // Does the gap between `lines` pinch to a point at event `e` — its two lines meeting there?
    let pinched = |e: usize, [lo, hi]: [usize; 2]| ns > 0 && height(e, lo) == height(e, hi);

    // **Where the lateral ends although both sides keep it** — a station secant walls run
    // through ([`Chart::slit_at`]): on the lateral's material side, the result has nothing
    // between two of the walls, so what is kept on the two sides meets on the line alone. The
    // lateral is cut there — no adjacency across, and a boundary on both sides even within one
    // component — so the line becomes an edge the other solid's two faces use too, and the
    // assembly sees the touch: four uses of one edge in one solid, or two solids apart. Run
    // through as one face, the lateral hides it under a shell every count calls closed.
    let slits: Vec<bool> = (0..chart.theta.len())
        .map(|j| chart.slit_at(jd, j, side, kind).map_err(reject))
        .collect::<Result<_, BoolError>>()?;
    let slit = |j: usize| slits.get(j).copied().unwrap_or(false);

    // ── 2. Emitted cells and their components (4-adjacency, θ-periodic). ──
    let emitted: Vec<bool> = reads.iter().map(|r| r.emit == Some(true)).collect();
    // A station with a wall that draws no piece is the lateral's own edge (`ChordSays::NoPiece`),
    // and `slit_at` reads it as one the lateral does not cross. Cells kept on both its sides in one
    // gap would be the lateral crossing it after all — the premise broken, refused rather
    // than walked through. A cell whose two walls are this one station (the circle opened there
    // alone) lies on both its sides.
    for (j, seg) in chart.theta.iter().enumerate() {
        let no_piece = seg.station.as_ref().is_some_and(|st| {
            st.chords
                .iter()
                .any(|c| matches!(c.says, super::ChordSays::NoPiece))
        });
        if !no_piece {
            continue;
        }
        for row in &by_gap {
            let kept_beside = |e: usize| {
                row.iter()
                    .any(|&ci| emitted[ci] && cells[ci].walls[e] == Some(j))
            };
            if kept_beside(0) && kept_beside(1) {
                return Err(split_disagrees());
            }
        }
    }
    let neighbours: Vec<Vec<usize>> = (0..cells.len())
        .map(|ci| {
            let c = &cells[ci];
            let mut out = Vec::new();
            let row = &by_gap[c.gap];
            if row.len() > 1 {
                let s = c.slab;
                let [x, y] = c.walls;
                // Across an event a gap goes on in the same place of the next slab's order — two
                // lines crossing there trade places, the gaps beside them do not — unless its two
                // lines meet there, where the cells touch at a point.
                if !y.is_some_and(slit) && !pinched((s + 1) % ns, c.lines) {
                    out.push(row[(s + 1) % row.len()]);
                }
                if !x.is_some_and(slit) && !pinched(s, c.lines) {
                    out.push(row[(s + row.len() - 1) % row.len()]);
                }
            }
            for g2 in [c.gap.wrapping_sub(1), c.gap + 1] {
                if g2 < n_gap {
                    out.push(by_gap[g2][c.slab]);
                }
            }
            out
        })
        .collect();
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
            for &y in &neighbours[x] {
                if emitted[y] && comp[y].is_none() {
                    comp[y] = Some(id);
                    stack.push(y);
                }
            }
        }
        members.sort_unstable();
        components.push(members);
    }

    // Endpoints (corners) of a directed edge.
    let ends = |e: Edge| -> (Corner, Corner) {
        match e {
            Edge::Arc { line, unit, ccw } => {
                let next = (unit + 1) % units;
                let (a, b) = ((unit, height(unit, line)), (next, height(next, line)));
                if ccw { (a, b) } else { (b, a) }
            }
            Edge::Ruling {
                station,
                lo,
                hi,
                up,
            } => {
                let (a, b) = (
                    (station, height(station, lo)),
                    (station, height(station, hi)),
                );
                if up { (a, b) } else { (b, a) }
            }
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
    // The class of a line, where this cylinder marks a circle on it — cut or whole.
    let class_on = |line: usize| -> Option<usize> {
        let c = chart.lines[line];
        (curved.split_rims.contains_key(&(k, c)) || curved.disk_labels.contains_key(&(k, c)))
            .then_some(c)
    };
    // A run along rim `c` ends at one of the rim's own nodes: every node it is cut at is an event
    // (`Chart::events`), and the run turns only where something cuts the rim.
    let rim_node_at = |c: usize, s: usize| -> Option<usize> {
        chart.rim_events.get(&c)?.iter().position(|&e| e == s)
    };

    // ── 3–7. Per component: boundary edges → cycles → runs → pieces → rings → bound. ──
    let mut faces: Vec<LocalFace> = Vec::new();
    let mut cells_of: Vec<Vec<usize>> = Vec::new();
    for (id, members) in components.iter().enumerate() {
        let in_c = |ci: usize| comp[ci] == Some(id);
        let mut edges: Vec<Edge> = Vec::new();
        for &ci in members {
            let c = &cells[ci];
            let u = c.slab;
            let below = c.gap.checked_sub(1).map(|g| by_gap[g][u]);
            if !below.is_some_and(in_c) {
                edges.push(Edge::Arc {
                    line: c.lines[0],
                    unit: u,
                    ccw: true,
                });
            }
            let above = (c.gap + 1 < n_gap).then(|| by_gap[c.gap + 1][u]);
            if !above.is_some_and(in_c) {
                edges.push(Edge::Arc {
                    line: c.lines[1],
                    unit: u,
                    ccw: false,
                });
            }
            // With one event the slab runs round to where it started: the cell is its own
            // neighbour on both sides, and no ruling bounds it.
            if ns > 1 {
                let row = &by_gap[c.gap];
                let right = row[(u + 1) % ns];
                let left = row[(u + ns - 1) % ns];
                let [x, y] = c.walls;
                let [lo, hi] = c.lines;
                // **The region ends at an event only on a ruling** — or at the point its gap
                // pinches to, which needs no edge. An event with no ruling across the gap (a rim
                // cut there, a ruling between other lines) carries no line the boundary could run
                // along — the chart has no vertical line for it.
                let ruling =
                    |w: Option<usize>| w.ok_or_else(|| reject(RejectReason::RulingBoundNotYet));
                let east = (u + 1) % ns;
                if !pinched(east, c.lines) && (!in_c(right) || y.is_some_and(slit)) {
                    ruling(y)?;
                    edges.push(Edge::Ruling {
                        station: east,
                        lo,
                        hi,
                        up: true,
                    });
                }
                if !pinched(u, c.lines) && (!in_c(left) || x.is_some_and(slit)) {
                    ruling(x)?;
                    edges.push(Edge::Ruling {
                        station: u,
                        lo,
                        hi,
                        up: false,
                    });
                }
            }
        }
        // Sorted by where each edge lies — a line's place in its slab, a ruling's lower height —
        // so the cycles start where they always have.
        edges.sort_by_key(|e| match *e {
            Edge::Arc { line, unit, ccw } => (0, chart.rank(unit, line), unit, usize::from(ccw)),
            Edge::Ruling {
                station, lo, up, ..
            } => (1, height(station, lo), station, usize::from(up)),
        });
        edges.dedup();

        // ── cycles, by the left-turn rule at a corner with two ways out (a hole touching
        //    the outer boundary at a corner: the sharper turn keeps each cycle simple) ──
        let mut used = vec![false; edges.len()];
        let mut cycles: Vec<Vec<Edge>> = Vec::new();
        for start in 0..edges.len() {
            if used[start] {
                continue;
            }
            let mut cyc = vec![edges[start]];
            used[start] = true;
            let (first, mut at) = ends(edges[start]);
            let mut last = edges[start];
            while at != first {
                let mut cands: Vec<(usize, usize)> = (0..edges.len())
                    .filter(|&j| !used[j] && ends(edges[j]).0 == at)
                    .map(|j| (departure_key(chart, last, edges[j]), j))
                    .collect();
                cands.sort_unstable();
                let Some(&(_, j)) = cands.first() else {
                    return Err(split_disagrees());
                };
                used[j] = true;
                cyc.push(edges[j]);
                last = edges[j];
                at = ends(edges[j]).1;
            }
            cycles.push(cyc);
        }

        // ── runs → pieces → rings ──
        // **Each cycle's winding about the axis, by its arc units.** An arc edge covers one
        // slab, `+1` walked counter-clockwise (`ccw` — the face above it) and `−1` against, and
        // a ruling covers none; a closed cycle returns to its station, so the sum is a whole
        // number of turns — `+1` for the lower boundary of a face that wraps the axis, `−1` for
        // its upper one, `0` for a cycle that does not wrap. No seam is read: the count is the
        // chart's own, wherever `θ = 0` falls.
        let winding = |cyc: &[Edge]| -> Result<i64, BoolError> {
            let units_walked: i64 = cyc
                .iter()
                .map(|e| match *e {
                    Edge::Arc { ccw: true, .. } => 1,
                    Edge::Arc { ccw: false, .. } => -1,
                    Edge::Ruling { .. } => 0,
                })
                .sum();
            let per_turn = units as i64;
            if units_walked % per_turn != 0 {
                return Err(split_disagrees());
            }
            Ok(units_walked / per_turn)
        };
        let mut bounds: Vec<(Bound, i64)> = Vec::new();
        let mut lowest_bottom: Option<usize> = None;
        // The bottom of the component's first cell — the lowest of its slab — is a boundary.
        let bottom = Edge::Arc {
            line: cells[members[0]].lines[0],
            unit: cells[members[0]].slab,
            ccw: true,
        };
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
                        let c = class_on(line).ok_or_else(split_disagrees)?;
                        let full = ns == 0 || run.len() == ns;
                        // An uncut rim is the assembly's closed edge — whole or nothing.
                        let Some(rim) = curved.split_rims.get(&(k, c)) else {
                            if !full {
                                return Err(split_disagrees());
                            }
                            whole_rim = Some((c, ccw));
                            continue;
                        };
                        // A cut rim's pieces are the arcs between its consecutive nodes — the
                        // split's, which the cap's ring holds too; the result's one cleaning pass
                        // drops a node neither ring turns at from both. The arc labels are not
                        // consulted for the edges. Each piece carries the rim's plane `c`, which
                        // the assembly checks against the cap's.
                        let m = rim.nodes.len();
                        // No node at all is a rim the split never cut. A rim cut at **one** node — a
                        // vertex alone on its circle — is one piece, the arc from that node back to
                        // itself, and only a run the whole way round reaches it: the run's two ends
                        // are that node, so a run ending part way round is the stages disagreeing.
                        if m == 0 || (m == 1 && !full) {
                            return Err(split_disagrees());
                        }
                        let (pa, pb) = if full {
                            (0, 0)
                        } else {
                            (
                                rim_node_at(c, ca.0).ok_or_else(split_disagrees)?,
                                rim_node_at(c, cb.0).ok_or_else(split_disagrees)?,
                            )
                        };
                        let mut p = pa;
                        loop {
                            let q = if ccw { (p + 1) % m } else { (p + m - 1) % m };
                            let (from, to) = if ccw {
                                (rim.nodes[p], rim.nodes[q])
                            } else {
                                (rim.nodes[q], rim.nodes[p])
                            };
                            pieces.push(Piece {
                                wall: Wall::Arc {
                                    cyl: k,
                                    ccw,
                                    plane: c,
                                },
                                ends: if ccw { (from, to) } else { (to, from) },
                            });
                            p = q;
                            if p == pb {
                                break;
                            }
                        }
                    }
                    Edge::Ruling { station, up, .. } => {
                        // The run's two heights on the station's generator, lower first.
                        let (h_lo, h_hi) = if up { (ca.1, cb.1) } else { (cb.1, ca.1) };
                        // The station's rulings are the event's; each piece states its own wall.
                        let h = |l: usize| height(station, l);
                        let mut segs: Vec<usize> = chart.events[station]
                            .rulings
                            .iter()
                            .copied()
                            .filter(|&j| {
                                let [a, b] = chart.theta[j].on;
                                h(a) >= h_lo && h(b) <= h_hi
                            })
                            .collect();
                        segs.sort_by_key(|&j| h(chart.theta[j].on[0]));
                        let tiled = !segs.is_empty()
                            && h(chart.theta[segs[0]].on[0]) == h_lo
                            && h(chart.theta[*segs.last().expect("non-empty")].on[1]) == h_hi
                            && segs.windows(2).all(|w| {
                                h(chart.theta[w[0]].on[1]) == h(chart.theta[w[1]].on[0])
                                    && chart.theta[w[0]].end[1] == chart.theta[w[1]].end[0]
                            });
                        if !tiled {
                            return Err(split_disagrees());
                        }
                        let ordered: Vec<usize> = if up {
                            segs
                        } else {
                            segs.into_iter().rev().collect()
                        };
                        for j in ordered {
                            let [e0, e1] = chart.theta[j].end;
                            let (wall, sd) = chart.ruling_name(jd, k, def, j).map_err(reject)?;
                            pieces.push(Piece {
                                wall: Wall::Ruling {
                                    cyl: k,
                                    side: sd,
                                    up,
                                    plane: wall,
                                },
                                ends: if up { (e0, e1) } else { (e1, e0) },
                            });
                        }
                    }
                }
            }
            let w = winding(cyc)?;
            // An uncut rim is a cycle by itself: anything else on its cycle meets the rim, and
            // a rim something meets is cut — the stages disagree, and the other runs would be
            // dropped with the rim kept whole.
            if whole_rim.is_some() && runs.len() > 1 {
                return Err(split_disagrees());
            }
            if let Some((c, ccw)) = whole_rim {
                bounds.push((Bound::Rim { plane: c, ccw }, w));
                continue;
            }
            // Two pieces close a cycle here, straight or not — a slit's two sides are a digon of
            // rulings, which the assembly's shell guard then names — and one piece does only as the
            // arc from a node back to itself. (Not `combinatorics::ring_floor`: its three for a
            // straight ring is the plane's fact, and would call the slit's digon a broken split.)
            let n = pieces.len();
            let self_loop = n == 1
                && matches!(pieces[0].wall, Wall::Arc { .. })
                && pieces[0].ends.0 == pieces[0].ends.1;
            if n < 2 && !self_loop {
                return Err(split_disagrees());
            }
            let mut nodes = Vec::with_capacity(n);
            let mut walls = Vec::with_capacity(n);
            for (p, piece) in pieces.iter().enumerate() {
                if piece.ends.1 != pieces[(p + 1) % n].ends.0 {
                    return Err(split_disagrees());
                }
                nodes.push(piece.ends.0);
                walls.push(piece.wall);
            }
            // The outer cycle of a component that does not wrap the axis is the one holding the
            // bottom of its lowest cell in a slab — below it the slab runs on without the
            // component, so it is never a hole's side.
            let holds_lowest_bottom = cyc.contains(&bottom);
            if holds_lowest_bottom && lowest_bottom.is_none() {
                lowest_bottom = Some(bounds.len());
            }
            bounds.push((Bound::Ring(Ring::new(nodes, walls)), w));
        }
        let (own, other) = reads[members[0]]
            .chamber
            .expect("an emitted cell has a chamber");
        let flip = !super::read_cell::keep_for(kind, side, own, other);
        // **The outer bound: the lower rim of a face that wraps, or the lowest-bottom cycle of
        // one that does not.** A connected region of the band has either no wrapping boundary or
        // exactly one of each sense — a second lower rim would cut the region in two — so
        // anything else is a walk this chart does not arrange.
        let up: Vec<usize> = (0..bounds.len()).filter(|&i| bounds[i].1 == 1).collect();
        let down = bounds.iter().filter(|b| b.1 == -1).count();
        if bounds.iter().any(|b| b.1.abs() > 1) {
            return Err(reject(RejectReason::ArcBoundNotYet));
        }
        let o = match (&up[..], down) {
            ([o], 1) => *o,
            ([], 0) => lowest_bottom.ok_or_else(split_disagrees)?,
            _ => return Err(reject(RejectReason::ArcBoundNotYet)),
        };
        let mut bounds: Vec<Bound> = bounds.into_iter().map(|(b, _)| b).collect();
        let outer = bounds.remove(o);
        let face = LocalFace {
            surf: ClassIx::Cyl(k),
            outer,
            inner: bounds,
            flip,
        };
        faces.push(face);
        cells_of.push(members.clone());
    }
    Ok(Walked {
        faces,
        cells_of,
        neighbours,
    })
}

#[cfg(test)]
#[path = "../../tests/regions.rs"]
mod tests;
