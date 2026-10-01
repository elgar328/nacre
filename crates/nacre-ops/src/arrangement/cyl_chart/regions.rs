//! **A lateral's result faces are the regions of its chart.**
//!
//! The plane side's stages, on the chart: cells → labels (read off the neighbouring classes)
//! → emitted cells → **connected components** → the boundary of each component → the boundary's
//! runs cut into the **neighbouring class's own pieces** (a rim's arcs between its cut nodes, a
//! wall's ruling pieces) → cycles → a `Bound`. A band, a panel, a chain rim, a hole and a
//! notched panel are one thing here: a component and its boundary cycles. There is no band /
//! panel dispatch, no both-rims-cut rule and no run ladder — a per-interval vocabulary
//! cannot spell a region that crosses a z-line
//! transversally in one sector and ends on it in another: the
//! corner it asks for is a node no class has.
//!
//! ★ **«Same engine» is the same stages, not the same code.** The planar DCEL's `walk_cells`
//! orders half-edges by angle at a vertex, `nest_cells` asks `point_in_ring` on the plane and
//! `label_cells` propagates from an unbounded root; the chart is an annulus (no unbounded
//! region, wrapping cycles) and rectilinear (four ways at a corner, cells known by
//! construction), and its labels are read rather than propagated — so it runs the stages on
//! its own grid.
//!
//! ★ **Vertices are where the neighbouring class has them.** A boundary run along a rim is cut
//! at the rim's nodes as the cleaned cap faces hold them ([`crate::draft::HeldRims`] — not the
//! split, whose nodes the cleaning pass may have dissolved) and a run along a ruling at the wall
//! class's piece ends (`Curved.rulings`); a station the boundary passes straight through gets no
//! vertex unless that class split its edge there. That is what keeps the lateral's edges welded
//! to the faces beside them, and what the refused corner violated.

use super::{Cell, CellRead, Chart, Lines};
use crate::arrangement::Curved;
use crate::combinatorics::NodeId;
use crate::combinatorics::Wall;
use crate::draft::{Bound, LocalFace, Ring};
use crate::planes::{ClassIx, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};

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

/// The walk's refusal where the split's own record is broken — the trace and the split disagree
/// about where a boundary runs: a run along a rim whose end is not one of the rim's nodes, a run
/// along a ruling the wall's pieces do not tile, a cycle whose pieces do not chain, a component
/// with no outer cycle, a run along an uncut rim that is not the whole circle. What the walk
/// cannot spell is named where it fires: a cut rim of one node and the cycles' winding not
/// classifying (`ArcBoundNotYet` — the assembly's own limits on chains and holes, which
/// `classify_cycles` states).
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
    lines: &Lines,
    cells: &[Cell],
    reads: &[CellRead<'_>],
    curved: &Curved,
    rims: &crate::draft::HeldRims,
) -> Result<Walked, BoolError> {
    let aliases = &curved.aliases;

    // ── 1. Stations: one per (wall, side), in one global θ order. ──
    let mut names: Vec<(usize, i8)> = Vec::new();
    let mut reps: Vec<NodeId> = Vec::new();
    let mut name_of_theta: Vec<usize> = Vec::with_capacity(chart.theta.len());
    for j in 0..chart.theta.len() {
        let n = chart.ruling_name(jd, k, def, j).map_err(reject)?;
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
        let (order, _) = crate::arrangement::circular_order(jd, k, def, &reps)
            .map_err(|e| reject(e.reason()))?;
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

    // ★ **Where the lateral ends although both sides keep it** — a station two walls state
    // ([`Chart::slit_at`]): the result's material between the two walls is gone, so what is
    // kept on the two sides meets on the line alone. The lateral is cut there — no adjacency
    // across, and a boundary on both sides even within one component — so the line becomes an
    // edge the other solid's two faces use too, and the assembly sees the touch: four uses of
    // one edge in one solid, or two solids apart. Run through as one face, the lateral hides it
    // under a shell every count calls closed.
    let slits: Vec<bool> = (0..chart.theta.len())
        .map(|j| chart.slit_at(jd, k, def, j, side, kind).map_err(reject))
        .collect::<Result<_, BoolError>>()?;
    let slit = |j: usize| slits.get(j).copied().unwrap_or(false);

    // ── 2. Emitted cells and their components (4-adjacency, θ-periodic). ──
    let emitted: Vec<bool> = reads.iter().map(|r| r.emit == Some(true)).collect();
    let neighbours: Vec<Vec<usize>> = (0..cells.len())
        .map(|ci| {
            let c = &cells[ci];
            let mut out = Vec::new();
            let row = &by_int[c.interval];
            if row.len() > 1 {
                let s = c.sector;
                let [x, y] = c.walls.unwrap_or([usize::MAX; 2]);
                if !slit(y) {
                    out.push(row[(s + 1) % row.len()]);
                }
                if !slit(x) {
                    out.push(row[(s + row.len() - 1) % row.len()]);
                }
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
            .find(|&c| curved.split_rims.contains_key(&(k, c)))
            .or_else(|| {
                cs.iter()
                    .copied()
                    .find(|&c| curved.disk_labels.contains_key(&(k, c)))
            })
    };
    // A station's canonical name on a line — [`Chart::station_on`]'s spelling, by name.
    // A rim station the run ends at is one the split placed, so a crossing that is not there is
    // the stages disagreeing ([`crate::arrangement::RulingNameFail::reason`]).
    let station_name = |c: usize, s: usize| -> Result<NodeId, BoolError> {
        let (wall, sd) = name_at_station(s);
        crate::arrangement::crossing_on_ruling(jd, def, c, wall, k, sd)
            .map(|n| aliases.canon_point(n))
            .map_err(|e| reject(e.reason()))
    };

    // ── 3–7. Per component: boundary edges → cycles → runs → pieces → rings → bound. ──
    let mut faces: Vec<LocalFace> = Vec::new();
    let mut cells_of: Vec<Vec<usize>> = Vec::new();
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
                    let [x, y] = c.walls.unwrap_or([usize::MAX; 2]);
                    if !in_c(right) || slit(y) {
                        edges.push(Edge::Ruling {
                            interval: c.interval,
                            station: b,
                            up: true,
                        });
                    }
                    if !in_c(left) || slit(x) {
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
            let mut prev = dir(edges[start]);
            while at != first {
                let mut cands: Vec<(usize, usize)> = (0..edges.len())
                    .filter(|&j| !used[j] && ends(edges[j]).0 == at)
                    .map(|j| {
                        let rel = (dir(edges[j]) + 4 - prev) % 4;
                        let rank = match rel {
                            1 => 0, // left
                            0 => 1, // straight
                            3 => 2, // right
                            _ => 3, // back
                        };
                        (rank, j)
                    })
                    .collect();
                cands.sort_unstable();
                let Some(&(_, j)) = cands.first() else {
                    return Err(split_disagrees());
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
                        let Some(rim) = rims.get(&(k, c)) else {
                            if !full {
                                return Err(split_disagrees());
                            }
                            whole_rim = Some((c, ccw));
                            continue;
                        };
                        // ★ A cut rim's pieces are the arcs between its consecutive nodes
                        // — the held table, which the assembly's arc join reads too; neither the
                        // split nor the arc labels is consulted for the edges.
                        let m = rim.nodes.len();
                        // One node cannot state a rim as arcs, and the cut-rim table keeps such a rim
                        // on purpose (`draft::held_rims`). No node at all is a rim the split never cut.
                        if m < 2 {
                            return Err(if m == 1 {
                                reject(RejectReason::ArcBoundNotYet)
                            } else {
                                split_disagrees()
                            });
                        }
                        let (pa, pb) = if full {
                            (0, 0)
                        } else {
                            let na = station_name(c, ca.1)?;
                            let nb = station_name(c, cb.1)?;
                            let pa = rim
                                .nodes
                                .iter()
                                .position(|&n| n == na)
                                .ok_or_else(split_disagrees)?;
                            let pb = rim
                                .nodes
                                .iter()
                                .position(|&n| n == nb)
                                .ok_or_else(split_disagrees)?;
                            (pa, pb)
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
                            return Err(split_disagrees());
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
            // The outer cycle of a component with no wrapping rim is the one holding the
            // lowest interval's bottom side — always a boundary, never a hole's.
            let holds_lowest_bottom = cyc
                .iter()
                .any(|e| matches!(e, Edge::Arc { line, ccw: true, .. } if *line == lowest));
            if holds_lowest_bottom && outer.is_none() {
                outer = Some(rings.len());
            }
            rings.push(Ring::new(nodes, walls));
        }
        let (own, other) = reads[members[0]]
            .chamber
            .expect("an emitted cell has a chamber");
        let flip = !super::read_cell::keep_for(kind, side, own, other);
        let face = if rims_lo.is_empty() && rims_hi.is_empty() {
            let o = outer.ok_or_else(split_disagrees)?;
            let mut rings = rings;
            let outer = rings.remove(o);
            LocalFace {
                surf: ClassIx::Cyl(k),
                outer: Bound::Ring(outer),
                inner: rings.into_iter().map(Bound::Ring).collect(),
                flip,
            }
        } else {
            crate::assembly::classify_cycles(k, flip, rings, rims_lo, rims_hi, rims)
                .map_err(|_| reject(RejectReason::ArcBoundNotYet))?
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
