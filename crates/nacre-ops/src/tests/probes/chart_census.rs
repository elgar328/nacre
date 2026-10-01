//! ★ **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement::cyl_chart`] mounts it with `#[path]` as `census`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use super::*;
/// **What the chart's line set looks like, beside what the hand-written roads answer today.**
///
/// It ships no capability: it builds the chart's two axes and *measures* them against the roads
/// that answer today.
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
    rims: &crate::draft::HeldRims,
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
        // The cell reader's inputs: the lines by parameter, and whose solid this class
        // is. ★ One class is one solid's surface — `planes` interns cylinders by handle and
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
        // rows' span ends — a *different spelling* of the axis parameter (the tracer's `footprint.span`)
        // than the classes' `axis_param` — so this says every lateral face's rim lands, with exact
        // `Rat` equality, on a ⊥ class the chart collected. Not a tautology: a rim on a class
        // `chart_of` skipped (no rational coefficients) or a span end the two spellings disagree on
        // would both fail here.
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

        // ★★★★★ **Every vertical line carries an answer.** A ruling without
        // one is a chart that can state where a wall crosses but not what changes across it — and
        // the census below would then be comparing a partial chart. ☑ Measured 862/862 across the
        // suite; the world sense declines none of them. ★ Except a **tangent** station:
        // the wall touches the cylinder along it and nothing changes across it — the face
        // ends there. Its answer is that it has none, and `read_cell` skips it.
        for t in &chart.theta {
            assert!(
                t.label.is_some() || t.side == 0,
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
        let Ok(cells) = chart.cells(jd, k, def) else {
            probe::push(probe::Row {
                z_lines: chart.z_lines.len(),
                theta: chart.theta.len(),
                rows: mine.len(),
                refused: true,
                ..probe::Row::empty()
            });
            continue;
        };

        // ── The cells, read off the horizontal lines. ──
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
                probe::cell_ends::push(probe::cell_ends::Row {
                    cells: cells.len(),
                    read_refused: 1,
                    ..probe::cell_ends::Row::default()
                });
                continue;
            }
        };
        let mut ends = probe::cell_ends::Row {
            cells: cells.len(),
            emitter_refused: emission.is_err(),
            end_swapped: chart.end_swapped,
            ..probe::cell_ends::Row::default()
        };
        // The premise of naming a station on a line — every (z-line, station)
        // pair of this chart asked for the station's canonical name there.
        {
            let mut names: Vec<(usize, i8)> = Vec::new();
            for j in 0..chart.theta.len() {
                if let Ok(n) = chart.ruling_name(jd, k, def, j) {
                    if !names.contains(&n) {
                        names.push(n);
                    }
                }
            }
            for l in &chart.z_lines {
                for &(wall, sd) in &names {
                    ends.station_pairs += 1;
                    if crate::arrangement::crossing_on_ruling(jd, def, l.class, wall, k, sd)
                        .is_err()
                    {
                        ends.station_name_failures += 1;
                    }
                }
            }
        }
        // ── The emitter's regions, held against the reads they came from. ──
        // The walk is the emitter's own code, so its faces are not compared with the emission
        // (that would be a tautology); what is asserted is what the *result* must satisfy given
        // the reads: every emitted cell lies in exactly one face, no unkept cell lies in any,
        // and two adjacent emitted cells lie in one face. A refusal here is the emitter's, made
        // on the same input.
        {
            match regions::walk(
                jd, k, def, kind, side, &chart, &lines, &cells, &reads, curved, rims,
            ) {
                Ok(w) => {
                    let mut owner: Vec<Option<usize>> = vec![None; cells.len()];
                    for (f, cs) in w.cells_of.iter().enumerate() {
                        for &ci in cs {
                            assert!(
                                owner[ci].replace(f).is_none(),
                                "a cell lies in two faces: cyl {k}"
                            );
                            assert_eq!(
                                reads[ci].emit,
                                Some(true),
                                "an unkept cell lies in a face: cyl {k}"
                            );
                        }
                    }
                    for (ci, r) in reads.iter().enumerate() {
                        assert_eq!(
                            owner[ci].is_some(),
                            r.emit == Some(true),
                            "an emitted cell lies in no face: cyl {k}"
                        );
                        for &cj in &w.neighbours[ci] {
                            if let (Some(a), Some(b)) = (owner[ci], owner[cj]) {
                                assert_eq!(a, b, "adjacent emitted cells in two faces: cyl {k}");
                            }
                        }
                    }
                    ends.emitted_faces = w.faces.len();
                }
                Err(_) => {
                    assert!(
                        emission.is_err(),
                        "the walk refused a class the emitter built: cyl {k}"
                    );
                }
            }
        }
        for (r, cell) in reads.iter().zip(&cells) {
            for e in &r.ends {
                match e {
                    End::Disk(_) => ends.end_disk += 1,
                    End::Exact(a) => {
                        ends.end_exact += 1;
                        ends.exact_run_arcs += a.len() - 1;
                    }
                    End::Other(_) => {
                        ends.end_other += 1;
                        // With stations placed by name the only `Other` left should be
                        // the single-cut sector (both walls one ruling — `arc_around`'s first
                        // return); anything else is counted apart so it can be named.
                        if cell.walls.is_some_and(|[x, y]| x == y) {
                            ends.end_other_single_cut += 1;
                        }
                    }
                    // An uncovered piece is an absence, not a read of a chamber.
                    End::Uncovered => ends.end_uncovered += 1,
                    End::NoCircle => ends.end_nocircle += 1,
                }
            }
            let sources = r
                .ends
                .iter()
                .filter(|e| matches!(e, End::Disk(_) | End::Exact(_)))
                .count();
            if sources == 0 && r.present {
                ends.src0_present += 1;
            }
            if r.present && r.ends.iter().any(|e| matches!(e, End::Other(_))) {
                ends.other_present += 1;
            }
            // A present cell with a line that carries no circle of this cylinder: the
            // `circle_on_class` premise at one end instead of both (`src0_present`).
            if r.present && r.ends.iter().any(|e| matches!(e, End::NoCircle)) {
                ends.nocircle_present += 1;
            }
            // ★ The one-mark contract on rim reads:
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
                    ends.arcs_read += 1;
                    match lateral {
                        0 => ends.arcs_no_mark += 1,
                        1 => {}
                        _ => ends.arcs_multi_mark += 1,
                    }
                    // ★ **Not «exactly one»:** a ⊥ cap through
                    // a wall boss's notch cuts the circle *inside the lateral's hole*, and the arc
                    // there carries **no** mark — the face is absent along it, and the cell reads
                    // absent. So the contract is «at most one, and none exactly where the cell is
                    // not there».
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
            ends.disagree += usize::from(r.disagree);
            ends.exist_disagree += usize::from(r.exist_disagree);
            // Both ends cut is exactly the panel road's population, where today `face_spans`
            // drops a sector for existence — the one count with a twin on the other side.
            if r.ends.iter().all(|e| matches!(e, End::Exact(_))) && !r.present {
                ends.exist_marks_false += 1;
            }
            match r.emit {
                Some(true) => ends.emit += 1,
                None => ends.emit_unknown += 1,
                Some(false) => {}
            }
        }
        assert_eq!(
            ends.end_disk
                + ends.end_exact
                + ends.end_other
                + ends.end_uncovered
                + ends.end_nocircle,
            2 * cells.len(),
            "every cell has two ends: cyl {k}"
        );
        // ★★★★★ **The horizontal lines always speak for a face that is there** — asserted where
        // the fact is made: `circle_on_class` leaves a circle on every ⊥ class within a lateral
        // face's span (Crosses inside, Grazes at the rims), and `emit_faces` labels every circle
        // cell and arc outside the keep filter — which is why no "which side of the wall" sign is
        // needed to read a cell. ★ **Not a bare zero:** a cut end the reader could not pair with
        // its rim (`End::Other`) would leave *no* end speaking, and `present` would fall back to
        // the row's span — which cannot see a hole. The true proposition is the emitter's: over
        // such a cell it builds nothing, it refuses the class by name.
        assert!(
            ends.src0_present == 0 || emission.is_err(),
            "a cell with a face has no label at either end and the emitter read it: cyl {k}"
        );
        // ★★★★★ **A cell with a face never has an end that says nothing** — promoted
        // from a count the suite measured 60 → **0** once stations were placed by name. What
        // `Other` still means is a whole-circle interval *beyond* a face's span whose cut rim's
        // arcs disagree about the far side (the other solid stands on one side of its own wall
        // there), and such a cell is absent. A present cell reading `Other` would be a chart with
        // a θ boundary it has no line for — named here, at the fact, not read from the other end.
        assert_eq!(
            ends.other_present, 0,
            "a cell with a face has an end that says nothing: cyl {k}"
        );
        // ★ Record-site assertions rather than counts (the suite measures them 0): a
        // reporting test sees only the rows recorded before
        // it, an assertion here sees every chart and names the offending test.
        assert_eq!(
            ends.end_swapped, 0,
            "a ruling arrived with z descending: cyl {k}"
        );
        // ★★★★★ **Over a present cell whose chamber could not be read, the emitter emits
        // nothing** — the true proposition. A present cell's chamber is empty when its speaking
        // sides contradict each other — an arrangement label defect, which the emitter refuses by
        // name rather than read from either side; the guard does not hide the defect, it says no
        // face is built over it. ★ Stated on `emit_unknown` (present cells only), not on
        // `disagree`, which also counts absent cells the emitter rightly ignores — that stays a
        // ledger column.
        assert!(
            ends.emit_unknown == 0 || emission.is_err(),
            "a present cell's chamber could not be read, yet the emitter emitted: cyl {k}"
        );
        // ★ The premise at one end: a present cell never meets a line with no circle of its own
        // cylinder (`circle_on_class` marks every ⊥ class within a face's span).
        assert_eq!(
            ends.nocircle_present, 0,
            "a cell with a face meets a line with no circle of its cylinder: cyl {k}"
        );

        let mut whole_circle = 0usize;
        let mut odd_k = 0usize;
        // ── The vertical answer, read. ──
        //
        // ★★★★★ **The consistency check is «sign-free», and deliberately so.** Asking "which side
        // of the wall is this sector" would need a third sign beside the two the label already
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
            // ★ **A wall's rulings need not come in pairs any more.** A band's lateral
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
                // ★★ **Existence before membership**, on the vertical axis. The
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
                // `paired` guard above), so a *per-ruling*
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
        probe::rulings::push(probe::rulings::Row {
            rulings,
            wall_flips,
            intervals_with_flip,
            closes,
            does_not_close,
            unpaired,
            grazing_rulings,
        });
        // ── The emitter, seen from the census: its faces are the regions of the chart,
        // checked above against the reads they came from (every emitted cell in exactly one
        // face, adjacent emitted cells in one face, no unkept cell in any).
        probe::cell_ends::push(ends);
    }
}
