//! **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement::cyl_chart`] mounts it with `#[path]` as `census`, which is what keeps
//! `super::` here meaning that module — the file's place in the tree does not change what it
//! belongs to.

use super::*;
/// **The chart's lines and cells, checked where they are made.** It ships no capability: it
/// rebuilds every cylinder class's chart beside the emitter, asserts the chart's invariants on
/// it, and records the counts `cyl_chart`'s tests read.
///
/// Called from the one place the plane arrangement's faces and [`Curved`] are both in hand. It
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
        // A chart that cannot be formed is skipped, not recorded: `emit_lateral`
        // refuses on the very same `chart_of`, so this is a refusal the emitter already named, not
        // an outcome the census introduces. Measured **0** across the suite either way.
        let Ok(chart) = chart_of(jd, cyls, k, curved) else {
            continue;
        };
        // How today's road sees the same cylinder: one row per lateral *face*, each with its own
        // clipped span. The chart has one line set for the whole class.
        let mine: Vec<&crate::bands::CylRow> = rows.iter().filter(|r| r.class == k).collect();
        // The boundary lines, for the merge counters; a refusal is a skip, as for the chart.
        let Ok(boundary) = boundary_lines(jd, k, plane_faces, curved, rows) else {
            continue;
        };
        // **The chart's lines cover its own boundary rule, asserted where the fact is made.**
        // The load-bearing half is the rows' ends — read off the model face's boundary (the face
        // row's `ends`) rather than off the classes the chart collects — so this says every lateral
        // face's end lands, with exact equality of its plane, on a line the chart collected. Not a
        // tautology: an end on a class `chart_of` skipped (no rational coefficients) or an end the
        // two readings disagree on would both fail here.
        for w in &boundary {
            assert!(
                chart.line_planes.contains(w),
                "a chart lost a line its own boundary rule names: cyl {k}, plane {w:?}"
            );
        }
        // And every circle of this cylinder a plane face is bounded by is on a line: the faces are
        // the other reading of the traces the lines come from.
        for lf in plane_faces {
            let crate::planes::ClassIx::Plane(c) = lf.surf else {
                continue;
            };
            let bounded = std::iter::once(&lf.outer)
                .chain(&lf.inner)
                .any(|b| matches!(b, crate::draft::Bound::Circle { cyl } if *cyl == k));
            assert!(
                !bounded || chart.lines.contains(&c),
                "a face bounded by a circle of this cylinder is on no line: cyl {k}, class {c}"
            );
        }
        let side = mine.first().map(|r| r.side);
        assert!(
            mine.iter().all(|r| Some(r.side) == side),
            "one cylinder class carries rows of both solids: cyl {k}"
        );

        // **Every vertical line carries an answer.** A ruling without
        // one is a chart that can state where a wall crosses but not what changes across it — and
        // the census below would then be comparing a partial chart. Measured 862/862 across the
        // suite; the world sense declines none of them. Except a **tangent** station:
        // the wall touches the cylinder along it and nothing changes across it — the face
        // ends there. Its answer is that it has none, and `read_cell` skips it. (A station's side
        // closed by a wall that draws no piece there is silent too, `ChordSays::NoPiece` — the
        // piece itself carries its answer, and this loop reads pieces.) The cell beside the
        // line outside the cylinder (`outer`) is formed with it — every half-edge of the class's
        // walk bounds a cell — and `slit_at` reads it where the lateral keeps the outside.
        for t in &chart.theta {
            assert!(
                t.label.is_some() || t.side == 0,
                "a ruling reached the chart with no label: cyl {k}, wall {}",
                t.wall
            );
            assert_eq!(
                t.outer.is_some(),
                t.label.is_some(),
                "a ruling's two cells came apart: cyl {k}, wall {}",
                t.wall
            );
        }

        // **The premise `cells()` is built on, checked where the fact is made.** A ruling's ends
        // are read as lines off their names; the point itself must lie on that line's plane — the
        // name's reading and the geometry's, two spellings that would part silently otherwise.
        for t in &chart.theta {
            for e in 0..2 {
                let on_plane = jd.world_rat(chart.lines[t.on[e]]).is_some_and(|w| {
                    crate::combinatorics::pierce_meet(jd, k, &cyls[k].def, t.end[e])
                        .and_then(crate::combinatorics::PiercePoint::into_exact)
                        .is_some_and(|(l, s)| {
                            nacre_exact::quad::plane_side(&w, &l, &s) == nacre_exact::Orient::Zero
                        })
                });
                assert!(
                    on_plane,
                    "a ruling ends off the line its name puts it on: cyl {k}, wall {}",
                    t.wall
                );
            }
        }

        let def = &cyls[k].def;
        // **A slab's order holds to its next event.** The order asked just after an event is the
        // order just before the next: lines that changed places in between crossed where the chart
        // has no event. (A pair missing both its crossings in one slab, or a touch, changes no
        // order and is not seen here.)
        let n = chart.events.len();
        for u in 0..n {
            let at = &chart.events[(u + 1) % n].at;
            let mut before = chart.slabs[u].clone();
            crate::combinatorics::try_sort_by(&mut before, |&a, &b| {
                crate::combinatorics::line_order(
                    jd,
                    def,
                    chart.lines[a],
                    chart.lines[b],
                    Some(at),
                    nacre_exact::quad::Toward::Before,
                )
            })
            .expect("the chart ordered these lines, so their order can be formed");
            assert_eq!(
                before, chart.slabs[u],
                "a slab's lines change places between two events: cyl {k}, slab {u}"
            );
        }
        // A line that is no circle: an oblique section.
        let oblique = |l: usize| {
            let w = chart.line_planes[l];
            !crate::planes::class_carries_circle(&[w[0], w[1], w[2]], &cyls[k].def.dir())
        };
        let cells = chart.cells();

        // ── The cells, read off the horizontal lines. ──
        //
        // `cyl_rows` refuses a class with no row before this census runs, so a class here always
        // has a side; stated as a check rather than an `unwrap` so the message names the class.
        let Some(side) = side else {
            panic!("a cylinder class reached the chart with no row: cyl {k}");
        };
        // A cell the reader refuses is a refusal of the whole class (the emitter's `?`); the
        // chart is still recorded, with the comparison empty and the refusal counted.
        let reads: Vec<CellRead<'_>> = match cells
            .iter()
            .map(|c| chart.read_cell(jd, k, def, kind, side, c, curved, rows))
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
            ..probe::cell_ends::Row::default()
        };
        // What the record-site assertions below count; no test reads them as columns.
        let (mut src0_present, mut nocircle_present, mut emit_unknown) = (0usize, 0usize, 0usize);
        // ── The emitter's regions, held against the reads they came from. ──
        // The walk is the emitter's own code, so its faces are not compared with the emission
        // (that would be a tautology); what is asserted is what the *result* must satisfy given
        // the reads: every emitted cell lies in exactly one face, no unkept cell lies in any,
        // and two adjacent emitted cells lie in one face. A refusal here is the emitter's, made
        // on the same input.
        {
            match regions::walk(jd, k, def, kind, side, &chart, &cells, &reads, curved) {
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
                    // **A face bounded along an oblique section has its partner across it** — a
                    // face of the section's class bounded by the same cylinder's curve there. An
                    // edge has two faces; a lateral ending along a section with nothing on the
                    // section's class is the two stages disagreeing about where the solid ends.
                    let partner = |c: usize| {
                        plane_faces.iter().any(|lf| {
                            lf.surf == crate::planes::ClassIx::Plane(c)
                                && std::iter::once(&lf.outer).chain(&lf.inner).any(|b| match b {
                                    crate::draft::Bound::Circle { cyl } => *cyl == k,
                                    _ => b.ring().is_some_and(|r| {
                                        r.walls.iter().any(|w| {
                                            matches!(*w, crate::combinatorics::Wall::Arc { cyl, .. } if cyl == k)
                                        })
                                    }),
                                })
                        })
                    };
                    let on_oblique = |c: usize| {
                        chart
                            .lines
                            .iter()
                            .position(|&x| x == c)
                            .is_some_and(oblique)
                    };
                    for f in &w.faces {
                        for b in std::iter::once(&f.outer).chain(&f.inner) {
                            let sections: Vec<usize> = match b {
                                crate::draft::Bound::Rim { plane, .. } => vec![*plane],
                                _ => b.ring().map_or(Vec::new(), |r| {
                                    r.walls
                                        .iter()
                                        .filter_map(|w| match *w {
                                            crate::combinatorics::Wall::Arc { plane, .. } => {
                                                Some(plane)
                                            }
                                            _ => None,
                                        })
                                        .collect()
                                }),
                            };
                            for c in sections.into_iter().filter(|&c| on_oblique(c)) {
                                assert!(
                                    partner(c),
                                    "a lateral face is bounded along an oblique section no face \
                                     of its class is bounded along: cyl {k}, class {c}"
                                );
                            }
                        }
                    }
                }
                Err(_) => {
                    assert!(
                        emission.is_err(),
                        "the walk refused a class the emitter built: cyl {k}"
                    );
                }
            }
        }
        for r in &reads {
            for e in &r.ends {
                match e {
                    End::Disk(_) => ends.end_disk += 1,
                    End::Arc(_) => ends.end_arc += 1,
                    // An uncovered piece is an absence, not a read of a chamber.
                    End::Uncovered => ends.end_uncovered += 1,
                    End::NoCircle => ends.end_nocircle += 1,
                }
            }
            let sources = r
                .ends
                .iter()
                .filter(|e| matches!(e, End::Disk(_) | End::Arc(_)))
                .count();
            if sources == 0 && r.present {
                src0_present += 1;
            }
            // A present cell with a line that carries no circle of this cylinder: the
            // `circle_on_class` premise at one end instead of both (`src0_present`).
            if r.present && r.ends.iter().any(|e| matches!(e, End::NoCircle)) {
                nocircle_present += 1;
            }
            // The one-mark contract on rim reads:
            // every cut end read carries at most one lateral mark of its own solid, and none only
            // where the cell is absent (`Seated` is a
            // planar face's word — `face_spans` skips it for the same reason). Counted here, not
            // in `read_cell`, so the emitter's read and the census's read are not counted twice.
            for e in &r.ends {
                let End::Arc(arc) = e else { continue };
                let lateral = arc
                    .marks
                    .iter()
                    .filter(|(s, kd)| {
                        *s == side && !matches!(kd, crate::arrangement::SegKind::Seated { .. })
                    })
                    .count();
                // **Not «exactly one»:** a ⊥ cap through a wall boss's notch cuts the circle
                // *inside the lateral's hole*, and the arc there carries **no** mark — the face is
                // absent along it, and the cell reads absent. So the contract is «at most one, and
                // none exactly where the cell is not there».
                assert!(
                    lateral <= 1,
                    "a cut end read carries {lateral} lateral marks of its own solid: {arc:?}"
                );
                assert!(
                    lateral == 1 || !r.present,
                    "a cut end read with no lateral mark of its own solid reads present: {arc:?}"
                );
            }
            ends.exist_disagree += usize::from(r.exist_disagree);
            // Both ends cut is exactly the panel road's population, where today `face_spans`
            // drops a sector for existence — the one count with a twin on the other side.
            if r.ends.iter().all(|e| matches!(e, End::Arc(_))) && !r.present {
                ends.exist_marks_false += 1;
            }
            match r.emit {
                Some(true) => ends.emit += 1,
                None => emit_unknown += 1,
                Some(false) => {}
            }
        }
        // **The horizontal lines always speak for a face that is there** — asserted where
        // the fact is made: `circle_on_class` leaves a circle on every ⊥ class within a lateral
        // face's span (Crosses inside, Grazes at the rims), and `emit_faces` labels every circle
        // cell and arc outside the keep filter — which is why no "which side of the wall" sign is
        // needed to read a cell. The true proposition is the emitter's: over a present cell with
        // no end speaking it builds nothing, it refuses the class by name.
        assert!(
            src0_present == 0 || emission.is_err(),
            "a cell with a face has no label at either end and the emitter read it: cyl {k}"
        );
        // **Over a present cell whose chamber could not be read, the emitter emits
        // nothing** — the true proposition. A present cell's chamber is empty when its speaking
        // sides contradict each other — an arrangement label defect, which the emitter refuses by
        // name rather than read from either side; the guard does not hide the defect, it says no
        // face is built over it. Stated on `emit_unknown` (present cells only), not on the
        // read's `disagree`, which also holds on absent cells the emitter rightly ignores.
        assert!(
            emit_unknown == 0 || emission.is_err(),
            "a present cell's chamber could not be read, yet the emitter emitted: cyl {k}"
        );
        // The premise at one end: a present cell never meets a line with no circle of its own
        // cylinder (`circle_on_class` marks every ⊥ class within a face's span).
        assert_eq!(
            nocircle_present, 0,
            "a cell with a face meets a line with no circle of its cylinder: cyl {k}"
        );

        let mut whole_circle = 0usize;
        let mut odd_k = 0usize;
        // ── The vertical answer, read. ──
        //
        // **The consistency check is «sign-free», and deliberately so.** Asking "which side
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
        let mut grazing_rulings = 0usize;
        let mut unpaired = 0usize;
        for i in 0..chart.lines.len().saturating_sub(1) {
            // An interval: a gap bounded by the same two lines all the way round, the band the
            // walk below goes round. A gap two crossing lines bound differently slab to slab is no
            // band.
            let pair = [chart.slabs[0][i], chart.slabs[0][i + 1]];
            if chart.slabs.iter().any(|o| [o[i], o[i + 1]] != pair) {
                continue;
            }
            let alive: Vec<&ThetaSeg> = (0..chart.events.len())
                .filter_map(|e| chart.wall_at(e, pair))
                .map(|j| &chart.theta[j])
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
            // **A wall's rulings need not come in pairs.** A band's lateral reaches both rulings of
            // every through-axis wall, so the walk below crosses each such wall twice; a panel's or
            // a chain's boundary may run along one ruling of a wall while the other is no face's
            // edge at all (a quarter boss on a corner). Around such an interval the walk crosses
            // that wall once and cannot return to its start — not a contradiction, a walk with an
            // open end — so the closure is asserted only where every wall is crossed an even number
            // of times, and the rest are counted (`unpaired`).
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
                // **Existence before membership**, on the vertical axis. The
                // chart holds every wall's ruling, whether or not *this* lateral face reaches it,
                // and a label only ever answers membership. `face_spans` is the one spelling of
                // the rule and it reads exactly this list.
                // **Asking whether the mark list is *empty* would be vacuous**: a `MergedRuling`
                // exists only because a lateral traced it — so the answer would be `0` because
                // nothing is being looked at, not because the population is empty. The signal
                // `face_spans` actually reads is the **kind**: a face running through says
                // `Transversal`, one whose boundary stops at the line says `Graze` and may not
                // reach the interval at all.
                // Measured over the lib suite: 12,247 rulings with a mark that is not `Graze`, **1,552
                // all `Graze`** — not empty.
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
                // **`label_cells`' final verification, on the cylinder chart.** The plane
                // side ends its labelling by checking every edge's flip relation and calling a
                // failure `LabelConflict`; this is that check, stated where no orientation is
                // needed — walk the circle, XOR what each wall changes, and come back to where you
                // started. Measured 201 intervals, **none** failing to close.
                // It is asserted (not merely counted) because a failure would mean the chart's
                // own labels contradict each other, which is a defect on this side and not a
                // disagreement with today's road.
                //
                // **What it cannot see, measured rather than assumed.** A walk that XORs is
                // blind to any error appearing an **even** number of times around the circle — and
                // every interval asserted here carries an even count of rulings per wall (the
                // `paired` guard above), so a *per-ruling*
                // systematic flip cancels itself and passes. Probed both ways:
                // adding one flip per ruling stays green, seeding the accumulator wrong goes red.
                // So this holds the labels against each other; what holds their absolute sense is
                // `ruling_probe::SIDE_CHECK`, which compares against content and is not a walk.
                assert!(
                    acc == [false; 2],
                    "the chart's labels do not close around an interval: cyl {k}, interval {i}"
                );
                closes += 1;
            }
        }

        probe::push(probe::Row {
            lines: chart.lines.len(),
            theta: chart.theta.len(),
            rows: mine.len(),
            events: chart.events.len(),
            cells: cells.len(),
            whole_circle,
            odd_k,
        });
        probe::rulings::push(probe::rulings::Row {
            rulings,
            wall_flips,
            intervals_with_flip,
            closes,
            unpaired,
            grazing_rulings,
        });
        // ── The emitter, seen from the census: its faces are the regions of the chart,
        // checked above against the reads they came from (every emitted cell in exactly one
        // face, adjacent emitted cells in one face, no unkept cell in any).
        probe::cell_ends::push(ends);
    }
}

/// **The plane classes where this cylinder's rim is a boundary the arrangement already made**:
/// (i) the circles a plane face emitted as a `Bound::Circle`, and (ii) the circles the arrangement
/// split, which emit arcs instead — the census's list of where a lateral face may end.
pub(crate) fn rim_classes(k: usize, plane_faces: &[LocalFace], curved: &Curved) -> Vec<usize> {
    let mut classes: Vec<usize> = Vec::new();
    let mut push = |c: usize| {
        if !classes.contains(&c) {
            classes.push(c);
        }
    };
    for lf in plane_faces {
        let crate::planes::ClassIx::Plane(c) = lf.surf else {
            continue;
        };
        if std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, crate::draft::Bound::Circle { cyl } if *cyl == k))
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
/// [`rim_classes`]) and every row's own ends, by their planes. A band-shaped emission may continue
/// across any other line, so [`emit_lateral`] merges there and nowhere else, and the census asserts
/// the chart's lines' planes ⊇ this list where it is made.
pub(crate) fn boundary_lines(
    jd: &Judge<'_, WorkingPlane>,
    k: usize,
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
) -> Result<Vec<[Rat; 4]>, BoolError> {
    let mut out: Vec<[Rat; 4]> = rim_classes(k, plane_faces, curved)
        .into_iter()
        .map(|c| {
            jd.world_rat(c)
                .ok_or_else(|| crate::reject(crate::RejectReason::WitnessNotRational))
        })
        .collect::<Result<_, _>>()?;
    for r in rows.iter().filter(|r| r.class == k) {
        out.extend(r.ends);
    }
    out.sort();
    out.dedup();
    Ok(out)
}
