/// The chart ledgers' column sums (`cyl_chart::probe`), printed.
#[test]
#[ignore = "measurement — prints the D-ladder ledger sums; run last, single-threaded"]
fn measure_d_ladder_ledgers() {
    let d1 = crate::cyl_chart::probe::ROWS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    let d2b = crate::cyl_chart::probe::d2b::ROWS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    {
        let g = *crate::combinatorics::hull_probe::ROWS.lock().unwrap();
        eprintln!(
            "HULL rings {} arc_rings {} below {} undecided {} broken {}",
            g.0, g.1, g.2, g.3, g.4
        );
        let t = *crate::combinatorics::hull_probe::TILTED.lock().unwrap();
        eprintln!("HULL irrational_extremum_arcs {t}");
    }
    let s1 = |f: fn(&crate::cyl_chart::probe::Row) -> usize| d1.iter().map(f).sum::<usize>();
    let s = |f: fn(&crate::cyl_chart::probe::d2b::Row) -> usize| d2b.iter().map(f).sum::<usize>();
    eprintln!(
        "ledger D1b: charts {} refused {} cells {}",
        d1.len(),
        d1.iter().filter(|r| r.refused).count(),
        s1(|r| r.cells),
    );
    eprintln!(
        "ledger D2b: rows {} refused_booleans {} cells {} end_swapped {} end_disk {} end_exact {} \
             end_other {} end_nocircle {} other_present {} src2_disagree {} src0_present {} exist_disagree {} \
             read_refused {} exist_marks_false {} emit {} emit_unknown {} nocircle_present {} \
             arcs_read {} exact_run_arcs {} arcs_no_mark {} arcs_multi_mark {} emitted_faces {}",
        d2b.len(),
        d2b.iter().filter(|r| r.emitter_refused).count(),
        s(|r| r.cells),
        s(|r| r.end_swapped),
        s(|r| r.end_disk),
        s(|r| r.end_exact),
        s(|r| r.end_other),
        s(|r| r.end_nocircle),
        s(|r| r.other_present),
        s(|r| r.src2_disagree),
        s(|r| r.src0_present),
        s(|r| r.exist_disagree),
        s(|r| r.read_refused),
        s(|r| r.exist_marks_false),
        s(|r| r.emit),
        s(|r| r.emit_unknown),
        s(|r| r.nocircle_present),
        s(|r| r.arcs_read),
        s(|r| r.exact_run_arcs),
        s(|r| r.arcs_no_mark),
        s(|r| r.arcs_multi_mark),
        s(|r| r.emitted_faces),
    );
    let hits = crate::arrangement::crossing_probe::HITS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    eprintln!(
        "ledger E3: ruling_crossings {} off_max {:e} side_disagree {}",
        hits.len(),
        hits.iter().flat_map(|h| h.off).fold(0.0_f64, f64::max),
        hits.iter().filter(|h| h.side_f64 != h.side).count(),
    );
    // The populations the region emitter moves, attributed to fixtures.
    eprintln!(
        "ledger D5-P4: station_pairs {} station_name_failures {}",
        s(|r| r.station_pairs),
        s(|r| r.station_name_failures),
    );
    eprintln!(
        "ledger D5-1a: end_other {} of which single_cut {} — by cause {:?}",
        s(|r| r.end_other),
        s(|r| r.end_other_single_cut),
        crate::cyl_chart::probe::other::COUNTS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone(),
    );
    {
        // The whole-circle disagreements, one line per (test, cyl, t, end, above, bits)
        // shape with its count — the population 1a leaves under `Other`.
        let whole = crate::cyl_chart::probe::other::WHOLE
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        let mut shapes: Vec<(String, usize)> = Vec::new();
        for w in &whole {
            let key = format!(
                "{} cyl {} t {} end {} above {} bits {:?}",
                w.test, w.cyl, w.t, w.end, w.above, w.bits
            );
            match shapes.iter_mut().find(|(k, _)| *k == key) {
                Some((_, n)) => *n += 1,
                None => shapes.push((key, 1)),
            }
        }
        for (k, n) in shapes {
            eprintln!("ledger D5-1a whole_disagree ×{n}: {k}");
        }
    }
    {
        // The tie probe's rows.
        let ties = crate::combinatorics::tie_probe::ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        let mut hist: Vec<(crate::combinatorics::tie_probe::Tie, usize)> = Vec::new();
        for (_, t) in &ties {
            match hist.iter_mut().find(|(k, _)| k == t) {
                Some((_, c)) => *c += 1,
                None => hist.push((*t, 1)),
            }
        }
        eprintln!(
            "ledger C2-P1: mixed abstentions {} by kind {hist:?}",
            ties.len()
        );
        let dec = crate::assembly::probe::deciding::ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        let mut tried: Vec<(usize, usize)> = Vec::new();
        for (_, t, ..) in dec.iter().filter(|r| r.3) {
            match tried.iter_mut().find(|(k, _)| *k == *t) {
                Some((_, c)) => *c += 1,
                None => tried.push((*t, 1)),
            }
        }
        tried.sort_unstable();
        let mut offered: Vec<(usize, usize)> = Vec::new();
        for (_, _, o, ..) in &dec {
            match offered.iter_mut().find(|(k, _)| *k == *o) {
                Some((_, c)) => *c += 1,
                None => offered.push((*o, 1)),
            }
        }
        offered.sort_unstable();
        eprintln!(
            "ledger C2-P4: deciding calls {} exhausted {} probes_tried histogram {tried:?} \
                 offered histogram {offered:?}",
            dec.iter().filter(|r| r.3).count(),
            dec.iter().filter(|r| !r.3).count()
        );
        for r in dec.iter().filter(|r| !r.3) {
            eprintln!(
                "ledger C2-P4 exhausted: {} offered {} ties {:?}",
                r.0, r.2, r.4
            );
        }
        {
            let rows = crate::combinatorics::witness_probe::NO_CANDIDATE
                .lock()
                .expect("the probe's lock is never held across a panic")
                .clone();
            let ans = *crate::combinatorics::witness_probe::ANSWERED
                .lock()
                .expect("the probe's lock is never held across a panic");
            eprintln!(
                "ledger C3-P1: cut caps with no candidate inside {} — answered by centre {} \
                     axis step {} chord point {}",
                rows.len(),
                ans[0],
                ans[1],
                ans[2]
            );
            for r in &rows {
                eprintln!("ledger C3-P1 no-candidate: {r}");
            }
        }
        eprintln!(
            "ledger C2b-P3: cylinder faces asked {}",
            *crate::combinatorics::cylinder_asks::COUNT
                .lock()
                .expect("the probe's lock is never held across a panic")
        );
        eprintln!(
            "ledger C2-P5: the nesting retry swallowed non-abstention errors {}",
            *crate::combinatorics::swallowed_probe::COUNT
                .lock()
                .expect("the probe's lock is never held across a panic")
        );
    }
    {
        let rows = crate::cyl_chart::probe::regions::ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        let sum = |f: fn(&crate::cyl_chart::probe::regions::Row) -> usize| {
            rows.iter().map(f).sum::<usize>()
        };
        eprintln!(
            "ledger D5: rows {} emitter_refused {} faces {} band_faces {} ring_faces {} \
                 emitted_cells {}",
            rows.len(),
            rows.iter().filter(|r| r.emitter_refused).count(),
            sum(|r| r.faces),
            sum(|r| r.band_faces),
            sum(|r| r.ring_faces),
            sum(|r| r.emitted_cells),
        );
        for r in rows.iter().filter(|r| r.emitter_refused) {
            eprintln!("ledger D5 refused: {} cyl {}", r.test, r.cyl);
        }
    }
}
