/// Every probe ledger's column sums, printed — one line per probe, labelled by the probe's name.
#[test]
#[ignore = "measurement — prints the probe ledger sums; run last, single-threaded"]
fn measure_probe_ledgers() {
    let chart = crate::arrangement::cyl_chart::probe::ROWS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    let ends = crate::arrangement::cyl_chart::probe::cell_ends::ROWS
        .lock()
        .expect("the probe's lock is never held across a panic")
        .clone();
    let s1 = |f: fn(&crate::arrangement::cyl_chart::probe::Row) -> usize| {
        chart.iter().map(f).sum::<usize>()
    };
    let s = |f: fn(&crate::arrangement::cyl_chart::probe::cell_ends::Row) -> usize| {
        ends.iter().map(f).sum::<usize>()
    };
    eprintln!(
        "ledger chart: charts {} refused {} cells {}",
        chart.len(),
        chart.iter().filter(|r| r.refused).count(),
        s1(|r| r.cells),
    );
    eprintln!(
        "ledger cell_ends: rows {} refused_booleans {} cells {} end_swapped {} end_disk {} end_exact {} \
             end_other {} end_nocircle {} other_present {} src2_disagree {} src0_present {} exist_disagree {} \
             read_refused {} exist_marks_false {} emit {} emit_unknown {} nocircle_present {} \
             arcs_read {} exact_run_arcs {} arcs_no_mark {} arcs_multi_mark {} emitted_faces {}",
        ends.len(),
        ends.iter().filter(|r| r.emitter_refused).count(),
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
        "ledger crossing: ruling_crossings {} off_max {:e} side_disagree {}",
        hits.len(),
        hits.iter().flat_map(|h| h.off).fold(0.0_f64, f64::max),
        hits.iter().filter(|h| h.side_f64 != h.side).count(),
    );
    // The populations the region emitter moves, attributed to fixtures.
    eprintln!(
        "ledger cell_ends stations: station_pairs {} station_name_failures {}",
        s(|r| r.station_pairs),
        s(|r| r.station_name_failures),
    );
    eprintln!(
        "ledger cell_ends: end_other {} of which single_cut {}",
        s(|r| r.end_other),
        s(|r| r.end_other_single_cut),
    );
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
            "ledger tie: mixed abstentions {} by kind {hist:?}",
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
            "ledger deciding: calls {} exhausted {} probes_tried histogram {tried:?} \
                 offered histogram {offered:?}",
            dec.iter().filter(|r| r.3).count(),
            dec.iter().filter(|r| !r.3).count()
        );
        for r in dec.iter().filter(|r| !r.3) {
            eprintln!(
                "ledger deciding exhausted: {} offered {} ties {:?}",
                r.0, r.2, r.4
            );
        }
        eprintln!(
            "ledger swallowed: the nesting retry swallowed non-abstention errors {}",
            *crate::combinatorics::swallowed_probe::COUNT
                .lock()
                .expect("the probe's lock is never held across a panic")
        );
    }
}
