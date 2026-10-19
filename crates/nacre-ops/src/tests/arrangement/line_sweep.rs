//! The line sweep's truth table: stations in, pieces or a refusal out, no geometry between.

use super::*;

/// Sweep `stations` — `(event, position, name)`, already in `(position, event)` order — and read
/// the pieces back as position pairs.
fn sweep(stations: &[(LineEvent, u8, usize)]) -> Result<Vec<([u8; 2], SegKind)>, DeclineKind> {
    let name = |k: usize| NodeId::three_planes(Canon3::three([0, 1, k + 2]));
    let events: Vec<(LineEvent, (u8, usize))> =
        stations.iter().map(|&(e, at, k)| (e, (at, k))).collect();
    let same: Vec<bool> = stations.windows(2).map(|w| w[0].1 == w[1].1).collect();
    sweep_line(&events, &same, |(_, k)| name(k), MAT).map(|pieces| {
        pieces
            .into_iter()
            .map(|([a, b], kind)| ([a.0, b.0], kind))
            .collect()
    })
}

const MAT: SegKind = SegKind::Transversal { mat: 1 };
const fn graze(body_above: bool) -> SegKind {
    SegKind::Graze { body_above }
}
const fn start(body_above: bool) -> LineEvent {
    LineEvent::GrazeStart { body_above }
}
const fn end(toggle: bool) -> LineEvent {
    LineEvent::GrazeEnd { toggle }
}
const T: LineEvent = LineEvent::Toggle;

/// **What a face leaves on a line, by the events alone.** Each row's pieces are written by
/// hand. The rows a pinched face makes are the ones with two events at one station under one name:
/// two runs meeting end to start (one graze, unbroken), a crossing at a run's end, and two
/// crossings of one point that cancel (the face goes on unbroken across it).
#[test]
fn a_line_sweep_reads_stations_not_events() {
    let ok = |rows: &[(LineEvent, u8, usize)], want: &[([u8; 2], SegKind)]| {
        assert_eq!(sweep(rows), Ok(want.to_vec()), "{rows:?}");
    };
    // Two crossings: present between them.
    ok(&[(T, 0, 0), (T, 1, 1)], &[([0, 1], MAT)]);
    // A run the boundary goes on across: a graze, then present to the next crossing.
    ok(
        &[(start(true), 0, 0), (end(true), 1, 1), (T, 2, 2)],
        &[([0, 1], graze(true)), ([1, 2], MAT)],
    );
    // A run the boundary turns back from: a graze and nothing after it.
    ok(
        &[(start(false), 0, 0), (end(false), 1, 1)],
        &[([0, 1], graze(false))],
    );
    // Two runs meeting end to start at one point under one name: one graze.
    ok(
        &[
            (start(true), 0, 0),
            (end(false), 1, 1),
            (start(true), 1, 1),
            (end(false), 2, 2),
        ],
        &[([0, 2], graze(true))],
    );
    // A crossing at a run's end, under the run's end's name: the run closes before the toggle.
    ok(
        &[
            (start(true), 0, 0),
            (end(false), 1, 1),
            (T, 1, 1),
            (T, 2, 2),
        ],
        &[([0, 1], graze(true)), ([1, 2], MAT)],
    );
    // A crossing at a run's start: the toggle before the run opens.
    ok(
        &[
            (T, 0, 0),
            (T, 1, 1),
            (start(false), 1, 1),
            (end(false), 2, 2),
        ],
        &[([0, 1], MAT), ([1, 2], graze(false))],
    );
    // Two crossings of one point cancel, and the pieces either side join.
    ok(
        &[(T, 0, 0), (T, 1, 1), (T, 1, 1), (T, 2, 2)],
        &[([0, 2], MAT)],
    );

    let refused = |rows: &[(LineEvent, u8, usize)], want: DeclineKind| {
        assert_eq!(sweep(rows), Err(want), "{rows:?}");
    };
    // One point under two names.
    refused(
        &[(T, 0, 0), (T, 1, 1), (T, 1, 2), (T, 2, 3)],
        DeclineKind::CoincidentFeatures,
    );
    // A run's end with no run open; a crossing inside a run; a run inside a run.
    refused(&[(end(false), 0, 0)], DeclineKind::RunSplit);
    refused(
        &[(start(true), 0, 0), (T, 1, 1), (end(false), 2, 2)],
        DeclineKind::RunSplit,
    );
    refused(
        &[
            (start(true), 0, 0),
            (start(true), 1, 1),
            (end(false), 2, 2),
            (end(false), 3, 3),
        ],
        DeclineKind::RunSplit,
    );
    // Still present, or a run still open, past the last station.
    refused(&[(T, 0, 0)], DeclineKind::OddParity);
    refused(&[(start(true), 0, 0)], DeclineKind::OddParity);
}
