//! The region walk's own rules, asked of hand-built charts.

use super::*;

/// Two lines crossing at event `1`: line `0` below line `1` in slab `0`, above it in slab `1`.
fn crossing() -> Chart {
    Chart {
        lines: vec![7, 9],
        slabs: vec![vec![0, 1], vec![1, 0]],
        ..Chart::default()
    }
}

fn arc(line: usize, unit: usize, ccw: bool) -> Edge {
    Edge::Arc { line, unit, ccw }
}

fn ruling(up: bool) -> Edge {
    Edge::Ruling {
        station: 1,
        lo: 0,
        hi: 1,
        up,
    }
}

/// **With one way out in each direction the order is left, straight, right, back** — the grid's
/// rule, which the corner's angular order reduces to. Arriving `+θ` along line `0` over slab `0`:
/// up the ruling, on along line `0`, down the ruling, back along line `0`.
#[test]
fn one_way_each_direction_is_left_straight_right_back() {
    let chart = crossing();
    let came = arc(0, 0, true);
    let key = |e| departure_key(&chart, came, e);
    let (up, on, down, back) = (
        ruling(true),
        arc(0, 1, true),
        ruling(false),
        arc(0, 0, false),
    );
    assert!(key(up) < key(on) && key(on) < key(down) && key(down) < key(back));
}

/// **Two lines leaving one way: the more left goes first.** Leaving event `1` `+θ` (into slab `1`)
/// the higher line, `0`, is the more left; leaving it `−θ` (into slab `0`) the lower one, `0`
/// again — read against an arrival up the ruling, which turns left to `−θ` and right to `+θ`.
#[test]
fn two_lines_leaving_one_way_go_by_the_slab_they_enter() {
    let chart = crossing();
    let key = |came, e| departure_key(&chart, came, e);
    let came = ruling(false); // arrived going down: back is up, left is `+θ`
    assert!(
        key(came, arc(0, 1, true)) < key(came, arc(1, 1, true)),
        "the higher line leaving +θ"
    );
    let came = ruling(true); // arrived going up: left is `−θ`
    assert!(
        key(came, arc(0, 0, false)) < key(came, arc(1, 0, false)),
        "the lower line leaving −θ"
    );
}

/// **A wedge closing at a crossing turns back along the other line.** Line `1` lies above line `0`
/// west of the crossing (slab `0`), and the region between them pinches to it. Arriving there
/// `+θ` along line `0` with the region above, the boundary goes on `−θ` along line `1` — the
/// sharpest left, ahead of the ruling up and of going straight on along either line.
#[test]
fn a_wedge_closing_at_a_crossing_turns_back_along_the_other_line() {
    let chart = crossing();
    let key = |e| departure_key(&chart, arc(0, 0, true), e);
    let back_along_1 = key(arc(1, 0, false));
    assert!(
        back_along_1 < key(ruling(true)),
        "back into the wedge before the ruling up"
    );
    assert!(
        back_along_1 < key(arc(0, 1, true)),
        "before going on along line 0"
    );
    assert!(
        back_along_1 < key(arc(1, 1, true)),
        "before going on along line 1"
    );
}
