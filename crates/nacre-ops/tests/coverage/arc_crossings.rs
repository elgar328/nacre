//! **A wall crossing an arc where the arc's two ends do not say so.**
//!
//! A circle meets a line in at most two points, so a straight wall can cross an arc edge of a
//! flat face twice — its two ends on one side of the wall, its middle bulging across — or once
//! on the way into an end lying on the wall itself. Neither is visible from the arc's two ends,
//! and neither needs a large arc: a slanted wall crosses a quarter-circle fillet twice.
//!
//! Every placement runs over both frames (the world and [`pythagorean_frame`]), `B` inside `A`'s
//! span and past both its caps, and all six booleans: each builds, valid, the volumes add up
//! against the operands' own, and `A ∩ B` is the area the wall cuts off — a circular segment
//! or a band, stated in closed form — over the height the two share.
//!
//! The second shape — one crossing into an end on the wall — is the three-quarter disk's cap read
//! against its own radius class `x = 0`, which every boolean of that disk traces; it is locked
//! there end to end, and where it is read (the walk's table in `curved_nesting`, the trace of a
//! half disk's cap against `x + y = 1` in `tests::arrangement::curved` — the booleans on that
//! pair stop at `RulingBoundNotYet` first).

use crate::common::*;
use crate::stated::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, Profile2d, RejectReason, boolean};

/// One family: `A` extruded over `a` (height 2), `B` the prism over `b`, and the area of
/// `a ∩ b`.
struct Family {
    name: &'static str,
    a: fn() -> Profile2d,
    b: &'static [[f64; 2]],
    common_area: f64,
}

/// `B`'s `(lift, height)` on the frame and the height it shares with `A` (which spans `0…2`).
const HEIGHTS: [(f64, f64, f64); 2] = [(0.5, 1.0, 1.0), (-0.5, 3.0, 2.0)];

/// One placement's six booleans, in the order `A ∪ B`, `B ∪ A`, `A − B`, `B − A`, `A ∩ B`,
/// `B ∩ A`: each built result as `(bodies, volume)` — validated here — or the refusal's reason;
/// and `(|A|, |B|)` as the operands state them.
type Six = ([Result<(usize, f64), RejectReason>; 6], (f64, f64));

/// Every placement of `f`, as `(label, shared height, six)`.
fn every_six(f: &Family) -> Vec<(String, f64, Six)> {
    let mut out = Vec::new();
    for tilted in [false, true] {
        for (lift, height, shared) in HEIGHTS {
            let label = format!(
                "{}, tilted = {tilted}, lift = {lift}, height = {height}",
                f.name
            );
            let mut six = Vec::new();
            let mut operands = (0.0, 0.0);
            for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                for swapped in [false, true] {
                    let (mut m, a, b) = a_solid_and_a_placed_prism(
                        |m| {
                            if tilted {
                                pythagorean_frame(m, Point3::from_array([0.0; 3]))
                            } else {
                                nacre_ops::SketchFrame::world(m, nacre_exact::Axis::Z)
                            }
                        },
                        (f.a)(),
                        f.b,
                        lift,
                        height,
                    );
                    operands = (volume(&m, a), volume(&m, b));
                    let (x, y) = if swapped { (b, a) } else { (a, b) };
                    six.push(match boolean(&mut m, kind, x, y) {
                        Ok(solids) => {
                            m.rebuild_adjacency();
                            let vs = nacre_validate::validate(&m);
                            assert!(vs.is_empty(), "{label}: {kind:?} {swapped}: {vs:?}");
                            Ok((solids.len(), solids.iter().map(|&s| volume(&m, s)).sum()))
                        }
                        Err(BoolError::Rejected { reason, .. }) => Err(reason),
                        Err(e) => panic!("{label}: {kind:?} {swapped}: {e:?}"),
                    });
                }
            }
            let six = six.try_into().expect("six booleans");
            out.push((label, shared, (six, operands)));
        }
    }
    out
}

/// A slot `[−2, 2] × [−1, 1]` whose right end is the half circle about `(2, 0)`.
fn slot() -> Profile2d {
    stated(vec![
        line(p2(-2.0, -1.0), p2(2.0, -1.0)),
        arc_turns(p2(2.0, 0.0), p2(2.0, -1.0), 2),
        line(p2(2.0, 1.0), p2(-2.0, 1.0)),
        line(p2(-2.0, 1.0), p2(-2.0, -1.0)),
    ])
    .expect("a slot")
    .remove(0)
}

/// The plate `[−3, 3]²` with the half disk of [`upper_half_disk`] cut through it.
fn plate_with_a_half_disk_hole() -> Profile2d {
    let mut p = stated(vec![
        line(p2(-3.0, -3.0), p2(3.0, -3.0)),
        line(p2(3.0, -3.0), p2(3.0, 3.0)),
        line(p2(3.0, 3.0), p2(-3.0, 3.0)),
        line(p2(-3.0, 3.0), p2(-3.0, -3.0)),
        arc_turns(p2(0.0, 0.0), p2(1.0, 0.0), 2),
        line(p2(-1.0, 0.0), p2(1.0, 0.0)),
    ])
    .expect("a plate with a hole");
    assert_eq!(p.len(), 1, "one profile, holed");
    p.remove(0)
}

/// The area of the unit circle's segment beyond a chord at distance `d` from its centre.
fn segment(d: f64) -> f64 {
    d.acos() - d * (1.0 - d * d).sqrt()
}

fn families() -> Vec<Family> {
    let pi = std::f64::consts::PI;
    let s = segment(0.5);
    // The fillet's circle (r 2) beyond `x − y = 6`: the centre `(1.5, −2)` is `2.5/√2` from it.
    let fd: f64 = 2.5 / 2.0f64.sqrt() / 2.0;
    vec![
        Family {
            name: "a half disk, the wall y = 0.5 across its arc twice",
            a: upper_half_disk,
            b: &[[-2.0, 0.5], [2.0, 0.5], [2.0, 2.0], [-2.0, 2.0]],
            common_area: s,
        },
        Family {
            name: "a half disk, the wall y = 0.5 keeping the band below it",
            a: upper_half_disk,
            b: &[[-2.0, -2.0], [2.0, -2.0], [2.0, 0.5], [-2.0, 0.5]],
            common_area: pi / 2.0 - s,
        },
        Family {
            name: "a slot, the wall x = 2.5 across its round end twice",
            a: slot,
            b: &[[2.5, -2.0], [4.0, -2.0], [4.0, 2.0], [2.5, 2.0]],
            common_area: s,
        },
        Family {
            name: "a plate, the wall y = 0.5 across its walls and its hole's arc twice",
            a: plate_with_a_half_disk_hole,
            b: &[[-4.0, 0.5], [4.0, 0.5], [4.0, 4.0], [-4.0, 4.0]],
            common_area: 6.0 * 2.5 - s,
        },
        Family {
            name: "a fillet, chamfered by the wall x − y = 6 across it twice",
            a: plate_with_a_fillet,
            b: &[[0.0, -6.0], [6.0, -6.0], [6.0, 0.0]],
            common_area: 4.0 * segment(fd),
        },
    ]
}

/// One placement's six booleans all built — `bodies(i)` bodies for boolean `i` — and their
/// volumes adding up against the operands', with `A ∩ B` the area the wall cuts off over the
/// height the two share.
fn built(label: &str, shared: f64, six: &Six, common_area: f64, bodies: impl Fn(usize) -> usize) {
    let (six, (va, vb)) = six;
    let v = |i: usize| match &six[i] {
        Ok((n, v)) => {
            assert_eq!(*n, bodies(i), "{label}: [{i}] bodies");
            *v
        }
        Err(reason) => panic!("{label}: [{i}] {reason:?}"),
    };
    let near = |x: f64, y: f64, what: &str| {
        assert!((x - y).abs() < 1e-9, "{label}: {what}: {x} vs {y}");
    };
    near(v(0), v(1), "the Fuse commutes");
    near(v(4), v(5), "Common commutes");
    near(v(0) + v(4), va + vb, "Fuse + Common = |A| + |B|");
    near(v(2) + v(4), *va, "(A − B) + Common = |A|");
    near(v(3) + v(4), *vb, "(B − A) + Common = |B|");
    near(v(4), common_area * shared, "Common");
}

/// ★ **Every family builds, all six booleans, in both frames and at both heights** — one body
/// each, except `B − A` where the plate's cut leaves the half disk's segment standing apart from
/// the frame around it (`B` inside `A`'s span); valid, the volumes adding up against the operands',
/// and `A ∩ B` the area the wall cuts off over the height the two share.
#[test]
fn a_wall_crossing_an_arc_twice_builds() {
    let mut ran = 0;
    for f in families() {
        for (label, shared, six) in every_six(&f) {
            let plate = f.name.starts_with("a plate");
            built(&label, shared, &six, f.common_area, |i| {
                if plate && i == 3 && shared == 1.0 {
                    2
                } else {
                    1
                }
            });
            ran += 6;
        }
    }
    assert_eq!(ran, 5 * 2 * 2 * 6, "the family");
}

/// The unit disk less its fourth quadrant — the arc from `(1, 0)` three quarter turns to
/// `(0, −1)`, and the two radii.
fn three_quarter_disk() -> Profile2d {
    stated(vec![
        line(p2(0.0, 0.0), p2(1.0, 0.0)),
        arc_turns(p2(0.0, 0.0), p2(1.0, 0.0), 3),
        line(p2(0.0, -1.0), p2(0.0, 0.0)),
    ])
    .expect("a three-quarter disk")
    .remove(0)
}

/// ★ **A three-quarter disk, whose arc winds its cap.** Its caps are three points — the centre and
/// the arc's two ends — whose chord triangle turns the other way round from the face, and its own
/// radius class `x = 0` meets the cap's arc once inside, at `(0, 1)`, on the way into the end
/// `(0, −1)` that lies on the line. Two boxes: the wall `x = 0.5` across the arc once (`A ∩ B` is
/// the upper half of the segment beyond `x = 0.5`) and the wall `y = 0.5` across it twice (the
/// whole segment beyond `y = 0.5`). Every boolean builds one body, the volumes adding up — except
/// on the tilted frame with `B` inside `A`'s span, where all six refuse by the chart's name
/// (`CylinderGateUndecided`: a lateral cell whose two ends read different chambers, recorded in
/// `todo.md`); that name is pinned here, cell by cell.
#[test]
fn a_three_quarter_disk_builds_but_on_the_tilted_frame_inside_its_span() {
    let s = segment(0.5);
    let sector = [
        Family {
            name: "a three-quarter disk, the wall x = 0.5 across its arc once",
            a: three_quarter_disk,
            b: &[[0.5, -2.0], [2.0, -2.0], [2.0, 2.0], [0.5, 2.0]],
            common_area: s / 2.0,
        },
        Family {
            name: "a three-quarter disk, the wall y = 0.5 across its arc twice",
            a: three_quarter_disk,
            b: &[[-2.0, 0.5], [2.0, 0.5], [2.0, 2.0], [-2.0, 2.0]],
            common_area: s,
        },
    ];
    let mut refused = 0;
    for f in &sector {
        // `every_six`'s order: the world frame then the tilted one, each `B` inside the span
        // then past both caps.
        for (k, (label, shared, six)) in every_six(f).into_iter().enumerate() {
            let tilted_inside = k == 2;
            if tilted_inside {
                for (i, r) in six.0.iter().enumerate() {
                    assert!(
                        matches!(r, Err(RejectReason::CylinderGateUndecided)),
                        "{label}: [{i}] the refusal moved — {r:?}"
                    );
                    refused += 1;
                }
                continue;
            }
            built(&label, shared, &six, f.common_area, |_| 1);
        }
    }
    assert_eq!(refused, 2 * 6, "the tilted cells inside the span");
}
