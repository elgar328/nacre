//! **One solid's edge on a ruling of the other's cylinder.**
//!
//! Two planes along a cylinder's axis meet in a line along it, so where a prism's corner lies on
//! the cylinder, the prism's edge there lies on a ruling — the whole line, not a point. Two
//! families put it there:
//!
//! - **two secant walls**: both of `B`'s walls at the corner cut `A`'s lateral, one entering `A`;
//! - **a tangent wall and a secant one — the keyhole**: a box as wide as the cylinder's diameter,
//!   its side walls tangent, its corners on the circle. The commonest shape of the two.
//!
//! Both families build.
//!
//! Plus the shapes whose honest answer is partly a refusal by name: both walls entering `A` (an
//! inward wedge — its `A − B` touches itself along the line), and the corner touching `A` from
//! outside (a line contact).
//!
//! Every family runs over both frames (the world and [`pythagorean_frame`]), the polygon as
//! written and rotated by one vertex (which reorders the plane classes), and four placements
//! of `B` along `A`'s axis ([`HEIGHTS`]).

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, Profile2d, RejectClass, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// Where `B` stands along `A`'s axis, `(lift, height)` — `A` spans `0…2` on the frame. On `A`'s
/// base cap; inside `A`'s span; past `A`'s top; past both caps.
const HEIGHTS: [(f64, f64); 4] = [(0.0, 1.0), (0.5, 1.0), (0.5, 2.0), (-0.5, 3.0)];

/// Every boolean of `A` (extruded over `a`) against the prism over each of `quads`, over both
/// frames, both vertex orders and every [`HEIGHTS`] placement, handed to `check` with a label.
/// Returns how many booleans ran.
fn every_placement(
    a: fn() -> Profile2d,
    quads: &[(&str, &[[f64; 2]])],
    mut check: impl FnMut(&str, &mut Model, Result<Vec<Handle<Solid>>, BoolError>),
) -> usize {
    let mut ran = 0;
    for tilted in [false, true] {
        for &(what, quad) in quads {
            for rotated in [false, true] {
                let mut quad = quad.to_vec();
                if rotated {
                    quad.rotate_left(1);
                }
                for (lift, height) in HEIGHTS {
                    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                        for swapped in [false, true] {
                            let (mut m, x, y) = a_solid_and_a_placed_prism(
                                |m| {
                                    if tilted {
                                        pythagorean_frame(m, Point3::from_array([0.0; 3]))
                                    } else {
                                        nacre_ops::SketchFrame::world(m, nacre_exact::Axis::Z)
                                    }
                                },
                                a(),
                                &quad,
                                lift,
                                height,
                            );
                            let (x, y) = if swapped { (y, x) } else { (x, y) };
                            let at = format!(
                                "{what}, tilted = {tilted}, rotated = {rotated}, \
                                 lift = {lift}, height = {height}, {kind:?} swapped = {swapped}"
                            );
                            let got = boolean(&mut m, kind, x, y);
                            check(&at, &mut m, got);
                            ran += 1;
                        }
                    }
                }
            }
        }
    }
    ran
}

/// Today's answer: refused, as a coverage limit. The class is asserted, not the variant —
/// reason names are engine vocabulary; a refusal that turns `Impossible` or `SuspectedDefect`, or
/// a placement that starts to build, is news here and in the todo item this pins.
fn refused_today(at: &str, _: &mut Model, got: Result<Vec<Handle<Solid>>, BoolError>) {
    match got {
        Err(BoolError::Rejected { reason, .. }) => assert_eq!(
            reason.class(),
            RejectClass::NotSupported,
            "{at}: refused as {reason:?}"
        ),
        other => panic!("{at}: the wall moved — {other:?}"),
    }
}

/// One placement's six booleans, in the order `A ∪ B`, `B ∪ A`, `A − B`, `B − A`, `A ∩ B`,
/// `B ∩ A`: each built result as `(bodies, volume)` — validated here — or the refusal's reason.
type Six = [Result<(usize, f64), nacre_ops::RejectReason>; 6];

/// Every placement of the prism over each of `quads` against `A` over `a`, as in
/// [`every_placement`], gathered into [`Six`] per placement with its label and `B`'s volume.
fn every_six(a: fn() -> Profile2d, quads: &[(&str, &[[f64; 2]])]) -> Vec<(String, f64, f64, Six)> {
    let mut out: Vec<(String, f64, f64, Six)> = Vec::new();
    let mut cur: Vec<Result<(usize, f64), nacre_ops::RejectReason>> = Vec::new();
    every_placement(a, quads, |at, m, got| {
        cur.push(match got {
            Ok(solids) => {
                m.rebuild_adjacency();
                let vs = nacre_validate::validate(m);
                assert!(vs.is_empty(), "{at}: {vs:?}");
                Ok((solids.len(), solids.iter().map(|&s| volume(m, s)).sum()))
            }
            Err(BoolError::Rejected { reason, .. }) => Err(reason),
            Err(e) => panic!("{at}: {e:?}"),
        });
        if cur.len() == 6 {
            let label = at.rsplit_once(", ").map_or(at, |(l, _)| l).to_string();
            let quad = quads
                .iter()
                .find(|(w, _)| at.starts_with(w))
                .expect("the label starts with the family")
                .1;
            let area = (0..quad.len())
                .map(|k| {
                    let (p, q) = (quad[k], quad[(k + 1) % quad.len()]);
                    p[0] * q[1] - q[0] * p[1]
                })
                .sum::<f64>()
                .abs()
                / 2.0;
            let height: f64 = label
                .split("height = ")
                .nth(1)
                .and_then(|h| h.parse().ok())
                .expect("the label states the height");
            let lift: f64 = label
                .split("lift = ")
                .nth(1)
                .and_then(|l| l.split(',').next())
                .and_then(|l| l.parse().ok())
                .expect("the label states the lift");
            let six: Six = std::mem::take(&mut cur).try_into().expect("six booleans");
            out.push((
                label,
                area * height,
                (lift + height).min(2.0) - lift.max(0.0),
                six,
            ));
        }
    });
    out
}

/// The volume identities a wrong answer `validate` cannot see would break, wherever the booleans
/// they need built: `Fuse + Common = |A| + |B|`, `(A − B) + Common = |A|`, `(B − A) + Common =
/// |B|`, and each pair of swapped operands alike.
fn identities(label: &str, b: f64, six: &Six) {
    let a = 2.0 * std::f64::consts::PI;
    let v = |i: usize| six[i].as_ref().ok().map(|&(_, v)| v);
    for (x, y) in [(0, 1), (4, 5)] {
        if let (Some(p), Some(q)) = (v(x), v(y)) {
            assert!((p - q).abs() < 1e-9, "{label}: {x} vs {y}: {p} vs {q}");
        }
    }
    for (i, want) in [(0, a + b), (2, a), (3, b)] {
        if let (Some(p), Some(c)) = (v(i), v(4)) {
            assert!(
                (p + c - want).abs() < 1e-9,
                "{label}: [{i}] + Common = {}, want {want}",
                p + c
            );
        }
    }
}

/// ★ **Two secant walls, one entering `A` — the corner on the rim's seam (`(1, 0)`, `θ = 0` on
/// both frames) and off it (`(0.8, −0.6)`).** Every boolean builds one body, valid, and the
/// volumes add up; on the seam, `A ∩ B` is the circular segment beyond the slanted wall
/// `0.8x + y = 0.8` (its distance from the axis `d = 0.8/√1.64`, area `acos d − d·√(1 − d²)`)
/// over the height the two solids share. A corner `0.001` off the lateral builds as well.
#[test]
fn a_corner_entering_a_lateral_builds() {
    let on_the_seam: &[[f64; 2]] = &[[-1.0, 1.6], [1.0, 0.0], [3.0, 0.0], [3.0, 3.6]];
    let off_the_seam: &[[f64; 2]] = &[[-1.0, 1.6], [0.8, -0.6], [3.0, -0.6], [3.0, 3.6]];
    let all = every_six(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
    );
    assert_eq!(all.len(), 2 * 2 * 2 * 4, "the family");
    let d = 0.8 / 1.64f64.sqrt();
    let segment = d.acos() - d * (1.0 - d * d).sqrt();
    for (label, b, shared, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match r {
                Ok((1, _)) => {}
                other => panic!("{label}: [{i}] {other:?}"),
            }
        }
        identities(label, *b, six);
        if label.starts_with("on the seam") {
            let common = six[4].as_ref().expect("built").1;
            assert!(
                (common - segment * shared).abs() < 1e-9,
                "{label}: Common {common}, want {}",
                segment * shared
            );
        }
    }
    let off: &[[f64; 2]] = &[[-1.0, 1.6], [1.0, 0.001], [3.0, 0.001], [3.0, 3.6]];
    every_placement(unit_disk, &[("0.001 off", off)], |at, _, got| {
        got.unwrap_or_else(|e| panic!("{at}: {e:?}"));
    });
}

/// ★ **The keyhole**: the box `[−1, 1] × [0, 3]` against the unit cylinder — its walls `x = ±1`
/// tangent, its wall `y = 0` through the axis, its corners `(±1, 0)` on the circle. Every boolean
/// builds one body, valid, and the volumes add up; `A ∩ B` is the half disk over the height the
/// two share. At the caps the rim's arc and the tangent wall leave each corner tangent and the
/// same way (ordered by the arc's bending), `B − A`'s cap has a cusp there, and the tangency's
/// verdict leaves the line — an edge of both solids' faces — to the assembly's structure.
#[test]
fn a_keyhole_builds() {
    let keyhole: &[[f64; 2]] = &[[-1.0, 0.0], [1.0, 0.0], [1.0, 3.0], [-1.0, 3.0]];
    let all = every_six(unit_disk, &[("the keyhole", keyhole)]);
    assert_eq!(all.len(), 2 * 2 * 4, "the family");
    let half_disk = std::f64::consts::PI / 2.0;
    for (label, b, shared, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match r {
                Ok((1, _)) => {}
                other => panic!("{label}: [{i}] {other:?}"),
            }
        }
        identities(label, *b, six);
        let common = six[4].as_ref().expect("built").1;
        assert!(
            (common - half_disk * shared).abs() < 1e-9,
            "{label}: Common {common}, want {}",
            half_disk * shared
        );
    }
}

/// ★ **A tangent wall on the edge of a flat through `A`'s own axis.** `A` is the half cylinder
/// over [`upper_half_disk`], so its flat `y = 0` meets its lateral in the ruling `(1, 0)` — `A`'s
/// own edge, where the lateral ends. `B`'s wall `x = 1` is tangent to `A`'s cylinder along that
/// edge, either running across it or ending on it (then `B`'s face `y = 0` lies on `A`'s flat).
/// Refused today, by name: a wall running across the line reaches the rulings split as a
/// stretch it does not arrange (`RulingBoundNotYet`) or the cap's winding as a doubling back
/// (`StraightAngle`); and 84 booleans stop at the arrangement's backstops — `LabelConflict` and,
/// past both caps, `RingOrientation` (`SuspectedDefect`), pinned here by count.
#[test]
fn a_tangent_wall_on_a_flats_edge_is_refused_today() {
    let across: &[[f64; 2]] = &[[-3.0, -2.0], [1.0, -2.0], [1.0, 0.5], [-3.0, 0.5]];
    let ending: &[[f64; 2]] = &[[-3.0, 0.0], [1.0, 0.0], [1.0, 0.5], [-3.0, 0.5]];
    let mut defects = 0;
    let ran = every_placement(
        upper_half_disk,
        &[("across the edge", across), ("ending on the edge", ending)],
        |at, m, got| match &got {
            Err(BoolError::Rejected { reason, .. })
                if reason.class() == RejectClass::SuspectedDefect =>
            {
                defects += 1;
            }
            _ => refused_today(at, m, got),
        },
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
    assert_eq!(defects, 84, "the backstop refusals");
}

/// ★ **Both walls entering `A` — an inward wedge**, its apex on the lateral at the seam and off
/// it. Near the line, `A − B` is two lobes meeting on it alone: where they join elsewhere (`B`
/// inside `A`'s span or reaching past one cap) that is one solid touching itself and is refused
/// (`NonManifoldResultEdge`, or `ArcBoundNotYet` on the seam with `B` past the top); where `B`
/// cuts `A` through (past both caps) it is two solids. **It is never one body** — the lateral run
/// through the line as one face would hide the touch under a closed shell. Common, `B − A` and
/// the lifted Fuse build one body; the volumes add up. With `B`'s base on `A`'s cap (lift 0) the
/// Fuse is refused `OpenResultShell` — the base cap's coplanar merge keeps the inner rim arc —
/// and `A − B` with the apex on the seam `ZeroLengthEdge`: `SuspectedDefect`, 20 booleans,
/// pinned by count.
#[test]
fn an_inward_wedge_on_a_lateral_builds_or_is_refused_by_name() {
    let on_the_seam: &[[f64; 2]] = &[[1.0, 0.0], [-3.0, 3.0], [-3.0, -3.0]];
    let off_the_seam: &[[f64; 2]] = &[[0.8, -0.6], [-3.0, 3.0], [-3.0, -3.0]];
    let all = every_six(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
    );
    assert_eq!(all.len(), 2 * 2 * 2 * 4, "the family");
    let mut defects = 0;
    for (label, b, _, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match (i, r) {
                (2, Ok((n, _))) => {
                    assert_eq!(*n, 2, "{label}: A − B is two bodies or none");
                    assert!(label.contains("lift = -0.5"), "{label}: A − B built");
                }
                (_, Ok((n, _))) => assert_eq!(*n, 1, "{label}: [{i}]"),
                (_, Err(reason)) if reason.class() == RejectClass::SuspectedDefect => {
                    assert!(label.contains("lift = 0,"), "{label}: [{i}] {reason:?}");
                    defects += 1;
                }
                (2, Err(reason)) => assert!(
                    matches!(
                        reason,
                        nacre_ops::RejectReason::NonManifoldResultEdge
                            | nacre_ops::RejectReason::ArcBoundNotYet
                    ),
                    "{label}: A − B refused as {reason:?}"
                ),
                (_, Err(reason)) => panic!("{label}: [{i}] {reason:?}"),
            }
        }
        identities(label, *b, six);
    }
    assert_eq!(defects, 20, "the rim placement's refusals");
}

/// ★ **The corner touching `A` from outside — a line contact.** Common would be empty, Cut the
/// operand, Fuse two bodies touching along the line. Refused today.
#[test]
fn a_line_contact_on_a_lateral_is_refused_today() {
    let touching: &[[f64; 2]] = &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]];
    let ran = every_placement(unit_disk, &[("a line contact", touching)], refused_today);
    assert_eq!(ran, 2 * 2 * 4 * 6, "the family");
}
