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
//! of `B` along `A`'s axis ([`Height`]).

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, Profile2d, RejectClass, RejectReason, boolean};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// Where `B` stands along `A`'s axis — `A` spans `0…2` on the frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Height {
    /// On `A`'s base cap.
    OnBase,
    /// Inside `A`'s span.
    Inside,
    /// Past `A`'s top.
    PastTop,
    /// Past both of `A`'s caps.
    PastBoth,
}

impl Height {
    const ALL: [Height; 4] = [
        Height::OnBase,
        Height::Inside,
        Height::PastTop,
        Height::PastBoth,
    ];

    /// `B`'s `(lift, height)` on the frame.
    fn span(self) -> (f64, f64) {
        match self {
            Height::OnBase => (0.0, 1.0),
            Height::Inside => (0.5, 1.0),
            Height::PastTop => (0.5, 2.0),
            Height::PastBoth => (-0.5, 3.0),
        }
    }

    /// The height `A` and `B` share.
    fn shared(self) -> f64 {
        let (lift, height) = self.span();
        (lift + height).min(2.0) - lift.max(0.0)
    }
}

/// One placement of `B` against `A` — what the six booleans of [`every_six`] share. The facts
/// travel as values; `Display` is the label, for messages only.
#[derive(Clone, Copy, Debug)]
struct Placement<'q> {
    family: &'q str,
    /// `B`'s profile as written (its area does not depend on the vertex order).
    quad: &'q [[f64; 2]],
    tilted: bool,
    rotated: bool,
    height: Height,
}

impl Placement<'_> {
    /// `|B|`.
    fn b_volume(&self) -> f64 {
        let q = self.quad;
        let area = (0..q.len())
            .map(|k| {
                let (p, r) = (q[k], q[(k + 1) % q.len()]);
                p[0] * r[1] - r[0] * p[1]
            })
            .sum::<f64>()
            .abs()
            / 2.0;
        area * self.height.span().1
    }
}

impl std::fmt::Display for Placement<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (lift, height) = self.height.span();
        write!(
            f,
            "{}, tilted = {}, rotated = {}, {:?} (lift = {lift}, height = {height})",
            self.family, self.tilted, self.rotated, self.height
        )
    }
}

/// One boolean of a [`Placement`]: the operation and whether `B` comes first.
struct Case<'a, 'q> {
    placement: &'a Placement<'q>,
    kind: BoolKind,
    swapped: bool,
}

impl std::fmt::Display for Case<'_, '_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}, {:?} swapped = {}",
            self.placement, self.kind, self.swapped
        )
    }
}

/// Every boolean of `A` (extruded over `a`) against the prism over each of `quads`, over both
/// frames, both vertex orders and every [`Height`], handed to `check` as a [`Case`]. Returns how
/// many booleans ran.
fn every_placement<'q>(
    a: fn() -> Profile2d,
    quads: &[(&'q str, &'q [[f64; 2]])],
    mut check: impl FnMut(&Case<'_, 'q>, &mut Model, Result<Vec<Handle<Solid>>, BoolError>),
) -> usize {
    let mut ran = 0;
    for tilted in [false, true] {
        for &(family, quad) in quads {
            for rotated in [false, true] {
                let mut written = quad.to_vec();
                if rotated {
                    written.rotate_left(1);
                }
                for height in Height::ALL {
                    let placement = Placement {
                        family,
                        quad,
                        tilted,
                        rotated,
                        height,
                    };
                    let (lift, tall) = height.span();
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
                                &written,
                                lift,
                                tall,
                            );
                            let (x, y) = if swapped { (y, x) } else { (x, y) };
                            let got = boolean(&mut m, kind, x, y);
                            let case = Case {
                                placement: &placement,
                                kind,
                                swapped,
                            };
                            check(&case, &mut m, got);
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
fn refused_today(c: &Case<'_, '_>, _: &mut Model, got: Result<Vec<Handle<Solid>>, BoolError>) {
    match got {
        Err(BoolError::Rejected { reason, .. }) => assert_eq!(
            reason.class(),
            RejectClass::NotSupported,
            "{c}: refused as {reason:?}"
        ),
        other => panic!("{c}: the wall moved — {other:?}"),
    }
}

/// One placement's six booleans, in the order `A ∪ B`, `B ∪ A`, `A − B`, `B − A`, `A ∩ B`,
/// `B ∩ A`: each built result as `(bodies, volume)` — validated here — or the refusal's reason.
type Six = [Result<(usize, f64), RejectReason>; 6];

/// Every placement of the prism over each of `quads` against `A` over `a`, as in
/// [`every_placement`], gathered into [`Six`] per [`Placement`].
fn every_six<'q>(
    a: fn() -> Profile2d,
    quads: &[(&'q str, &'q [[f64; 2]])],
) -> Vec<(Placement<'q>, Six)> {
    let mut out: Vec<(Placement<'q>, Six)> = Vec::new();
    let mut cur: Vec<Result<(usize, f64), RejectReason>> = Vec::new();
    every_placement(a, quads, |c, m, got| {
        cur.push(match got {
            Ok(solids) => {
                m.rebuild_adjacency();
                let vs = nacre_validate::validate(m);
                assert!(vs.is_empty(), "{c}: {vs:?}");
                Ok((solids.len(), solids.iter().map(|&s| volume(m, s)).sum()))
            }
            Err(BoolError::Rejected { reason, .. }) => Err(reason),
            Err(e) => panic!("{c}: {e:?}"),
        });
        if cur.len() == 6 {
            let six: Six = std::mem::take(&mut cur).try_into().expect("six booleans");
            out.push((*c.placement, six));
        }
    });
    out
}

/// The volume identities a wrong answer `validate` cannot see would break, wherever the booleans
/// they need built: `Fuse + Common = |A| + |B|`, `(A − B) + Common = |A|`, `(B − A) + Common =
/// |B|`, and each pair of swapped operands alike.
fn identities(p: &Placement<'_>, six: &Six) {
    let (a, b) = (2.0 * std::f64::consts::PI, p.b_volume());
    let v = |i: usize| six[i].as_ref().ok().map(|&(_, v)| v);
    for (x, y) in [(0, 1), (4, 5)] {
        if let (Some(r), Some(q)) = (v(x), v(y)) {
            assert!((r - q).abs() < 1e-9, "{p}: {x} vs {y}: {r} vs {q}");
        }
    }
    for (i, want) in [(0, a + b), (2, a), (3, b)] {
        if let (Some(r), Some(c)) = (v(i), v(4)) {
            assert!(
                (r + c - want).abs() < 1e-9,
                "{p}: [{i}] + Common = {}, want {want}",
                r + c
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
    for (p, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match r {
                Ok((1, _)) => {}
                other => panic!("{p}: [{i}] {other:?}"),
            }
        }
        identities(p, six);
        if p.family == "on the seam" {
            let common = six[4].as_ref().expect("built").1;
            let want = segment * p.height.shared();
            assert!(
                (common - want).abs() < 1e-9,
                "{p}: Common {common}, want {want}"
            );
        }
    }
    let off: &[[f64; 2]] = &[[-1.0, 1.6], [1.0, 0.001], [3.0, 0.001], [3.0, 3.6]];
    every_placement(unit_disk, &[("0.001 off", off)], |c, _, got| {
        got.unwrap_or_else(|e| panic!("{c}: {e:?}"));
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
    for (p, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match r {
                Ok((1, _)) => {}
                other => panic!("{p}: [{i}] {other:?}"),
            }
        }
        identities(p, six);
        let common = six[4].as_ref().expect("built").1;
        let want = half_disk * p.height.shared();
        assert!(
            (common - want).abs() < 1e-9,
            "{p}: Common {common}, want {want}"
        );
    }
}

/// ★ **A tangent wall on the edge of a flat through `A`'s own axis.** `A` is the half cylinder
/// over [`upper_half_disk`], so its flat `y = 0` meets its lateral in the ruling `(1, 0)` — `A`'s
/// own edge, where the lateral ends. `B`'s wall `x = 1` is tangent to `A`'s cylinder along that
/// edge, either running across it or ending on it (then `B`'s face `y = 0` lies on `A`'s flat).
/// Refused today, by name — and the name is the placement's, never the operation's or the vertex
/// order's: every boolean stops before the operation decides anything.
///
/// - **Across**, `A` first: `RulingBoundNotYet`, at every height and in both frames.
/// - **Across**, `B` first: the arrangement's backstops — `RingOrientation` past both caps,
///   `LabelConflict` on the base cap. Inside the span and past the top, **the frame decides the
///   name**: `LabelConflict` in the world frame, the cap's winding doubling back
///   (`StraightAngle`) in the tilted one.
/// - **Ending**, either order and frame: `StraightAngle` inside the span and past the top,
///   `LabelConflict` on the base cap and past both caps.
///
/// The backstop cells ([`backstop`], `SuspectedDefect`, 84 booleans) are pinned by variant, so a
/// backstop that moves to another placement is news; the rest by class.
#[test]
fn a_tangent_wall_on_a_flats_edge_is_refused_today() {
    let across: &[[f64; 2]] = &[[-3.0, -2.0], [1.0, -2.0], [1.0, 0.5], [-3.0, 0.5]];
    let ending: &[[f64; 2]] = &[[-3.0, 0.0], [1.0, 0.0], [1.0, 0.5], [-3.0, 0.5]];
    let mut defects = 0;
    let ran = every_placement(
        upper_half_disk,
        &[("across the edge", across), ("ending on the edge", ending)],
        |c, m, got| match backstop(c) {
            Some(want) => {
                match &got {
                    Err(BoolError::Rejected { reason, .. }) => {
                        assert_eq!(*reason, want, "{c}: the backstop moved");
                    }
                    other => panic!("{c}: the backstop moved — {other:?}"),
                }
                defects += 1;
            }
            None => refused_today(c, m, got),
        },
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
    assert_eq!(defects, 84, "the backstop refusals");
}

/// The half-cylinder family's backstop cells — where [`a_tangent_wall_on_a_flats_edge_is_refused_today`]
/// stops at `SuspectedDefect` today, and with which reason. Neither the operation nor the vertex
/// order enters.
fn backstop(c: &Case<'_, '_>) -> Option<RejectReason> {
    use Height::*;
    let p = c.placement;
    match (p.family, c.swapped, p.height, p.tilted) {
        ("ending on the edge", _, OnBase | PastBoth, _) => Some(RejectReason::LabelConflict),
        ("across the edge", true, PastBoth, _) => Some(RejectReason::RingOrientation),
        ("across the edge", true, OnBase, _)
        | ("across the edge", true, Inside | PastTop, false) => Some(RejectReason::LabelConflict),
        _ => None,
    }
}

/// ★ **Both walls entering `A` — an inward wedge**, its apex on the lateral at the seam and off
/// it. Near the line, `A − B` is two lobes meeting on it alone: where they join elsewhere (`B`
/// inside `A`'s span or reaching past one cap) that is one solid touching itself and is refused
/// (`NonManifoldResultEdge`, or `ArcBoundNotYet` on the seam with `B` past the top); where `B`
/// cuts `A` through (past both caps) it is two solids. **It is never one body** — the lateral run
/// through the line as one face would hide the touch under a closed shell. Common, `B − A` and
/// the Fuse build one body; the volumes add up. With `B`'s base on `A`'s cap (lift 0) and the apex
/// on the seam, `A − B` is refused `ZeroLengthEdge`: `SuspectedDefect`, 4 booleans (the frames),
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
    for (p, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match (i, r) {
                (2, Ok((n, _))) => {
                    assert_eq!(*n, 2, "{p}: A − B is two bodies or none");
                    assert_eq!(p.height, Height::PastBoth, "{p}: A − B built");
                }
                (_, Ok((n, _))) => assert_eq!(*n, 1, "{p}: [{i}]"),
                (_, Err(reason)) if reason.class() == RejectClass::SuspectedDefect => {
                    assert!(
                        i == 2
                            && p.family == "on the seam"
                            && p.height == Height::OnBase
                            && *reason == RejectReason::ZeroLengthEdge,
                        "{p}: [{i}] {reason:?}"
                    );
                    defects += 1;
                }
                (2, Err(reason)) => assert!(
                    match reason {
                        RejectReason::NonManifoldResultEdge => true,
                        RejectReason::ArcBoundNotYet => {
                            p.family == "on the seam" && p.height == Height::PastTop
                        }
                        _ => false,
                    },
                    "{p}: A − B refused as {reason:?}"
                ),
                (_, Err(reason)) => panic!("{p}: [{i}] {reason:?}"),
            }
        }
        identities(p, six);
    }
    assert_eq!(defects, 4, "the rim placement's refusals");
}

/// ★ **The corner touching `A` from outside — a line contact.** Common would be empty, Cut the
/// operand, Fuse two bodies touching along the line. Refused today.
#[test]
fn a_line_contact_on_a_lateral_is_refused_today() {
    let touching: &[[f64; 2]] = &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]];
    let ran = every_placement(unit_disk, &[("a line contact", touching)], refused_today);
    assert_eq!(ran, 2 * 2 * 4 * 6, "the family");
}

/// ★ **A tangent face that ends on the line at one height and crosses it at another.** `B` is
/// an L-shaped profile in the plane `x = −1` extruded to `x = 1`: its caps `x = ±1` are tangent
/// to `A`'s unit cylinder, each an L whose concave corner puts one stretch of its boundary on the
/// line `(±1, 0)` — beside `B`'s wall `y = 0`, a secant through `A`'s axis — while below that
/// stretch the line runs through the cap's interior. That is the face the tangency verdict's
/// `line_is_an_edge` keeps out (`!straddles`: an edge on the line does not make the contact an
/// edge where the face runs across it). Today the rulings split refuses it first — the cap's
/// segment crosses the line (`RulingBoundNotYet`), both vertex orders, all six booleans — so the
/// verdict's guard is not reached; when this moves, the guard speaks.
#[test]
fn a_tangent_face_across_its_edge_on_the_line_is_refused_today() {
    use nacre_math::Vector3;
    use nacre_ops::{OpOutput, Operation, SketchPlane, apply};
    let l: &[[f64; 2]] = &[
        [-1.0, 0.5],
        [3.0, 0.5],
        [3.0, 1.5],
        [0.0, 1.5],
        [0.0, 1.0],
        [-1.0, 1.0],
    ];
    let mut ran = 0;
    for rotated in [false, true] {
        let mut l = l.to_vec();
        if rotated {
            l.rotate_left(1);
        }
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            for swapped in [false, true] {
                let mut m = Model::new();
                let extrude = |m: &mut Model, frame, profile, dist| {
                    let Ok(OpOutput::Extrude { solid, .. }) = apply(
                        m,
                        &Operation::Extrude {
                            frame,
                            profile,
                            dist,
                        },
                    ) else {
                        panic!("extrude")
                    };
                    m.rebuild_adjacency();
                    solid
                };
                let world = nacre_ops::SketchFrame::world(&m, nacre_exact::Axis::Z);
                let a = extrude(&mut m, world, unit_disk(), 2.0);
                let side = datum_frame(
                    &mut m,
                    SketchPlane::from_axes(
                        Point3::from_array([-1.0, 0.0, 0.0]),
                        Vector3::from_array([0.0, 1.0, 0.0]),
                        Vector3::from_array([0.0, 0.0, 1.0]),
                    ),
                );
                let profile =
                    Profile2d::polygon(l.iter().map(|q| p2(q[0], q[1])).collect()).expect("an L");
                let b = extrude(&mut m, side, profile, 2.0);
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let at = format!("rotated = {rotated}, {kind:?} swapped = {swapped}");
                match boolean(&mut m, kind, x, y) {
                    Err(BoolError::Rejected { reason, .. }) => assert_eq!(
                        reason,
                        RejectReason::RulingBoundNotYet,
                        "{at}: the wall moved"
                    ),
                    other => panic!("{at}: the wall moved — {other:?}"),
                }
                ran += 1;
            }
        }
    }
    assert_eq!(ran, 2 * 6, "the family");
}
