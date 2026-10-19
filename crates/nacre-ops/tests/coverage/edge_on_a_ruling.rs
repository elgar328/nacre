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
//! Plus the shapes whose corner stands on `A`'s rim: both walls entering `A` (an inward wedge —
//! its `A − B` touches itself along the line, a refusal by name) and a corner touching `A` from
//! outside (a line contact, a tangent side), which builds — see
//! [`a_vertex_on_a_cap_circle_builds_or_is_refused_by_name`]. And a V-notch whose apex lies on
//! the lateral — a box round `A`, the two walls opening away from the axis: its `B − A` touches
//! itself along the line, a refusal by name — see
//! [`a_notch_whose_apex_is_on_a_lateral_builds_or_is_refused_by_name`].
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
                let v = solids.iter().map(|&s| volume(m, s)).sum();
                drawn(c, m, &solids, v);
                Ok((solids.len(), v))
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

/// **Every built result is drawn** — what `boolean`'s own exit asserts in the library's tests
/// (`tess_census`), which these front-door tests do not pass through. The mesh exists, is
/// watertight, and its signed volume is the solids' own up to what sampling cuts off a curved face:
/// at most the face's area times the sagitta the budget allows, the smaller of `TessConfig::tol`
/// and `r(1 − cos(max_angle / 2))` — the two budgets `TessConfig` states. The slack comes from the
/// budget, not from the fixture: `tol` alone would allow `0.126` on the unit cylinder's lateral,
/// more than a groove filled in on a cap costs (`0.1`).
fn drawn(c: &Case<'_, '_>, m: &Model, solids: &[Handle<Solid>], volume: f64) {
    let cfg = nacre_tess::TessConfig::default();
    let t = nacre_tess::tessellate(m, &cfg)
        .unwrap_or_else(|e| panic!("{c}: the mesher refuses it: {e:?}"));
    assert_eq!(non_watertight(&t), 0, "{c}: the mesh leaks");
    let mut slack = 0.0;
    for &s in solids {
        let solid = m.solid(s);
        for &shell in std::iter::once(&solid.outer).chain(&solid.cavities) {
            for &f in &m.shell(shell).faces {
                if let nacre_geom::Surface::Cylinder(cyl) = m.surface_cache(m.face(f).surface) {
                    let turn = (cfg.max_angle_deg.to_radians() / 2.0).cos();
                    let sagitta = cfg.tol.min(cyl.radius() * (1.0 - turn));
                    let area = nacre_props::face_props(m, f)
                        .expect("a lateral's area")
                        .area
                        .abs();
                    slack += area * sagitta;
                }
            }
        }
    }
    let got = mesh_volume(&t);
    assert!(
        (got - volume).abs() <= slack + 1e-9,
        "{c}: mesh volume {got}, solids {volume}, slack {slack}"
    );
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

/// **Two secant walls, one entering `A` — the corner on the rim's seam (`(1, 0)`, `θ = 0` on
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

/// **The keyhole**: the box `[−1, 1] × [0, 3]` against the unit cylinder — its walls `x = ±1`
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

/// **A tangent wall on the edge of a flat through `A`'s own axis.** `A` is the half cylinder
/// over [`upper_half_disk`], so its flat `y = 0` meets its lateral in the ruling `(1, 0)` — `A`'s
/// own edge, where the lateral ends. `B`'s wall `x = 1` is tangent to `A`'s cylinder along that
/// edge, either running across it or ending on it (then `B`'s face `y = 0` lies on `A`'s flat).
/// Refused today, by name — and the name is the placement's, never the operation's or the vertex
/// order's: every boolean stops before the operation decides anything.
///
/// - **Across**: `TangentLineInAnotherPlane`, from the population gate, in either operand order,
///   at every height and in both frames. `A`'s flat holds the tangent line and `B`'s wall runs
///   across it, so the line is no edge both faces end on; around it six regions meet, which the
///   tangency's verdict does not read, and the gate refuses the pair before anything is arranged.
/// - **Ending**, either order and frame: the line *is* such an edge, so the gate passes it to the
///   structure, and the arrangement does not build that structure yet — `StraightAngle` inside the
///   span and past the top, `LabelConflict` on the base cap and past both caps.
///
/// The backstop cells ([`backstop`], `SuspectedDefect`, 48 booleans) and the across family's name
/// are pinned by variant, so a refusal that moves is news; the rest by class.
#[test]
fn a_tangent_wall_on_a_flats_edge_is_refused_today() {
    let across: &[[f64; 2]] = &[[-3.0, -2.0], [1.0, -2.0], [1.0, 0.5], [-3.0, 0.5]];
    let ending: &[[f64; 2]] = &[[-3.0, 0.0], [1.0, 0.0], [1.0, 0.5], [-3.0, 0.5]];
    let mut defects = 0;
    let ran = every_placement(
        upper_half_disk,
        &[("across the edge", across), ("ending on the edge", ending)],
        |c, m, got| {
            let pinned = match (c.placement.family, backstop(c)) {
                ("across the edge", _) => Some(RejectReason::TangentLineInAnotherPlane),
                (_, Some(want)) => {
                    defects += 1;
                    Some(want)
                }
                _ => None,
            };
            match pinned {
                Some(want) => match &got {
                    Err(BoolError::Rejected { reason, .. }) => {
                        assert_eq!(*reason, want, "{c}: the refusal moved");
                    }
                    other => panic!("{c}: the refusal moved — {other:?}"),
                },
                None => refused_today(c, m, got),
            }
        },
    );
    assert_eq!(ran, 2 * 2 * 2 * 4 * 6, "the family");
    assert_eq!(defects, 48, "the backstop refusals");
}

/// The half-cylinder family's backstop cells — where [`a_tangent_wall_on_a_flats_edge_is_refused_today`]
/// stops at `SuspectedDefect` today, and with which reason: the wall ending on the edge, on the
/// base cap or past both caps. Neither the operation nor the vertex order enters.
fn backstop(c: &Case<'_, '_>) -> Option<RejectReason> {
    use Height::*;
    let p = c.placement;
    match (p.family, p.height) {
        ("ending on the edge", OnBase | PastBoth) => Some(RejectReason::LabelConflict),
        _ => None,
    }
}

/// **Both walls entering `A` — an inward wedge**, its apex on the lateral at the seam and off
/// it. Near the line, `A − B` is two lobes meeting on it alone: where they join elsewhere (`B`
/// inside `A`'s span or reaching past one cap) that is one solid touching itself and is refused
/// (`NonManifoldResultEdge`); where `B`
/// cuts `A` through (past both caps) it is two solids. **It is never one body** — the lateral run
/// through the line as one face would hide the touch under a closed shell. Common, `B − A` and
/// the Fuse build one body; the volumes add up. The apex on the seam names the touch as it does
/// off it — the lateral's chart meets the apex as a corner its region walk passes twice, wherever
/// `θ = 0` falls, and the shell guard speaks. No placement answers `SuspectedDefect`.
#[test]
fn an_inward_wedge_on_a_lateral_builds_or_is_refused_by_name() {
    let on_the_seam: &[[f64; 2]] = &[[1.0, 0.0], [-3.0, 3.0], [-3.0, -3.0]];
    let off_the_seam: &[[f64; 2]] = &[[0.8, -0.6], [-3.0, 3.0], [-3.0, -3.0]];
    let all = every_six(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
    );
    assert_eq!(all.len(), 2 * 2 * 2 * 4, "the family");
    for (p, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match (i, r) {
                (2, Ok((n, _))) => {
                    assert_eq!(*n, 2, "{p}: A − B is two bodies or none");
                    assert_eq!(p.height, Height::PastBoth, "{p}: A − B built");
                }
                (_, Ok((n, _))) => assert_eq!(*n, 1, "{p}: [{i}]"),
                (2, Err(reason)) => assert_eq!(
                    *reason,
                    RejectReason::NonManifoldResultEdge,
                    "{p}: A − B refused as {reason:?}"
                ),
                (_, Err(reason)) => panic!("{p}: [{i}] {reason:?}"),
            }
        }
        identities(p, six);
    }
}

/// **A V-notch whose apex lies on the lateral** — `B` a box round `A` with a notch opening
/// away from the axis, its two walls meeting on a ruling of `A` at the seam and off it: the
/// inward wedge's walls turned inside out. `A` lies inside `B` but for that line, so `B − A` keeps
/// the material **outside** the cylinder on both sides of the line and none in the notch between
/// the walls — two lobes meeting on the line alone, one solid touching itself, refused
/// (`NonManifoldResultEdge`). The lateral keeps the outside here, so whether it runs on across the
/// line is asked of the chamber between the walls outside the cylinder (`cyl_chart::Chart::slit_at`);
/// read off the chamber inside, which is `A`'s and gone, the answer is right by accident, and
/// read off no chamber the lateral runs through the line as one face over a body touching itself
/// — `validate` clean and the volumes right. Everything else builds: `A ∩ B` is `A` over the
/// height the two share, `A − B` the rest of `A` (two bodies with `B` inside `A`'s span, none
/// with `B` past both caps).
#[test]
fn a_notch_whose_apex_is_on_a_lateral_builds_or_is_refused_by_name() {
    let on_the_seam: &[[f64; 2]] = &[
        [-2.0, -2.0],
        [3.0, -2.0],
        [3.0, -1.0],
        [1.0, 0.0],
        [3.0, 1.0],
        [3.0, 2.0],
        [-2.0, 2.0],
    ];
    // The same box turned by `(x, y) → (0.8x + 0.6y, −0.6x + 0.8y)`: the apex at `(0.8, −0.6)`.
    let off_the_seam: &[[f64; 2]] = &[
        [-2.8, -0.4],
        [1.2, -3.4],
        [1.8, -2.6],
        [0.8, -0.6],
        [3.0, -1.0],
        [3.6, -0.2],
        [-0.4, 2.8],
    ];
    let all = every_six(
        unit_disk,
        &[("on the seam", on_the_seam), ("off the seam", off_the_seam)],
    );
    assert_eq!(all.len(), 2 * 2 * 2 * 4, "the family");
    for (p, six) in &all {
        for (i, r) in six.iter().enumerate() {
            match (i, r) {
                (2, Ok((n, _))) => {
                    let lobes = match p.height {
                        Height::OnBase | Height::PastTop => 1,
                        Height::Inside => 2,
                        Height::PastBoth => 0,
                    };
                    assert_eq!(*n, lobes, "{p}: A − B");
                }
                (3, Err(reason)) => assert_eq!(
                    *reason,
                    RejectReason::NonManifoldResultEdge,
                    "{p}: B − A refused as {reason:?}"
                ),
                (3, Ok(got)) => panic!("{p}: B − A touches itself on the line, built: {got:?}"),
                (_, Ok((n, _))) => assert_eq!(*n, 1, "{p}: [{i}]"),
                (_, Err(reason)) => panic!("{p}: [{i}] {reason:?}"),
            }
        }
        identities(p, six);
        let common = six[4].as_ref().expect("built").1;
        let want = std::f64::consts::PI * p.height.shared();
        assert!(
            (common - want).abs() < 1e-9,
            "{p}: Common {common}, want {want}"
        );
    }
}

/// **A body touching itself along a ruling is refused; two bodies touching there build.** The
/// prism's two walls meet on a ruling of `A`'s cylinder and reach it from the side `A`'s lateral
/// does not keep — into the bore of a plate, or from outside a cylinder — so the lateral runs on
/// across the line and the line is an edge of the walls alone. Where that edge's body is the
/// lateral's own, joined elsewhere — a wedge across the bore into the plate beyond, a C-shaped
/// prism round the cylinder through it — the edge lies inside the solid's own face, and the result
/// is refused (`SelfTouchingResult`, the self-touch verdict's curved arm,
/// `combinatorics::ruling_meets_lateral`). Where it is another body — a short wedge inside the
/// bore, a corner touching the cylinder — the answer is the two bodies. Each refusal sits beside
/// a twin that must build, as the self-touch propositions do (`tests/probes/self_touch.rs`): a
/// check that refused every such touch would pass half of them.
#[test]
fn a_body_touching_itself_along_a_ruling_is_refused() {
    fn bored_plate() -> Profile2d {
        let square = nacre_ops::Ring2d::polygon_decimal(vec![
            p2(-2.0, -2.0),
            p2(2.0, -2.0),
            p2(2.0, 2.0),
            p2(-2.0, 2.0),
        ])
        .expect("a square");
        let bore = nacre_ops::Ring2d::circle(p2(0.0, 0.0), 1.0).expect("a circle");
        let mut plates = nacre_ops::from_paths(vec![square, bore]).expect("a plate with a bore");
        assert_eq!(plates.len(), 1, "one plate");
        plates.remove(0)
    }
    let tau = 2.0 * std::f64::consts::PI;
    type Case = (
        &'static str,
        fn() -> Profile2d,
        &'static [[f64; 2]],
        Option<(usize, f64)>,
    );
    let cases: [Case; 4] = [
        (
            "a wedge across the bore into the plate",
            bored_plate,
            &[[1.0, 0.0], [-3.0, 1.0], [-3.0, -1.0]],
            None,
        ),
        (
            "a wedge inside the bore",
            bored_plate,
            &[[1.0, 0.0], [0.0, 0.3], [0.0, -0.3]],
            Some((2, 32.0 - tau + 0.3)),
        ),
        (
            "a C round the cylinder and through it",
            unit_disk,
            &[
                [1.0, 0.0],
                [3.0, -2.0],
                [4.0, 0.0],
                [4.0, 3.0],
                [-0.6, 3.0],
                [-0.6, -2.0],
                [-0.3, -2.0],
                [-0.3, 2.5],
                [3.0, 2.5],
                [3.0, 2.0],
            ],
            None,
        ),
        (
            "a corner touching the cylinder",
            unit_disk,
            &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]],
            Some((2, tau + 6.0)),
        ),
    ];
    for tilted in [false, true] {
        for &(family, a, b, want) in &cases {
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
                    b,
                    0.5,
                    1.0,
                );
                let (x, y) = if swapped { (y, x) } else { (x, y) };
                let at = format!("{family}, tilted = {tilted}, swapped = {swapped}");
                match (boolean(&mut m, BoolKind::Fuse, x, y), want) {
                    (Ok(solids), Some((bodies, v))) => {
                        m.rebuild_adjacency();
                        let vs = nacre_validate::validate(&m);
                        assert!(vs.is_empty(), "{at}: {vs:?}");
                        assert_eq!(solids.len(), bodies, "{at}: bodies");
                        let got: f64 = solids.iter().map(|&s| volume(&m, s)).sum();
                        assert!((got - v).abs() < 1e-9, "{at}: volume {got}, want {v}");
                        let t = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                            .unwrap_or_else(|e| panic!("{at}: the mesher refuses it: {e:?}"));
                        assert_eq!(non_watertight(&t), 0, "{at}: the mesh leaks");
                    }
                    (Err(BoolError::Rejected { reason, .. }), None) => assert_eq!(
                        reason,
                        RejectReason::SelfTouchingResult,
                        "{at}: refused as {reason:?}"
                    ),
                    (got, _) => panic!("{at}: {got:?}"),
                }
            }
        }
    }
}

/// **A vertex lying on a cap's circle, and nothing else crossing it.** The prism's corner
/// sits on `A`'s rim — from inside the disk (a wedge whose apex alone reaches the lateral, on the
/// seam and off it) or from outside, with one side tangent there or with both sides leaving it (a
/// line contact) — so the circle is cut at exactly that one node: one arc from the node back to
/// itself. The planar arrangement takes that ring
/// (`combinatorics::ring_floor`, the self-loop OCCT and Parasolid keep for a closed edge), the
/// lateral's chart reads a rim of one node through that arc for every sector
/// (`Chart::arc_around`), and the region walk spells it as that one arc (`cyl_chart::regions`).
/// The wedge's union past `A`'s top leaves the cap **pinched** — the groove's tip on the rim —
/// and the closed edge sits in a loop of several: its half-edge runs the way its wall says, not
/// the way its one vertex could (`assembly::reconstruct`), and every reader takes it as the whole
/// turn (`Model::derive_edge_curve`); the mesher draws the pinch (`every_six`'s [`drawn`]).
/// Standing on `A`'s base, the wedge's floor merges into the cap: the coplanar merge leaves the
/// arc from the apex back to itself — a ring by the same floor — and the cleaning closes it as the
/// whole circle, so the union is `A`.
///
/// The cap ring an outside corner leaves is pinched at the node on the rim — the circle and the
/// prism's section meet there alone — and on either frame its winding is read where the circle
/// bulges past that node, so it is no refusal (`combinatorics::loop_winding`). Past `A`'s top,
/// `B`'s edge on the line of contact is stated as a ruling along `A`'s lateral and as a plane
/// segment above the rim — one line spelled two ways, which the straight-angle cleaning reads as
/// one (`combinatorics::SharedLines`), so `B` keeps no vertex at the rim point. Past both caps
/// the union is two bodies touching along that line, and the cylinder's every corner lies on it;
/// which body is inside which is read through its caps, each the arc from the corner back to
/// itself — a closed arc, whose span is the whole turn (`combinatorics::arc_span`).
///
/// What builds is checked against the input: `A ∩ B` is the wedge's section over the height the
/// two share (an outside corner meets `A` in a line only — empty), and [`identities`] carries
/// that to the union and the differences; the body counts are asserted beside. Two answers are
/// fixed whatever builds: the wedge's `A − B` never builds — the groove touches the lateral along
/// a line, a solid touching itself, `NonManifoldResultEdge` at every height — and no refusal is a
/// suspected defect. Everything else builds.
#[test]
fn a_vertex_on_a_cap_circle_builds_or_is_refused_by_name() {
    use std::collections::BTreeMap;
    let shapes: &[(&str, &[[f64; 2]])] = &[
        ("wedge on the seam", &[[1.0, 0.0], [0.5, 0.3], [0.5, -0.3]]),
        (
            "wedge off the seam",
            &[[0.8, -0.6], [0.5, -0.5], [0.4, -0.1]],
        ),
        ("tangent side", &[[1.0, 0.0], [1.0, -1.0], [2.0, -0.5]]),
        (
            "line contact",
            &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]],
        ),
    ];
    let tau = 2.0 * std::f64::consts::PI;
    let close = |got: f64, want: f64| (got - want).abs() < 1e-9;
    let mut tally: BTreeMap<(&str, usize, String, String), usize> = BTreeMap::new();
    for (p, six) in every_six(unit_disk, shapes) {
        identities(&p, &six);
        let (_, h) = p.height.span();
        let shared = p.height.shared();
        let section = p.b_volume() / h;
        let wedge = p.family.starts_with("wedge");
        // The pieces of the wedge outside `A`: none, the part above its top, or both ends.
        let outside = match p.height {
            Height::OnBase | Height::Inside => 0,
            Height::PastTop => 1,
            Height::PastBoth => 2,
        };
        for (i, r) in six.iter().enumerate() {
            match r {
                Ok((n, v)) => {
                    let (bodies, volume) = match (i, wedge) {
                        (0 | 1, true) => (1, tau + section * (h - shared)),
                        (0 | 1, false) => (2, tau + section * h),
                        (2, true) => panic!("{p}: [2] the wedge's `A − B` touches itself, built"),
                        (2, false) => (1, tau),
                        (3, true) => (outside, section * (h - shared)),
                        (3, false) => (1, section * h),
                        (_, true) => (1, section * shared),
                        (_, false) => (0, 0.0),
                    };
                    assert_eq!(*n, bodies, "{p}: [{i}] bodies");
                    assert!(close(*v, volume), "{p}: [{i}] volume {v}, want {volume}");
                }
                Err(reason) => {
                    assert_ne!(
                        reason.class(),
                        RejectClass::SuspectedDefect,
                        "{p}: [{i}] {reason:?}"
                    );
                    *tally
                        .entry((
                            p.family,
                            i,
                            format!("{:?}", p.height),
                            format!("{reason:?}"),
                        ))
                        .or_default() += 1;
                }
            }
        }
    }
    let mut want: BTreeMap<(&str, usize, String, String), usize> = BTreeMap::new();
    let mut pin = |family, ops: &[usize], heights: &[Height], reason: &str, each: usize| {
        for &i in ops {
            for h in heights {
                want.insert((family, i, format!("{h:?}"), reason.to_string()), each);
            }
        }
    };
    for wedge in ["wedge on the seam", "wedge off the seam"] {
        pin(wedge, &[2], &Height::ALL, "NonManifoldResultEdge", 4);
    }
    assert_eq!(tally, want);
}

/// **A tangent face that ends on the line at one height and crosses it at another.** `B` is
/// an L-shaped profile in the plane `x = −1` extruded to `x = 1`: its caps `x = ±1` are tangent
/// to `A`'s unit cylinder, each an L whose concave corner puts one stretch of its boundary on the
/// line `(±1, 0)` — beside `B`'s wall `y = 0`, a secant through `A`'s axis — while below that
/// stretch the line runs through the cap's interior. That is the face the tangency verdict's
/// `line_is_an_edge` keeps out (`!runs_through`: an edge on the line does not make the contact an
/// edge where the face runs across it) — and that guard is what speaks here: the row is no edge,
/// so the population gate refuses it (`TangentLineInAnotherPlane`), both vertex orders, all six
/// booleans, before the rulings split meets the cap's segment across the line.
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
                        RejectReason::TangentLineInAnotherPlane,
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

/// **A pinched cap as an operand** — the union a wedge with its apex on the rim leaves past the
/// cylinder's top ([`nacre_ops::fixtures::wedge_on_a_rim`]), whose top cap passes the apex twice,
/// taken into a second boolean against a box: far from the pinch on top, across it, and apart.
/// Every one builds, checked as the other families are — `validate`, the mesh, the volume
/// identities over the six booleans.
///
/// The cap's loop names the apex once, by its own definition (`combinatorics::loop_triples`);
/// the closed rim edge is the whole turn to every reader on the way — the span on its circle
/// (`combinatorics::theta_between`), and a cut line through the apex, which the arc crosses
/// once more at the other root (`combinatorics::ring_against_plane`); and the apex is pinned on a
/// line by the plane of its own that cuts that line (`combinatorics::pin_for`). The wedge's two
/// walls meet on a ruling of the cylinder, one station of its chart; taken from the box across the
/// pinch, the lateral keeps the material outside the cylinder, and whether it runs on across that
/// line is asked of the chamber between the walls out there (`cyl_chart::Chart::slit_at`) — the
/// one between the chords, inside, is the wedge's own and gone.
#[test]
fn a_pinched_cap_as_an_operand_builds() {
    use nacre_ops::fixtures::{cuboid, wedge_on_a_rim};
    let tools: &[(&str, [f64; 3], [f64; 3])] = &[
        ("box far on top", [-0.5, -0.25, 1.5], [0.0, 0.25, 3.0]),
        ("box across the pinch", [0.3, -0.5, 1.5], [1.5, 0.5, 3.0]),
        ("box apart", [5.0, 5.0, 0.0], [6.0, 6.0, 1.0]),
    ];
    for &(tool, lo, hi) in tools {
        let mut six: Vec<f64> = Vec::new();
        let (mut vw, mut vt) = (0.0, 0.0);
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            for swapped in [false, true] {
                let mut m = Model::new();
                let w = wedge_on_a_rim(&mut m, [[1.0, 0.0], [0.5, 0.3], [0.5, -0.3]], 0.5, 2.0);
                let t = cuboid(&mut m, Point3::from_array(lo), Point3::from_array(hi));
                m.rebuild_adjacency();
                (vw, vt) = (volume(&m, w), volume(&m, t));
                let (x, y) = if swapped { (t, w) } else { (w, t) };
                let at = format!("{tool}, {kind:?} swapped = {swapped}");
                let solids =
                    boolean(&mut m, kind, x, y).unwrap_or_else(|e| panic!("{at}: refused — {e:?}"));
                m.rebuild_adjacency();
                let vs = nacre_validate::validate(&m);
                assert!(vs.is_empty(), "{at}: {vs:?}");
                let t = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                    .unwrap_or_else(|e| panic!("{at}: the mesher refuses it: {e:?}"));
                assert_eq!(non_watertight(&t), 0, "{at}: the mesh leaks");
                six.push(solids.iter().map(|&s| volume(&m, s)).sum());
            }
        }
        // `W ∪ T`, `T ∪ W`, `W − T`, `T − W`, `W ∩ T`, `T ∩ W`.
        let close = |a: f64, b: f64| (a - b).abs() < 1e-9;
        for (x, y) in [(0, 1), (4, 5)] {
            assert!(
                close(six[x], six[y]),
                "{tool}: [{x}] {} vs [{y}] {}",
                six[x],
                six[y]
            );
        }
        for (i, want) in [(0, vw + vt), (2, vw), (3, vt)] {
            assert!(
                close(six[i] + six[4], want),
                "{tool}: [{i}] + Common = {}, want {want}",
                six[i] + six[4]
            );
        }
    }
}

/// **An outside corner as tall as the cylinder** — the tangent side and the line contact of
/// [`a_vertex_on_a_cap_circle_builds_or_is_refused_by_name`], `B` standing on `A`'s base and
/// reaching exactly its top, on both frames. The touch runs the cylinder's whole height, so every
/// corner of `A`'s body lies on the line of contact and both its caps are the closed arc from that
/// corner back to itself. All six build — valid and meshed: the union is two bodies touching along
/// the line, `A ∩ B` is empty, and the differences are the operands.
#[test]
fn an_outside_corner_as_tall_as_the_cylinder_builds() {
    let shapes: &[(&str, &[[f64; 2]])] = &[
        ("tangent side", &[[1.0, 0.0], [1.0, -1.0], [2.0, -0.5]]),
        (
            "line contact",
            &[[1.0, 0.0], [3.0, -2.0], [4.0, 0.0], [3.0, 2.0]],
        ),
    ];
    let tau = 2.0 * std::f64::consts::PI;
    for tilted in [false, true] {
        for &(family, quad) in shapes {
            let q = quad;
            let area = (0..q.len())
                .map(|k| {
                    let (p, r) = (q[k], q[(k + 1) % q.len()]);
                    p[0] * r[1] - r[0] * p[1]
                })
                .sum::<f64>()
                .abs()
                / 2.0;
            let b = area * 2.0;
            let mut i = 0;
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
                        unit_disk(),
                        quad,
                        0.0,
                        2.0,
                    );
                    let (x, y) = if swapped { (y, x) } else { (x, y) };
                    let at = format!("{family}, tilted = {tilted}, {kind:?} swapped = {swapped}");
                    let solids = boolean(&mut m, kind, x, y)
                        .unwrap_or_else(|e| panic!("{at}: refused — {e:?}"));
                    m.rebuild_adjacency();
                    let vs = nacre_validate::validate(&m);
                    assert!(vs.is_empty(), "{at}: {vs:?}");
                    let t = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                        .unwrap_or_else(|e| panic!("{at}: the mesher refuses it: {e:?}"));
                    assert_eq!(non_watertight(&t), 0, "{at}: the mesh leaks");
                    let v: f64 = solids.iter().map(|&s| volume(&m, s)).sum();
                    let (bodies, want) = match i {
                        0 | 1 => (2, tau + b),
                        2 => (1, tau),
                        3 => (1, b),
                        _ => (0, 0.0),
                    };
                    assert_eq!(solids.len(), bodies, "{at}: bodies");
                    assert!((v - want).abs() < 1e-9, "{at}: volume {v}, want {want}");
                    i += 1;
                }
            }
        }
    }
}
