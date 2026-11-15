//! Unit tests of the corner reader's statement doors.

use super::*;
use crate::realize::{push_plane_realized, push_plane_through_realized};
use crate::{OpOutput, Operation, apply};
use nacre_exact::{Angle, Axis, Rat};
use nacre_math::Point3;
use nacre_topo::{Motion, MotionNode, Orientation};

fn q(n: i128, d: i128) -> Rat {
    Rat::new(n, d).expect("a rational")
}

fn f64s(p: [Rat; 3]) -> Point3 {
    Point3::from_array(p.map(|r| r.to_f64()))
}

/// Push the plane through `points` under `motion`, its cache figure through `world` — the same
/// points as the motion places them, so the figure faces the way the truth does.
fn plane(
    m: &mut Model,
    points: [[Rat; 3]; 3],
    world: [Point3; 3],
    motion: Option<Handle<MotionNode>>,
) -> Handle<Surface> {
    let figure =
        nacre_geom::Plane::through_points(world[0], world[1], world[2]).expect("not collinear");
    push_plane_realized(m, figure, points, motion, Orientation::Forward).0
}

/// Points whose plane's canonical integers need more than `i128` (they run to about 146 bits).
fn wide_points() -> [[Rat; 3]; 3] {
    [
        [q(1, 1_000_003), q(2, 1_000_033), q(3, 1_000_037)],
        [q(5, 1_000_039), q(7, 1_000_081), q(1, 1_000_099)],
        [q(11, 1_000_117), q(13, 1_000_121), q(17, 1_000_133)],
    ]
}

/// `p` turned about the `z` axis through the origin, in `f64` — only the cache figure reads it.
fn turned_f64(p: [Rat; 3], deg: f64) -> Point3 {
    let (s, c) = deg.to_radians().sin_cos();
    let [x, y, z] = p.map(|r| r.to_f64());
    Point3::from_array([c * x - s * y, s * x + c * y, z])
}

fn turn_z(deg: i128) -> Motion {
    Motion::Rotate {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
    }
}

/// A row's answer: read, or the cause it failed by.
type Answer = Result<(), CornerFail>;

/// **A plane's narrow world coefficients fail by their cause** — geometry when the plane has no
/// name (moved or not) or its chain does not fold, width when the name is `Wide` or a folding
/// chain overflows, and geometry first when both hold. Each row asserts the premise that makes it
/// that cause before its answer is read; the answers are gathered so a wrong rule shows every row
/// it moves.
#[test]
fn a_planes_narrow_world_coefficients_fail_by_their_cause() {
    let mut m = Model::new();
    let zero = Rat::from_int(0);
    let one = Rat::from_int(1);
    let mut rows: Vec<(&str, Answer, Answer)> = Vec::new();
    let mut read = |m: &Model, label: &'static str, h: Handle<Surface>, want| {
        let got = world_narrow(m, h).map(|c| {
            assert_eq!(
                Some(&c),
                m.world_plane_name(h).as_ref().and_then(|n| n.narrow()),
                "{label}: the coefficients are the world name's narrow projection"
            );
        });
        rows.push((label, got, want));
    };

    // Narrow and unmoved: the world name itself.
    let flat = [[zero, zero, one], [one, zero, one], [zero, one, one]];
    let h = plane(&mut m, flat, flat.map(f64s), None);
    read(&m, "narrow", h, Ok(()));

    // Wide and unmoved: the world holds it, `Rat` does not.
    let wide = wide_points();
    let h = plane(&mut m, wide, wide.map(f64s), None);
    assert!(
        m.surface_motion(h).is_none() && m.surface_name(h).is_some_and(|n| n.narrow().is_none()),
        "premise: an unmoved plane with a wide name"
    );
    read(&m, "unmoved wide", h, Err(CornerFail::Width));

    // No name: one vertex of an unmoved box and two of a box turned 37° about `z` meet in no
    // one frame.
    let a = crate::fixtures::cuboid(&mut m, [0.0; 3], [1.0; 3]);
    let b = crate::fixtures::cuboid(&mut m, [3.0, 0.0, 0.0], [4.0, 1.0, 1.0]);
    let Ok(OpOutput::Transform { solid: b }) = apply(
        &mut m,
        &Operation::Transform {
            solid: b,
            isometry: nacre_exact::Isometry::rotation(nacre_exact::Rotation {
                axis: Axis::Z,
                pivot: [zero; 3],
                angle: Angle::from_deg(Rat::from_int(37)).expect("an angle"),
            }),
        },
    ) else {
        panic!("the turn");
    };
    let verts = |m: &Model, s: Handle<Solid>| -> Vec<Handle<nacre_topo::Vertex>> {
        let mut out = Vec::new();
        for &f in &m.shell(m.solid(s).outer).faces {
            for he in &m.face(f).outer.half_edges {
                let v = m.edge(he.edge).vertices[0];
                if !out.contains(&v) {
                    out.push(v);
                }
            }
        }
        out
    };
    let (va, vb) = (verts(&m, a), verts(&m, b));
    let mut through = [va[0], vb[0], vb[1]];
    through.sort_by_key(|v| v.index());
    let figure = nacre_geom::Plane::through_points(
        m.vertex_point(through[0]),
        m.vertex_point(through[1]),
        m.vertex_point(through[2]),
    )
    .expect("not collinear");
    let h = push_plane_through_realized(&mut m, figure, through, None, Orientation::Forward).0;
    assert!(
        m.surface_name(h).is_none(),
        "premise: a plane through vertices of two frames has no name"
    );
    read(&m, "no name", h, Err(CornerFail::Unstated));

    // The same statement carried by a translation: moved, on a chain that folds, and still no
    // name — the chain is not the cause.
    let nudge = m.push_motion(
        Motion::Translate {
            offset: [q(1, 3), zero, zero],
        },
        None,
    );
    let step = nacre_math::Vector3::from_array([1.0 / 3.0, 0.0, 0.0]);
    let figure = nacre_geom::Plane::through_points(
        m.vertex_point(through[0]) + step,
        m.vertex_point(through[1]) + step,
        m.vertex_point(through[2]) + step,
    )
    .expect("not collinear");
    let h =
        push_plane_through_realized(&mut m, figure, through, Some(nudge), Orientation::Forward).0;
    assert!(
        m.surface_name(h).is_none() && m.chain_folds(nudge),
        "premise: a moved plane with no name, on a chain that folds"
    );
    read(&m, "moved no name, folded", h, Err(CornerFail::Unstated));

    // A chain that does not fold: a narrow plane `z = 1` turned 30° about `z`.
    let turn = m.push_motion(turn_z(30), None);
    let h = plane(&mut m, flat, flat.map(|p| turned_f64(p, 30.0)), Some(turn));
    assert!(
        !m.chain_folds(turn) && m.surface_name(h).is_some_and(|n| n.narrow().is_some()),
        "premise: a narrow name on a chain that does not fold"
    );
    read(&m, "unfolded chain", h, Err(CornerFail::Unstated));

    // A wide name carried by a translation: the chain folds, the name has no narrow vessel.
    let shift = [q(1, 3), q(2, 7), q(5, 11)];
    let slide = m.push_motion(Motion::Translate { offset: shift }, None);
    let moved = |p: [Rat; 3]| {
        Point3::from_array(std::array::from_fn(|k| p[k].to_f64() + shift[k].to_f64()))
    };
    let h = plane(&mut m, wide, wide.map(moved), Some(slide));
    assert!(
        m.chain_folds(slide) && m.surface_name(h).is_some_and(|n| n.narrow().is_none()),
        "premise: a wide name on a chain that folds"
    );
    read(&m, "moved wide, folded", h, Err(CornerFail::Width));

    // A folding chain that overflows: `z = 1/p` (named `[0, 0, p, −1]`) lifted by `1/q`; the
    // canonical form clears the denominator `q` and needs `p·q`, past `i128`.
    let (p, qq): (i128, i128) = ((1 << 64) - 59, (1 << 64) - 83);
    assert!(
        p.checked_mul(qq).is_none(),
        "premise: p·q does not fit i128"
    );
    let low = [
        [zero, zero, q(1, p)],
        [one, zero, q(1, p)],
        [zero, one, q(1, p)],
    ];
    let lift = m.push_motion(
        Motion::Translate {
            offset: [zero, zero, q(1, qq)],
        },
        None,
    );
    let lifted = |pt: [Rat; 3]| {
        Point3::from_array([
            pt[0].to_f64(),
            pt[1].to_f64(),
            pt[2].to_f64() + 1.0 / qq as f64,
        ])
    };
    let h = plane(&mut m, low, low.map(lifted), Some(lift));
    let name = *m
        .surface_name(h)
        .and_then(|n| n.narrow())
        .expect("premise: a narrow name");
    assert!(
        m.chain_folds(lift) && m.chain_plane_coeffs(lift, name).is_none(),
        "premise: a chain that folds and overflows carrying the name"
    );
    read(&m, "folded overflow", h, Err(CornerFail::Width));

    // Both causes: a wide name on a chain that does not fold — geometry first.
    let h = plane(&mut m, wide, wide.map(|p| turned_f64(p, 30.0)), Some(turn));
    assert!(
        !m.chain_folds(turn) && m.surface_name(h).is_some_and(|n| n.narrow().is_none()),
        "premise: a wide name on a chain that does not fold"
    );
    read(&m, "moved wide, unfolded", h, Err(CornerFail::Unstated));

    let wrong: Vec<_> = rows.iter().filter(|(_, got, want)| got != want).collect();
    assert!(
        wrong.is_empty(),
        "rows answering by the wrong cause: {wrong:?}"
    );
}

/// **A vertex's reach holds its proven interval, and only a proven one is finite.**
///
/// The sieve drops a candidate before the exact test whenever reaches miss a box, so a reach
/// narrower than the truth's interval drops a real contact in silence. A `Realized` reach must
/// contain `coord ± bound` on every axis — an ordinary half-ulp, a zero bound (an exact
/// coordinate), and a radius far below an `f64`'s — and a cache that proves nothing must
/// keep every candidate.
#[test]
fn a_reach_holds_the_proven_interval_and_nothing_else_is_finite() {
    let coord = nacre_math::Point3::from_array([0.1, -3.0, 1.0e10]);
    let bound = [
        nacre_exact::Mag::of(0.1).times(nacre_exact::Mag::pow2(-53)),
        nacre_exact::Mag::ZERO,
        nacre_exact::Mag::pow2(-2000),
    ];
    let (lo, hi) = cache_reach(&nacre_topo::PointCache::Realized { coord, bound });
    for k in 0..3 {
        let (c, r) = (coord.as_array()[k], bound[k].upper_f64());
        assert!(
            lo[k] < c && c < hi[k],
            "axis {k}: the coordinate itself is inside"
        );
        assert!(
            c - lo[k] >= r && hi[k] - c >= r,
            "axis {k}: the bound {r:e} is inside"
        );
        assert!(
            lo[k].is_finite() && hi[k].is_finite(),
            "axis {k}: a proven reach is finite"
        );
    }
    for cache in [
        nacre_topo::PointCache::Ceiling { coord },
        nacre_topo::PointCache::Unrealized { coord },
    ] {
        let (lo, hi) = cache_reach(&cache);
        assert!(
            lo.iter().all(|x| *x == f64::NEG_INFINITY) && hi.iter().all(|x| *x == f64::INFINITY),
            "{cache:?} proves no bound, so its reach is the whole line"
        );
    }
}

/// **A vertex's rational reach moves an end `Rat` cannot hold outward, and holds the truth.**
///
/// An exact coordinate's reach ends at `±2⁻¹⁰²²` (the floor of the widening), and a coordinate
/// near `1e-30` ends past `Rat`'s denominators too: each end moves out to `±2⁻⁶⁴` or `0` — all four
/// ways, by the end's sign and side. An end `Rat` holds is that end, exactly. Too large an end, or a
/// cache that proves nothing, has no rational box.
#[test]
fn a_rational_reach_moves_an_unheld_end_outward() {
    let tiny = q(1, 1i128 << 64);
    let zero = Rat::from_int(0);
    let minus_tiny = zero.checked_sub(tiny).expect("a rational");
    let realized = |c: [f64; 3]| nacre_topo::PointCache::Realized {
        coord: Point3::from_array(c),
        bound: [nacre_exact::Mag::ZERO; 3],
    };
    let small = [0.0, 1.0e-30, -1.0e-30];
    let (lo, hi) = rational_reach(&realized(small)).expect("a finite reach");
    assert_eq!(
        (lo, hi),
        ([minus_tiny, zero, minus_tiny], [tiny, tiny, zero]),
        "every small end moves outward"
    );
    let held = [1.5, -3.0, 0.25];
    let cache = realized(held);
    let (lo, hi) = rational_reach(&cache).expect("a finite reach");
    let (flo, fhi) = cache_reach(&cache);
    for k in 0..3 {
        assert_eq!(
            Some(lo[k]),
            Rat::try_from_f64(flo[k]),
            "axis {k}: the low end held"
        );
        assert_eq!(
            Some(hi[k]),
            Rat::try_from_f64(fhi[k]),
            "axis {k}: the high end held"
        );
    }
    for c in [small, held] {
        let (lo, hi) = rational_reach(&realized(c)).expect("a finite reach");
        // Every end here is an `f64` exactly, so the comparison in `f64` is exact.
        for k in 0..3 {
            assert!(
                lo[k].to_f64() < c[k] && c[k] < hi[k].to_f64(),
                "axis {k}: the coordinate {} is inside",
                c[k]
            );
        }
    }
    assert_eq!(
        rational_reach(&realized([1.0e300, 0.0, 0.0])),
        None,
        "an end past `Rat` has no rational box"
    );
    let coord = Point3::from_array(held);
    for cache in [
        nacre_topo::PointCache::Ceiling { coord },
        nacre_topo::PointCache::Unrealized { coord },
    ] {
        assert_eq!(rational_reach(&cache), None, "{cache:?} proves no box");
    }
}

/// The faces of `solid` whose plane's world name is `coeffs`.
fn faces_on(m: &Model, solid: Handle<nacre_topo::Solid>, coeffs: &[Rat; 4]) -> Vec<Handle<Face>> {
    let s = m.solid(solid);
    std::iter::once(&s.outer)
        .chain(s.cavities.iter())
        .flat_map(|&sh| m.shell(sh).faces.iter().copied())
        .filter(|&fh| world_narrow(m, m.face(fh).surface).ok().as_ref() == Some(coeffs))
        .collect()
}

/// **A corner no exact spelling names is read on its proven box, by the one door.** The named wall
/// fused with a turned box ([`crate::tests::named_wall_with_unspelled_corners`]) carries, on its
/// faces on `x = 9` (the box runs through the wall's height, so it parts the face in two), corners
/// of that world plane, a turned wall of no world statement and a cap: no
/// frame holds the three, so the exact spelling stops at their shape. Each such corner comes back
/// as a box whose cause is that shape and whose carriers name `x = 9`, so the face's own class is
/// shown to hold it; the face's other corners keep their rational spelling.
#[test]
fn a_corner_no_exact_spelling_names_is_read_on_its_box() {
    let mut m = Model::new();
    let (_, tool) = crate::tests::named_wall_with_unspelled_corners(&mut m);
    let wall = [q(1, 1), q(0, 1), q(0, 1), q(-9, 1)];
    let faces = faces_on(&m, tool, &wall);
    assert_eq!(faces.len(), 2, "the box parts the wall's face on x = 9");
    let (mut boxed, mut rational) = (0, 0);
    for &fh in &faces {
        let face = m.face(fh);
        for he in &face.outer.half_edges {
            let c = corner_of(&m, face, &face.outer, he).expect("every corner of the face is read");
            assert!(c.on_plane(&wall), "the corner is shown on its face's plane");
            match c {
                Corner::Bounded { cause, .. } => {
                    assert_eq!(
                        cause,
                        CornerFail::Shape,
                        "no frame holds the corner's planes"
                    );
                    boxed += 1;
                }
                Corner::Rational(_) => rational += 1,
                Corner::Pierce(..) | Corner::Round { .. } => panic!("a straight face of planes"),
            }
        }
    }
    assert_eq!(boxed, 4, "the box's two walls cross the face at both caps");
    assert!(
        rational >= 2,
        "the wall's own corners stay exact ({rational})"
    );
}

/// A box corner about `centre` with reach `rho`, its carriers `carriers`.
fn bounded(centre: [Rat; 3], rho: Rat, carriers: [Option<[Rat; 4]>; 3]) -> Corner {
    Corner::Bounded {
        centre,
        rho,
        carriers,
        cause: CornerFail::Shape,
    }
}

/// **A box proves a side of the strip only when its disk is clear of it.** The plane `x = 0` holds
/// the axis of a radius-2 cylinder along `z`, so its strip is `|y| < 2`, with `Plus` on `−y`
/// (`n × m = −ŷ`). A box clear on either side says so — its centre may sit off the plane, since
/// its projection is what is asked — and a box touching the strip, or spanning one of its rulings,
/// proves no side: the true point may be anywhere in it.
#[test]
fn a_box_proves_a_strip_side_only_clear_of_the_strip() {
    use nacre_exact::StripSide::{Minus, Plus};
    let plane = [q(1, 1), q(0, 1), q(0, 1), q(0, 1)];
    let (o, m) = ([q(0, 1); 3], [q(0, 1), q(0, 1), q(1, 1)]);
    let r2 = nacre_exact::BigRat::from_rat(q(4, 1));
    let side = |c: [Rat; 3], rho: Rat| bounded(c, rho, [None; 3]).strip_side(&plane, &o, &m, &r2);
    let one = q(1, 1);
    assert_eq!(
        side([q(0, 1), q(5, 1), one], one),
        Some(Minus),
        "clear on +y"
    );
    assert_eq!(
        side([q(0, 1), q(-5, 1), one], one),
        Some(Plus),
        "clear on −y"
    );
    assert_eq!(
        side([q(1, 2), q(5, 1), one], one),
        Some(Minus),
        "off the plane, projected"
    );
    assert_eq!(
        side([q(0, 1), q(5, 2), one], one),
        None,
        "spanning the ruling y = 2"
    );
    assert_eq!(side([q(0, 1), q(0, 1), one], one), None, "inside the strip");
    assert_eq!(
        side([q(0, 1), q(3, 1), one], one),
        None,
        "touching the ruling"
    );
}

/// **A box reaches past a station when its ball does** — the axis question, about the same box:
/// centre at `z = 1`, reach `1`, so it sweeps `z ∈ [0, 2]`.
#[test]
fn a_box_reaches_past_a_station_by_its_ball() {
    use nacre_exact::Orient::{Negative, Positive};
    let (o, m) = ([q(0, 1); 3], [q(0, 1), q(0, 1), q(1, 1)]);
    let c = bounded([q(0, 1), q(5, 1), q(1, 1)], q(1, 1), [None; 3]);
    assert!(
        !c.reaches(&o, &m, q(3, 1), Positive),
        "z = 3 is past its reach"
    );
    assert!(c.reaches(&o, &m, q(3, 2), Positive), "z = 3/2 is within it");
    assert!(
        c.reaches(&o, &m, q(3, 1), Negative),
        "its centre is below z = 3"
    );
    assert!(
        !c.reaches(&o, &m, q(-1, 1), Negative),
        "nothing of it is below z = −1"
    );
    assert!(
        c.reaches(&o, &m, q(1, 2), Negative),
        "its reach is below z = 1/2"
    );
}

/// **A box is on a plane exactly when one of its corner's planes is that plane** — the box cannot
/// show a zero, so the corner's definition answers.
#[test]
fn a_box_is_on_a_plane_only_through_its_carriers() {
    let x0 = [q(1, 1), q(0, 1), q(0, 1), q(0, 1)];
    let x1 = [q(1, 1), q(0, 1), q(0, 1), q(-1, 1)];
    let c = bounded([q(0, 1); 3], q(1, 1), [None, Some(x0), None]);
    assert!(c.on_plane(&x0), "a carrier names x = 0");
    assert!(
        !c.on_plane(&x1),
        "no carrier names x = 1, though the box reaches it"
    );
}
