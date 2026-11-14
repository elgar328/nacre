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
