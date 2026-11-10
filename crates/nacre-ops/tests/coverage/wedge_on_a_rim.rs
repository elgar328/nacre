//! **A wedge whose two slanted faces meet on a bore's rim** — the scoreboard.
//!
//! `A` is the plate `[0,40]²×[0,10]` less a through bore of radius 5 about `(20, 10)`. `B` is the
//! `Common` of two prisms whose slanted faces `x = z + 13` and `y = 11 − z/2` meet on the line
//! through `P = (23, 6, 10)`, a point of the bore's top rim; the bottom-rim twin is the same wedge
//! mirrored in `z = 5`, meeting the bottom rim at `(23, 6, 0)`. At the rim's height the wedge is a
//! quadrant with its corner at `P`, and the rim's tangent there is `±(4, 3)`; on the other side of
//! the rim's plane the corner moves into the bore. So which quadrant the wedge takes decides
//! everything at `P`:
//!
//! - **outside** (`x ≥ z + 13`, `y ≤ 11 − z/2`): the quadrant lies outside the bore but for `P`.
//!   `A − B` keeps the plate on both sides of the quadrant, and below the rim the bore and the
//!   wedge are one hole: two sectors of material meeting at `P` alone — a pinch,
//!   `NonManifoldVertex`. `A ∪ B` fills that hole's corner and its link at `P` is one circle: one
//!   body. `B − A` is the wedge's part beyond the plate and its part in the bore, meeting at `P`
//!   alone: two bodies.
//! - **inside** (`x ≤ z + 13`, `y ≥ 11 − z/2`): the quadrant lies in the bore, so the wedge touches
//!   the plate at `P` alone — `A ∪ B` is the pinch.
//! - **across** (the two others): the quadrant crosses the rim, and every boolean builds.
//!
//! Each wedge either ends inside the plate (`z ∈ [5, 15]` on the top rim's side) or runs through
//! it (`z ∈ [−5, 15]`), where `B − A` also splits at the plate's far face unless the bore joins the
//! two parts. Every cell is one boolean on a fresh model: built ones are held to their body count,
//! `validate` and the mesh, refused ones to their name, and the volumes to the boolean identities
//! with `A = 16000 − 250π` and `B` in closed form — its section is the product of two widths
//! linear in `z` over the one height range both prisms share, so Simpson's rule is exact.
//!
//! **The bridged family** stands a bar across the bore's far side on `A`'s top face: its floor cuts
//! the top rim at `(16, 13)` and `(24, 13)` while the lateral runs on under it, so the lateral's
//! boundary cycle through `P` carries a whole turn of a rim cut at three nodes — and with the
//! bore's seam along `+y`, `(16, 13)` comes before `P` in the rim's order, so that turn does not
//! start at the rim's first node.
//!
//! **One shape answers wrong today**: the outside wedge ending inside the plate on the top rim,
//! fused — one body, `validate` clean, its volumes right, and the mesher refuses its lateral
//! (`DegenerateRing`; `SelfTouchingBoundary` with the bridged family's seam). The bottom-rim twin
//! meshes. todo «띠 옆면의 링이 림 위 한 점을 두 번 지나는 `Fuse` 는 메셔가 거절한다».

use crate::common::*;
use crate::stated::*;
use nacre_exact::{Axis, Isometry, Rat};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpOutput, Operation, RejectReason, SketchFrame, apply, boolean,
};
use nacre_store::Handle;
use nacre_tess::TessError;
use nacre_topo::{Model, Solid};

#[derive(Clone, Copy, Debug, PartialEq)]
enum Rim {
    Top,
    Bottom,
}

/// Which quadrant of the rim's plane the wedge takes at `P`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Quadrant {
    Outside,
    Inside,
    /// `x ≥ z + 13`, `y ≥ 11 − z/2`.
    AcrossUp,
    /// `x ≤ z + 13`, `y ≤ 11 − z/2`.
    AcrossDown,
}

/// The wedge's far end, in the top rim's frame: ending inside the plate at `z = 5`, or running
/// through it to `z = −5`.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Depth {
    Inside,
    Through,
}

#[derive(Clone, Copy, Debug)]
struct Case {
    rim: Rim,
    quadrant: Quadrant,
    depth: Depth,
    bridged: bool,
}

/// What one cell answers.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Want {
    /// This many bodies, `validate` clean, meshed.
    Built(usize),
    /// Refused by this name.
    Refused(RejectReason),
    /// This many bodies, `validate` clean, and the mesher refuses them by this name.
    Undrawn(usize, TessError),
}

/// A decimal the sketch lifts exactly.
fn q(x: f64) -> Rat {
    Rat::from_decimal(x).expect("a short decimal")
}

/// A polygon extruded `depth` along the frame's normal, then moved by `shift`. The corners are
/// put counter-clockwise first.
fn prism(
    m: &mut Model,
    axis: Axis,
    mut pts: Vec<[f64; 2]>,
    depth: f64,
    shift: [f64; 3],
) -> Handle<Solid> {
    let n = pts.len();
    let twice_area: f64 = (0..n)
        .map(|i| {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            a[0] * b[1] - b[0] * a[1]
        })
        .sum();
    if twice_area < 0.0 {
        pts.reverse();
    }
    let p2 = |p: [f64; 2]| Point2::from_array(p);
    let edges = (0..n)
        .map(|i| line(p2(pts[i]), p2(pts[(i + 1) % n])))
        .collect();
    let profile = stated(edges).expect("a polygon").remove(0);
    let frame = SketchFrame::world(m, axis);
    let Ok(OpOutput::Extrude { solid, .. }) = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: depth,
        },
    ) else {
        panic!("the prism extrudes")
    };
    let Ok(OpOutput::Transform { solid }) = apply(
        m,
        &Operation::Transform {
            solid,
            isometry: Isometry::translation(shift.map(q)),
        },
    ) else {
        panic!("the prism moves")
    };
    solid
}

/// The wedge's two widths at height `t` of the top rim's frame — `[x_lo, x_hi]` and `[y_lo,
/// y_hi]`, each linear in `t`.
fn widths(quadrant: Quadrant, t: f64) -> ([f64; 2], [f64; 2]) {
    let (x_out, y_out) = match quadrant {
        Quadrant::Outside => (true, true),
        Quadrant::Inside => (false, false),
        Quadrant::AcrossUp => (true, false),
        Quadrant::AcrossDown => (false, true),
    };
    let x = if x_out {
        [t + 13.0, 34.0]
    } else {
        [6.0, t + 13.0]
    };
    let y = if y_out {
        [1.0, 11.0 - t / 2.0]
    } else {
        [11.0 - t / 2.0, 19.0]
    };
    (x, y)
}

/// The wedge's height range in the top rim's frame.
fn span(depth: Depth) -> [f64; 2] {
    match depth {
        Depth::Inside => [5.0, 15.0],
        Depth::Through => [-5.0, 15.0],
    }
}

fn wedge(m: &mut Model, c: Case) -> Handle<Solid> {
    let z = |t: f64| match c.rim {
        Rim::Top => t,
        Rim::Bottom => 10.0 - t,
    };
    let [a, b] = span(c.depth);
    let (xa, ya) = widths(c.quadrant, a);
    let (xb, yb) = widths(c.quadrant, b);
    // In the world `zx` frame `u = z`, `v = x`, extruded along `+y`; in `yz`, `u = y`, `v = z`,
    // along `+x`.
    let along_y = prism(
        m,
        Axis::Y,
        vec![[z(a), xa[0]], [z(a), xa[1]], [z(b), xb[1]], [z(b), xb[0]]],
        20.0,
        [0.0, 0.5, 0.0],
    );
    let along_x = prism(
        m,
        Axis::X,
        vec![[ya[0], z(a)], [ya[1], z(a)], [yb[1], z(b)], [yb[0], z(b)]],
        33.0,
        [3.0, 0.0, 0.0],
    );
    boolean(m, BoolKind::Common, along_y, along_x).expect("the wedge")[0]
}

/// `B` in closed form: Simpson's rule over its section, a quadratic in the height.
fn wedge_volume(c: Case) -> f64 {
    let [a, b] = span(c.depth);
    let section = |t: f64| {
        let (x, y) = widths(c.quadrant, t);
        (x[1] - x[0]) * (y[1] - y[0])
    };
    (b - a) / 6.0 * (section(a) + 4.0 * section((a + b) / 2.0) + section(b))
}

fn plate(m: &mut Model, bridged: bool) -> Handle<Solid> {
    let slab = nacre_ops::fixtures::cuboid(m, [0.0; 3], [40.0, 40.0, 10.0]);
    let seam = if bridged {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, -1.0, 0.0]
    };
    let bore = nacre_ops::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([20.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array(seam),
        5.0,
        12.0,
    )
    .solid;
    let mut a = boolean(m, BoolKind::Cut, slab, bore).expect("the bored plate")[0];
    if bridged {
        // A bridge over the bore's far side, standing on the top face: its floor `z = 10` cuts
        // the top rim at `(16, 13)` and `(24, 13)` while the lateral runs on under it.
        let bridge = nacre_ops::fixtures::cuboid(m, [10.0, 13.0, 10.0], [30.0, 18.0, 12.0]);
        a = boolean(m, BoolKind::Fuse, a, bridge).expect("the bridge")[0];
    }
    a
}

/// What each of `[A ∪ B, A − B, A ∩ B, B ∪ A, B − A, B ∩ A]` answers, derived in the module doc.
fn expected(c: Case) -> [Want; 6] {
    use Want::{Built, Refused, Undrawn};
    let pinch = Refused(RejectReason::NonManifoldVertex);
    match c.quadrant {
        Quadrant::Outside => {
            // The mesher's refusal (module doc): its name follows where the bore's seam falls.
            let fuse = match (c.rim, c.depth, c.bridged) {
                (Rim::Top, Depth::Inside, false) => Undrawn(1, TessError::DegenerateRing),
                (Rim::Top, Depth::Inside, true) => Undrawn(1, TessError::SelfTouchingBoundary),
                _ => Built(1),
            };
            [fuse, pinch, Built(1), fuse, Built(2), Built(1)]
        }
        Quadrant::Inside => {
            // Run through, the wedge's part below the plate is apart from the part above it.
            let rest = if c.depth == Depth::Through { 2 } else { 1 };
            [pinch, Built(1), Built(1), pinch, Built(rest), Built(1)]
        }
        // Across the rim the wedge's parts above and below the plate are joined through the bore
        // when it leans into it on the far side too, and apart when it does not.
        Quadrant::AcrossUp => [Built(1); 6],
        Quadrant::AcrossDown => {
            let rest = if c.depth == Depth::Through { 2 } else { 1 };
            let mut w = [Built(1); 6];
            w[4] = Built(rest);
            w
        }
    }
}

fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    for rim in [Rim::Top, Rim::Bottom] {
        for quadrant in [
            Quadrant::Outside,
            Quadrant::Inside,
            Quadrant::AcrossUp,
            Quadrant::AcrossDown,
        ] {
            for depth in [Depth::Inside, Depth::Through] {
                out.push(Case {
                    rim,
                    quadrant,
                    depth,
                    bridged: false,
                });
            }
        }
    }
    for depth in [Depth::Inside, Depth::Through] {
        out.push(Case {
            rim: Rim::Top,
            quadrant: Quadrant::Outside,
            depth,
            bridged: true,
        });
    }
    out
}

/// **Every cell, every family** — the mismatches gathered and reported together, so a change
/// shows every cell it moves.
#[test]
fn a_wedge_meeting_a_rim_at_a_point_builds_or_is_refused_by_name() {
    let mut wrong: Vec<String> = Vec::new();
    for c in cases() {
        let want = expected(c);
        // Volumes of the six cells that build, for the identities.
        let mut vol: [Option<f64>; 6] = [None; 6];
        let mut a_volume = 0.0;
        for (i, &cell) in want.iter().enumerate() {
            let (kind, swap) = (
                [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common][i % 3],
                i >= 3,
            );
            let mut m = Model::new();
            let a = plate(&mut m, c.bridged);
            a_volume = volume(&m, a);
            let b = wedge(&mut m, c);
            let b_volume = volume(&m, b);
            if (b_volume - wedge_volume(c)).abs() > 1e-9 * b_volume {
                wrong.push(format!("{c:?}: B is {b_volume}, not {}", wedge_volume(c)));
            }
            let (x, y) = if swap { (b, a) } else { (a, b) };
            let got = match boolean(&mut m, kind, x, y) {
                Err(BoolError::Rejected { reason, .. }) => Want::Refused(reason),
                Err(e) => panic!("{c:?} {kind:?} swap={swap}: {e:?}"),
                Ok(out) => {
                    let v = nacre_validate::validate(&m);
                    if !v.is_empty() {
                        wrong.push(format!("{c:?} {kind:?} swap={swap}: invalid {v:?}"));
                    }
                    vol[i] = Some(out.iter().map(|&s| volume(&m, s)).sum());
                    match nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()) {
                        Ok(_) => Want::Built(out.len()),
                        Err(e) => Want::Undrawn(out.len(), e),
                    }
                }
            };
            if got != cell {
                wrong.push(format!(
                    "{c:?} {kind:?} swap={swap}: got {got:?}, want {cell:?}"
                ));
            }
        }
        if !c.bridged
            && (a_volume - (16000.0 - 250.0 * std::f64::consts::PI)).abs() > 1e-9 * 16000.0
        {
            wrong.push(format!("{c:?}: A is {a_volume}"));
        }
        let b_volume = wedge_volume(c);
        let near = |p: f64, q: f64| (p - q).abs() <= 1e-9 * a_volume;
        let mut identity = |what: &str, lhs: Option<f64>, rhs: Option<f64>| {
            if let (Some(l), Some(r)) = (lhs, rhs)
                && !near(l, r)
            {
                wrong.push(format!("{c:?}: {what}: {l} vs {r}"));
            }
        };
        let sum = |p: Option<f64>, q: Option<f64>| Some(p? + q?);
        for o in [0, 3] {
            identity(
                "A ∪ B + A ∩ B = A + B",
                sum(vol[o], vol[o + 2]),
                Some(a_volume + b_volume),
            );
        }
        identity("A − B = A − A ∩ B", vol[1], vol[2].map(|v| a_volume - v));
        identity("B − A = B − A ∩ B", vol[4], vol[5].map(|v| b_volume - v));
    }
    assert!(
        wrong.is_empty(),
        "{} cells:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}
