//! **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement`] mounts it with `#[path]` as `disk_side_probe`, which is what keeps
//! `super::` here meaning that module — the file's place in the tree does not change what it
//! belongs to.

use super::*;
use crate::ledger::Ledger;

/// The sentence the disk-side check panics with — one spelling, so a commuting oracle `KNOWN`
/// entry can name the site by the constant the check prints.
pub(crate) const NOT_OWN_SOLID: &str = "a cut circle's disk side does not carry its own solid";

/// One entry per arc.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Row {
    /// The figure beside the arc named a side, so the rule was checked against the geometry here.
    pub(crate) checked: bool,
    /// Both cells are bounded, each on its own side — the strongest form of the check.
    pub(crate) both_bounded: bool,
    /// The class's `frame_sign` was `-1`, so the factor decided this arc.
    pub(crate) frame_negative: bool,
}

pub(crate) static ROWS: Ledger<Row> = Ledger::new();

/// Check the disk-side rule against the geometry, and record the row. `witness` is
/// [`witness`]'s answer for the arc: which half-edge's cell lies on the disk side, and whether
/// both cells are bounded. Where it speaks it must **agree with the derived rule** — the check that
/// keeps the derivation honest.
pub(crate) fn record(witness: Option<(bool, bool)>, rule_says_even: bool, frame_sign: i8) {
    if let Some((even_is_disk, _)) = witness {
        assert_eq!(
            even_is_disk, rule_says_even,
            "the disk-side rule and the figure beside the arc disagree (frame_sign {frame_sign})"
        );
    }
    ROWS.push(Row {
        checked: witness.is_some(),
        both_bounded: witness.is_some_and(|w| w.1),
        frame_negative: frame_sign < 0,
    });
}

/// A node's point, realized for the watcher below; `None` where it has no exact description.
fn node_point(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
) -> Option<[f64; 3]> {
    if let Some(p) = crate::combinatorics::node_coords_rat(jd, n) {
        return Some(p.map(|x| x.to_f64()));
    }
    let crate::combinatorics::NodeKind::Pierce { cyl, .. } = n.kind() else {
        return None;
    };
    let (line, s) =
        crate::combinatorics::pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?.into_exact()?;
    let (b, d, t) = (
        line.base().map(|x| x.to_f64()),
        line.dir().map(|x| x.to_f64()),
        s.to_f64_approx(),
    );
    Some(core::array::from_fn(|k| b[k] + t * d[k]))
}

/// A cylinder's section on one plane in `f64`, by the angle about the axis — the watcher's figure.
struct SectionFrame {
    o: [f64; 3],
    m: [f64; 3],
    e: [[f64; 3]; 2],
    r: f64,
    plane: [f64; 4],
}

impl SectionFrame {
    fn of(def: &nacre_topo::CylinderDef, plane: [f64; 4]) -> Self {
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let unit = |v: [f64; 3]| {
            let l = dot(v, v).sqrt();
            v.map(|x| x / l)
        };
        let m = unit(def.dir().map(|x| x.to_f64()));
        let r0 = def.ref_dir().map(|x| x.to_f64());
        let along = dot(r0, m);
        let e1 = unit(core::array::from_fn(|k| r0[k] - along * m[k]));
        let e2 = [
            m[1] * e1[2] - m[2] * e1[1],
            m[2] * e1[0] - m[0] * e1[2],
            m[0] * e1[1] - m[1] * e1[0],
        ];
        Self {
            o: def.origin().map(|x| x.to_f64()),
            m,
            e: [e1, e2],
            r: def.radius_f64(),
            plane,
        }
    }

    fn angle(&self, p: [f64; 3]) -> f64 {
        let q: [f64; 3] = core::array::from_fn(|k| p[k] - self.o[k]);
        let dot = |u: [f64; 3]| q[0] * u[0] + q[1] * u[1] + q[2] * u[2];
        dot(self.e[1]).atan2(dot(self.e[0]))
    }

    /// The point of the plane at angle `th` and distance `rho` off the axis.
    fn at(&self, th: f64, rho: f64) -> [f64; 3] {
        let (c, s) = (th.cos(), th.sin());
        let b: [f64; 3] =
            core::array::from_fn(|k| self.o[k] + rho * (c * self.e[0][k] + s * self.e[1][k]));
        let w = self.plane;
        let t = -(w[0] * b[0] + w[1] * b[1] + w[2] * b[2] + w[3])
            / (w[0] * self.m[0] + w[1] * self.m[1] + w[2] * self.m[2]);
        core::array::from_fn(|k| b[k] + t * self.m[k])
    }
}

/// **Which of an arc's two cells lies on the disk side, read beside the arc** — `Some((even, two))`:
/// whether the even half-edge's cell does, and whether both cells are bounded by contours that hold
/// opposite sides. `None` where the figure cannot say.
///
/// A cell's far corners cannot answer this: only a closed curve parts the plane, and where the
/// circle is drawn in part a cell reaches round the arc's end to the other side of it. So the two
/// points just off the arc's middle are asked of the arc's two contours, flattened. Contours that
/// hold opposite points are two bounded cells, each on its own point. Contours that hold the same
/// point are a bounded cell and the outer contour of its component, which encloses it: the bounded
/// one is on that point and the outer one's cell on the other. No sign of the walk or of the frame
/// is read, which is what lets this watch the rule that is written with them.
pub(super) fn witness(
    jd: &Judge<'_, WorkingPlane>,
    edges: &ClassEdges<'_>,
    cells: &[Cell],
    even: usize,
    wc: usize,
) -> Option<(bool, bool)> {
    let tau = std::f64::consts::TAU;
    let w = jd.world_rat(wc)?.map(|x| x.to_f64());
    let drop = (0..3).max_by(|&a, &b| w[a].abs().total_cmp(&w[b].abs()))?;
    let flat = |p: [f64; 3]| {
        let (u, v) = ((drop + 1) % 3, (drop + 2) % 3);
        [p[u], p[v]]
    };
    let span_of = |frame: &SectionFrame, i: usize| -> Option<(f64, f64)> {
        let [a, b] = edges.pieces[i].end;
        let ta = frame.angle(node_point(jd, edges.cyls, a)?);
        let tb = frame.angle(node_point(jd, edges.cyls, b)?);
        Some((
            ta,
            if a == b {
                tau
            } else {
                (tb - ta).rem_euclid(tau)
            },
        ))
    };
    let contour = |he: usize| -> Option<(usize, Vec<[f64; 2]>)> {
        let c = cells.iter().position(|q| q.half_edges.contains(&he))?;
        let mut out = Vec::new();
        for &h in &cells[c].half_edges {
            let HalfEdgeKind::Piece(i) = edges.kind(h) else {
                return None;
            };
            out.push(flat(node_point(jd, edges.cyls, edges.origin(h))?));
            if let OnCurve::Arc { cyl } = edges.pieces[i].on {
                let frame = SectionFrame::of(&edges.cyls.get(cyl)?.def, w);
                let (ta, span) = span_of(&frame, i)?;
                let steps = 2 * ((span / (tau / 360.0)).ceil() as usize).max(2);
                for k in 1..steps {
                    let f = k as f64 / steps as f64;
                    let f = if h % 2 == 0 { f } else { 1.0 - f };
                    out.push(flat(frame.at(ta + f * span, frame.r)));
                }
            }
        }
        Some((c, out))
    };
    let holds = |ring: &[[f64; 2]], q: [f64; 2]| {
        let mut inside = false;
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            if (a[1] > q[1]) != (b[1] > q[1])
                && q[0] < a[0] + (q[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0])
            {
                inside = !inside;
            }
        }
        inside
    };
    let area = |ring: &[[f64; 2]]| {
        let twice: f64 = (0..ring.len())
            .map(|k| {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                a[0] * b[1] - a[1] * b[0]
            })
            .sum();
        0.5 * twice.abs()
    };
    let i = even / 2;
    let OnCurve::Arc { cyl } = edges.pieces[i].on else {
        return None;
    };
    let frame = SectionFrame::of(&edges.cyls.get(cyl)?.def, w);
    let (ta, span) = span_of(&frame, i)?;
    let mid = ta + 0.5 * span;
    // Past the half-degree chords' sag (`9.5e-6·r`), or a neighbouring curve's chord could step over
    // the point.
    let delta = 2e-5 * frame.r;
    let (near, far) = (
        flat(frame.at(mid, frame.r - delta)),
        flat(frame.at(mid, frame.r + delta)),
    );
    let (ce, re) = contour(even)?;
    let (co, ro) = contour(even + 1)?;
    if ce == co {
        return None;
    }
    let (e_near, e_far) = (holds(&re, near), holds(&re, far));
    let (o_near, o_far) = (holds(&ro, near), holds(&ro, far));
    if e_near == e_far || o_near == o_far {
        return None;
    }
    if e_near != o_near {
        return Some((e_near, true));
    }
    let (ae, ao) = (area(&re), area(&ro));
    if (ae - ao).abs() <= 1e-9 * ae.max(ao) {
        return None;
    }
    Some((e_near != (ae > ao), false))
}
