//! **Which way a hole's arcs actually run, realized** — the independent oracle for
//! [`cycle_on_class`]'s winding rule.
//!
//! ★★★★★ **The rule is derived from signs and nothing downstream reads it yet.** A reversed arc
//! would put the band's answer inside the hole and the hole's on the band — the exactly-opposite
//! answer — and today's fixtures stop at `loop_winding` before `label_cells` could notice. So every
//! extent a hole carves is **realized here** and a lock judges it against the fixture's own
//! geometry. The kernel reads no coordinate to choose an arc; this reads one afterwards, to check.

use super::{Judge, NodeId, WorkingPlane, combinatorics};
use std::sync::Mutex;

/// The direction, from the circle's centre, of the **midpoint of the stated counter-clockwise
/// arc** — one entry per extent a cycle has carved out of a circle in this binary, whether the
/// face grazes there (the cycle's own arc) or is absent (a hole's interior, a panel's outside).
/// Beside it the cycle's kind and the cylinder's origin, so a reader can pick **its own**
/// fixture's holes out of a ledger every test in the binary writes to.
pub(crate) static MIDS: Mutex<Vec<Mid>> = Mutex::new(Vec::new());

#[derive(Clone, Copy, Debug)]
pub(crate) struct Mid {
    pub dir: [f64; 3],
    pub kind: combinatorics::CycleKind,
    pub origin: [f64; 3],
}

pub(crate) fn record(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    kind: combinatorics::CycleKind,
    arc: [NodeId; 2],
) {
    let (Some(pa), Some(pb)) = (
        combinatorics::pierce_point(jd, cyl, def, arc[0]),
        combinatorics::pierce_point(jd, cyl, def, arc[1]),
    ) else {
        return;
    };
    let o = def.origin().map(|x| x.to_f64());
    let raw = def.dir().map(|x| x.to_f64());
    let ml = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
    let m: [f64; 3] = core::array::from_fn(|i| raw[i] / ml);
    // The circle's centre is not `origin` — that is a point on the axis — so the axial part
    // comes off first.
    let perp = |p: [f64; 3]| -> [f64; 3] {
        let v: [f64; 3] = core::array::from_fn(|i| p[i] - o[i]);
        let h = v[0] * m[0] + v[1] * m[1] + v[2] * m[2];
        core::array::from_fn(|i| v[i] - h * m[i])
    };
    let va = perp(pa);
    let vb = perp(pb);
    let ra = (va[0] * va[0] + va[1] * va[1] + va[2] * va[2]).sqrt();
    let u: [f64; 3] = core::array::from_fn(|i| va[i] / ra);
    // `w` completes a right-handed frame with `u` about `m`, so +90° about `m` takes `u` to
    // `w` — which is the direction θ increases in.
    let w = [
        m[1] * u[2] - m[2] * u[1],
        m[2] * u[0] - m[0] * u[2],
        m[0] * u[1] - m[1] * u[0],
    ];
    // The counter-clockwise sweep from `a` to `b`, then half of it. Written as a rotation
    // rather than as `va + vb`, which vanishes when the arc is exactly a half — and the hole
    // this measures **is** exactly a half.
    let mut phi = (vb[0] * w[0] + vb[1] * w[1] + vb[2] * w[2])
        .atan2(vb[0] * u[0] + vb[1] * u[1] + vb[2] * u[2]);
    if phi <= 0.0 {
        phi += core::f64::consts::TAU;
    }
    let (c, s) = ((phi / 2.0).cos(), (phi / 2.0).sin());
    MIDS.lock()
        .expect("the probe's lock is never held across a panic")
        .push(Mid {
            dir: core::array::from_fn(|i| c * u[i] + s * w[i]),
            kind,
            origin: o,
        });
}
