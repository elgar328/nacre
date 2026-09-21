use super::{Judge, NodeId, WorkingPlane, combinatorics};
use std::sync::Mutex;

#[derive(Clone, Copy, Debug)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct Hit {
    /// Read through `Debug` in the probe's messages.
    #[allow(dead_code)]
    pub point: [f64; 3],
    /// Distances to the class plane, the face plane, and the cylinder's surface.
    pub off: [f64; 3],
    /// `sign((x − o) · (m̂ × n_fc))` of the realized point — the f64 twin of
    /// [`crate::combinatorics::ruling_side`].
    pub side_f64: i8,
    /// The side the ring's edge carried.
    pub side: i8,
    /// `true` for a crossing on an arc, whose carried side is 0 — the ruling-side
    /// twin check is the rulings' alone.
    pub arc: bool,
}

pub(crate) static HITS: Mutex<Vec<Hit>> = Mutex::new(Vec::new());

#[allow(clippy::too_many_arguments)]
pub(crate) fn record(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    cyl: usize,
    wc: usize,
    fc: usize,
    side: i8,
    arc: bool,
    id: NodeId,
) {
    let Some(p) = combinatorics::pierce_point(jd, cyl, def, id) else {
        return;
    };
    let coeffs = |c: usize| -> Option<[f64; 4]> {
        combinatorics::class_coeffs_rat(jd, c).map(|w| w.map(|x| x.to_f64()))
    };
    let (Some(w), Some(v)) = (coeffs(wc), coeffs(fc)) else {
        return;
    };
    let plane_off = |w: [f64; 4]| -> f64 {
        let n = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        (w[0] * p[0] + w[1] * p[1] + w[2] * p[2] + w[3]).abs() / n
    };
    let o = def.origin().map(|x| x.to_f64());
    let raw = def.dir().map(|x| x.to_f64());
    let ml = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
    let m: [f64; 3] = core::array::from_fn(|i| raw[i] / ml);
    let d: [f64; 3] = core::array::from_fn(|i| p[i] - o[i]);
    let h = d[0] * m[0] + d[1] * m[1] + d[2] * m[2];
    let perp: [f64; 3] = core::array::from_fn(|i| d[i] - h * m[i]);
    let rho = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
    let cyl_off = (rho - def.radius_f64()).abs();
    let c = [
        m[1] * v[2] - m[2] * v[1],
        m[2] * v[0] - m[0] * v[2],
        m[0] * v[1] - m[1] * v[0],
    ];
    let dot = c[0] * d[0] + c[1] * d[1] + c[2] * d[2];
    let side_f64 = if dot > 0.0 {
        1
    } else if dot < 0.0 {
        -1
    } else {
        0
    };
    HITS.lock()
        .expect("the probe's lock is never held across a panic")
        .push(Hit {
            point: p,
            off: [plane_off(w), plane_off(v), cyl_off],
            side_f64,
            side,
            arc,
        });
}
