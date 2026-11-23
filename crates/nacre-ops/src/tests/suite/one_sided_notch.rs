//! **A lateral whose boundary runs along one of its rulings and turns back to one side** — a notch
//! cut into the unit rod, its flat wall through the line `(1, 0)`, against a wedge whose apex runs
//! along that line. Over the notch the lateral ends on the line and the wall holds it; above and
//! below the lateral runs on across it. A plane through the line reads a ruling piece on each side
//! of the run and none along it (`rulings::cycle_on_wall`).
//!
//! Every body is read in polar form about the apex: a section of the rod is the rays `θ` from
//! `(1, 0)` measured off the inward direction `(−1, 0)`, each reaching the circle at `2 cos θ`, so
//! the rays `[p, q]` sweep the area `F(q) − F(p)`, `F(θ) = θ + sin θ cos θ`. A notch whose wall holds
//! the line keeps an interval of those rays, and the volumes are closed forms slab by slab.
//!
//! A cell is built or refused as a result touching itself along the line. The refusals and the
//! bodies built are counted per family, the counts read off the slabs' cross-sections apart from
//! the kernel: a result touches itself where two of its sectors at the apex are one body.

use super::scoreboard::{OPS, Op, run};
use super::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

/// The area the rays `[p, q]` from the apex sweep in the unit disk.
fn swept(p: f64, q: f64) -> f64 {
    let f = |t: f64| t + t.sin() * t.cos();
    if p < q { f(q) - f(p) } else { 0.0 }
}

/// A direction's angle off the inward direction, in `(−π, π]`.
fn off_inward(v: [f64; 2]) -> f64 {
    let t = v[1].atan2(v[0]) - PI;
    if t <= -PI { t + 2.0 * PI } else { t }
}

/// The prism over the polygon `pts` from height `z` for `h`.
pub(super) fn prism(m: &mut Model, pts: &[[f64; 2]], z: f64, h: f64) -> Handle<Solid> {
    let frame = datum_frame(
        m,
        SketchPlane::from_axes(
            Point3::from_array([0.0, 0.0, z]),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
        ),
    );
    let profile =
        Profile2d::polygon(pts.iter().map(|q| p2(q[0], q[1])).collect()).expect("a polygon");
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    )
    .expect("the prism extrudes") else {
        unreachable!()
    };
    solid
}

/// The unit rod about the `z` axis over `z ∈ [z0, z1]`.
fn rod(m: &mut Model, z0: f64, z1: f64) -> Handle<Solid> {
    crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([0.0, 0.0, z0]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
        z1 - z0,
    )
    .solid
}

fn cut(m: &mut Model, a: Handle<Solid>, pts: &[[f64; 2]], z: f64, h: f64) -> Handle<Solid> {
    let t = prism(m, pts, z, h);
    let out = crate::boolean(m, BoolKind::Cut, a, t).expect("the notch cuts");
    assert_eq!(out.len(), 1, "the notched rod is one body");
    out[0]
}

/// A wedge with its apex on the line and its walls along `a` and `b`, counterclockwise.
#[derive(Clone, Copy)]
struct Wedge {
    a: [f64; 2],
    b: [f64; 2],
}

impl Wedge {
    /// Its cross-section, long enough to pass every rod's section.
    fn section(self) -> [[f64; 2]; 3] {
        let at = |v: [f64; 2]| [1.0 + 8.0 * v[0], 8.0 * v[1]];
        [[1.0, 0.0], at(self.a), at(self.b)]
    }

    fn area(self) -> f64 {
        32.0 * (self.a[0] * self.b[1] - self.a[1] * self.b[0])
    }

    /// What it takes from a section holding the rays `[lo, hi]`.
    fn in_rays(self, [lo, hi]: [f64; 2]) -> f64 {
        swept(off_inward(self.a).max(lo), off_inward(self.b).min(hi))
    }

    /// The same wedge reflected in `y = 0`.
    fn mirrored(self) -> Self {
        let flip = |v: [f64; 2]| [v[0], -v[1]];
        Wedge {
            a: flip(self.b),
            b: flip(self.a),
        }
    }
}

const OUTSIDE: Wedge = Wedge {
    a: [1.0, 1.0],
    b: [0.0, 1.0],
};
const INTO_ABOVE: Wedge = Wedge {
    a: [0.0, 1.0],
    b: [-1.0, 1.0],
};
const ACROSS: Wedge = Wedge {
    a: [-1.0, 1.0],
    b: [-1.0, -1.0],
};
const INTO_BELOW: Wedge = Wedge {
    a: [-1.0, -1.0],
    b: [0.0, -1.0],
};
const OVER_TOP_TANGENT: Wedge = Wedge {
    a: [1.0, 1.0],
    b: [-1.0, 1.0],
};
const OVER_BOTTOM_TANGENT: Wedge = Wedge {
    a: [-1.0, -1.0],
    b: [1.0, -1.0],
};

const WEDGES: [Wedge; 6] = [
    OUTSIDE,
    INTO_ABOVE,
    ACROSS,
    INTO_BELOW,
    OVER_TOP_TANGENT,
    OVER_BOTTOM_TANGENT,
];

/// Over a notch's end, inside the first notch, between the notches, and over both.
const SPANS: [[f64; 2]; 4] = [[0.5, 1.5], [1.25, 1.75], [2.25, 2.75], [1.5, 3.5]];

/// One boolean against its closed-form volume: the number of bodies, or `None` for a result that
/// touches itself along the line, the one refusal these bodies may answer with.
pub(super) fn scored(
    label: &str,
    m: &mut Model,
    op: Op,
    a: Handle<Solid>,
    b: Handle<Solid>,
    want: f64,
) -> Option<usize> {
    match run(m, op, a, b) {
        Ok(out) => {
            let bad = nacre_validate::validate(m);
            assert!(bad.is_empty(), "{label}: {bad:?}");
            let v: f64 = out
                .iter()
                .map(|&s| nacre_props::mass_props(m, s).expect("props").volume)
                .sum();
            assert!((v - want).abs() < 1e-9, "{label}: {v}, want {want}");
            mesh_covers_faces(label, m, &out);
            Some(out.len())
        }
        Err(BoolError::Rejected {
            reason: RejectReason::NonManifoldResultEdge | RejectReason::SelfTouchingResult,
            ..
        }) => None,
        Err(e) => panic!("{label}: {e:?}"),
    }
}

/// The four volumes of `A ∪ B`, `A − B`, `B − A`, `A ∩ B`.
pub(super) fn four(a: f64, b: f64, common: f64) -> [f64; 4] {
    [a + b - common, a - common, b - common, common]
}

const BELOW_Y0: &[[f64; 2]] = &[[-3.0, -3.0], [3.0, -3.0], [3.0, 0.0], [-3.0, 0.0]];
const ABOVE_Y0: &[[f64; 2]] = &[[-3.0, 0.0], [3.0, 0.0], [3.0, 3.0], [-3.0, 3.0]];
/// The cap `x − y > 1`.
const CAP_LOW: &[[f64; 2]] = &[[-2.0, -3.0], [4.0, -3.0], [4.0, 3.0]];
/// The cap `x + y > 1`.
const CAP_HIGH: &[[f64; 2]] = &[[4.0, -3.0], [4.0, 3.0], [-2.0, 3.0]];
/// The cap `x + 2y > 1`.
const CAP_SECOND: &[[f64; 2]] = &[[5.0, -2.0], [5.0, 5.0], [-3.0, 5.0], [-3.0, 2.0]];

/// `−atan(1/2)`, the ray of the wall `x + 2y = 1` on the disk's side.
fn second_wall() -> f64 {
    -(0.5f64.atan())
}

/// A body's stretches of height and the rays each keeps: `kept` over the first notch, the second
/// notch's where there is one, the half turn elsewhere.
fn slabs(kept: [f64; 2], second: bool) -> Vec<([f64; 2], [f64; 2])> {
    let whole = [-FRAC_PI_2, FRAC_PI_2];
    let mut out = vec![([0.0, 1.0], whole), ([1.0, 2.0], kept)];
    if second {
        out.extend([
            ([2.0, 3.0], whole),
            ([3.0, 4.0], [second_wall(), FRAC_PI_2]),
            ([4.0, 5.0], whole),
        ]);
    } else {
        out.push(([2.0, 5.0], whole));
    }
    out
}

/// **A notch whose wall holds the line, against every wedge on it.** Past the run the lateral
/// straddles a plane through the line and along it the lateral ends there, so the plane states the
/// run and the ruling piece stops at its two ends. `kept` is the rays over the first notch, `second`
/// whether the body has the second, and `counts` the refused cells and the bodies of the built ones.
fn a_wedge_on_the_line(
    build: fn(&mut Model) -> Handle<Solid>,
    kept: [f64; 2],
    second: bool,
    counts: (usize, usize),
) {
    let parts = slabs(kept, second);
    let va: f64 = parts
        .iter()
        .map(|(z, r)| (z[1] - z[0]) * swept(r[0], r[1]))
        .sum();
    let (mut refused, mut bodies) = (0, 0);
    for w in WEDGES {
        for span in SPANS {
            let common: f64 = parts
                .iter()
                .map(|(z, r)| (span[1].min(z[1]) - span[0].max(z[0])).max(0.0) * w.in_rays(*r))
                .sum();
            let wants = four(va, w.area() * (span[1] - span[0]), common);
            for (op, want) in OPS.into_iter().zip(wants) {
                let mut m = Model::new();
                let a = build(&mut m);
                let b = prism(&mut m, &w.section(), span[0], span[1] - span[0]);
                let label = format!("{kept:?}, {:?}..{:?}, {span:?}, {op:?}", w.a, w.b);
                match scored(&label, &mut m, op, a, b, want) {
                    Some(n) => bodies += n,
                    None => refused += 1,
                }
            }
        }
    }
    assert_eq!((refused, bodies), counts);
}

/// The rod over `[0, 5]` less `tool` over `[1, 2]`.
fn notched(m: &mut Model, tool: &[[f64; 2]]) -> Handle<Solid> {
    let r = rod(m, 0.0, 5.0);
    cut(m, r, tool, 1.0, 1.0)
}

#[test]
fn a_wedge_on_the_line_of_the_half_notch_below_builds() {
    a_wedge_on_the_line(|m| notched(m, BELOW_Y0), [-FRAC_PI_2, 0.0], false, (7, 89));
}

#[test]
fn a_wedge_on_the_line_of_the_half_notch_above_builds() {
    a_wedge_on_the_line(|m| notched(m, ABOVE_Y0), [0.0, FRAC_PI_2], false, (7, 89));
}

#[test]
fn a_wedge_on_the_line_of_the_low_cap_builds() {
    a_wedge_on_the_line(
        |m| notched(m, CAP_LOW),
        [-FRAC_PI_2, FRAC_PI_4],
        false,
        (3, 91),
    );
}

#[test]
fn a_wedge_on_the_line_of_the_high_cap_builds() {
    a_wedge_on_the_line(
        |m| notched(m, CAP_HIGH),
        [-FRAC_PI_4, FRAC_PI_2],
        false,
        (3, 91),
    );
}

#[test]
fn a_wedge_on_the_line_of_two_notches_builds() {
    let build = |m: &mut Model| {
        let c = notched(m, BELOW_Y0);
        cut(m, c, CAP_SECOND, 3.0, 1.0)
    };
    a_wedge_on_the_line(build, [-FRAC_PI_2, 0.0], true, (8, 89));
}

/// **A window in the lateral with a side on the line** — a slot from the axis out through the
/// lateral, `y ∈ [0, 1/2]` or `[−1/2, 0]`, whose flat `y = 0` holds the line, against a wedge on
/// the side of the line away from the slot, where the rod's section is the whole disk's.
#[test]
fn a_window_on_the_line_builds() {
    let area = {
        let f = |y: f64| 0.5 * (y * (1.0 - y * y).sqrt() + y.asin());
        f(0.5) - f(0.0)
    };
    let va = 5.0 * PI - area;
    let slots: [(&[[f64; 2]], [Wedge; 2]); 2] = [
        (
            &[[0.0, 0.0], [3.0, 0.0], [3.0, 0.5], [0.0, 0.5]],
            [INTO_BELOW, OVER_BOTTOM_TANGENT],
        ),
        (
            &[[0.0, -0.5], [3.0, -0.5], [3.0, 0.0], [0.0, 0.0]],
            [INTO_ABOVE, OVER_TOP_TANGENT],
        ),
    ];
    let (mut refused, mut bodies) = (0, 0);
    for (slot, wedges) in slots {
        for w in wedges {
            for span in SPANS {
                let h = span[1] - span[0];
                let wants = four(va, w.area() * h, w.in_rays([-FRAC_PI_2, FRAC_PI_2]) * h);
                for (op, want) in OPS.into_iter().zip(wants) {
                    let mut m = Model::new();
                    let r = rod(&mut m, 0.0, 5.0);
                    let a = cut(&mut m, r, slot, 1.0, 1.0);
                    let b = prism(&mut m, &w.section(), span[0], h);
                    let label = format!("{slot:?}, {:?}..{:?}, {span:?}, {op:?}", w.a, w.b);
                    match scored(&label, &mut m, op, a, b, want) {
                        Some(n) => bodies += n,
                        None => refused += 1,
                    }
                }
            }
        }
    }
    assert_eq!((refused, bodies), (0, 64));
}

/// **A window whose sides the tangent walls hold** — the rod over `[−1/2, 3/2]` fused with the box
/// `[−1, 1] × [0, 3] × [0, 1]`, or its mirror over `y ≤ 0`: the box buries the strip of the lateral
/// on its side, and the strip's two rulings lie in the box's walls `x = ±1`, which touch the rod
/// there. The wedges stand on the line `(1, 0)` away from the box, or beside it past its wall.
#[test]
fn a_window_held_by_the_tangent_walls_builds() {
    let va = 6.0 + FRAC_PI_2 + PI;
    let (mut refused, mut bodies) = (0, 0);
    for mirrored in [false, true] {
        let (y0, y1) = if mirrored { (-3.0, 0.0) } else { (0.0, 3.0) };
        for w in [OUTSIDE, INTO_BELOW, OVER_BOTTOM_TANGENT] {
            let w = if mirrored { w.mirrored() } else { w };
            for span in [[0.25, 0.75], [-0.25, 0.5], [0.5, 1.25], [-0.25, 1.25]] {
                let h = span[1] - span[0];
                let wants = four(va, w.area() * h, w.in_rays([-FRAC_PI_2, FRAC_PI_2]) * h);
                for (op, want) in OPS.into_iter().zip(wants) {
                    let mut m = Model::new();
                    let bx = crate::fixtures::cuboid(&mut m, [-1.0, y0, 0.0], [1.0, y1, 1.0]);
                    let r = rod(&mut m, -0.5, 1.5);
                    let u = crate::boolean(&mut m, BoolKind::Fuse, bx, r).expect("the window");
                    assert_eq!(u.len(), 1, "the window is one body");
                    let b = prism(&mut m, &w.section(), span[0], h);
                    let label = format!(
                        "mirrored = {mirrored}, {:?}..{:?}, {span:?}, {op:?}",
                        w.a, w.b
                    );
                    match scored(&label, &mut m, op, u[0], b, want) {
                        Some(n) => bodies += n,
                        None => refused += 1,
                    }
                }
            }
        }
    }
    assert_eq!((refused, bodies), (6, 82));
}
