//! **A lateral's boundary running along a section and leaving it** — the unit rod with a hole that
//! has a stretch of its boundary on a plane across the axis, read against that plane
//! (`cycle_on_class`, `combinatorics::past_run_ends`). Each body is cut from the rod by rectangular
//! prisms:
//!
//! - **stepped** — the half `y < 0` over `z ∈ [1, 2]`, and `y < −1/2` over `[2, 3]`. On `z = 2` the
//!   boundary runs along the arc between the two walls and goes on to the far side of the plane at
//!   each end of it;
//! - **panel, half notch** — a panel `x > −1/2` of the rod with the same half notch;
//! - **panel, notch off the axis** — that panel with `y < −1/2` taken over `[1, 2]`;
//! - **cape** — `y < 0` taken over `[2, 3]` and over `[1, 2]` except `|x| < 2/5`: material standing
//!   into the hole, whose top arc lies on `z = 2` with the hole past both of its ends.
//!
//! Each is also the bore it leaves in a block (`[−2, 2]² × [0, 5]` less the body), whose lateral
//! faces the other way.
//!
//! **Derived, not read back.** A body's section on each stretch of height is a set of disjoint
//! rectangles met with the unit disk, so every common part is a sum of disk–rectangle areas in
//! closed form ([`disk_in_rect`]). The body counts come from a slab model of the same sections
//! outside the kernel. OCCT agrees with the volume and area of every cell of the seven boxes that
//! reach the body, as built, to the six digits it prints.

use super::one_sided_notch::{four, prism, scored};
use super::scoreboard::OPS;
use super::*;
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};

/// `[x0, y0, x1, y1]`.
type Rect = [f64; 4];

const WHOLE: Rect = [-2.0, -2.0, 2.0, 2.0];
const BELOW: Rect = [-3.0, -3.0, 3.0, 0.0];
const BELOW_HALF: Rect = [-3.0, -3.0, 3.0, -0.5];
const PAST: Rect = [-3.0, -3.0, -0.5, 3.0];
const LEFT: Rect = [-3.0, -3.0, -0.4, 0.0];
const RIGHT: Rect = [0.4, -3.0, 3.0, 0.0];

/// The stretches of height a body's section is constant on.
const SLABS: [[f64; 2]; 4] = [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 5.0]];

/// A body: the prisms cut from the rod (`rect`, from `z`, for `h`), and the rectangles its section
/// keeps of the disk on each of [`SLABS`].
struct Body {
    cuts: &'static [(Rect, f64, f64)],
    kept: [&'static [Rect]; 4],
}

const STEPPED: Body = Body {
    cuts: &[(BELOW, 1.0, 1.0), (BELOW_HALF, 2.0, 1.0)],
    kept: [
        &[WHOLE],
        &[[-2.0, 0.0, 2.0, 2.0]],
        &[[-2.0, -0.5, 2.0, 2.0]],
        &[WHOLE],
    ],
};
const PANEL: Rect = [-0.5, -2.0, 2.0, 2.0];
const PANEL_HALF: Body = Body {
    cuts: &[(PAST, -1.0, 7.0), (BELOW, 1.0, 1.0)],
    kept: [&[PANEL], &[[-0.5, 0.0, 2.0, 2.0]], &[PANEL], &[PANEL]],
};
const PANEL_OFF: Body = Body {
    cuts: &[(PAST, -1.0, 7.0), (BELOW_HALF, 1.0, 1.0)],
    kept: [&[PANEL], &[[-0.5, -0.5, 2.0, 2.0]], &[PANEL], &[PANEL]],
};
const CAPE: Body = Body {
    cuts: &[(BELOW, 2.0, 1.0), (LEFT, 1.0, 1.0), (RIGHT, 1.0, 1.0)],
    kept: [
        &[WHOLE],
        &[[-2.0, 0.0, 2.0, 2.0], [-0.4, -2.0, 0.4, 0.0]],
        &[[-2.0, 0.0, 2.0, 2.0]],
        &[WHOLE],
    ],
};

/// The area of the unit disk inside `r`: `∫ (top − bottom) dx` over the pieces of `x` between which
/// each of the two is a side of the rectangle or the circle.
fn disk_in_rect([x0, y0, x1, y1]: Rect) -> f64 {
    let s = |x: f64| (1.0 - x * x).max(0.0).sqrt();
    let g = |x: f64| 0.5 * (x * s(x) + x.clamp(-1.0, 1.0).asin());
    let (x0, x1) = (x0.max(-1.0), x1.min(1.0));
    if x0 >= x1 || y0 >= y1 {
        return 0.0;
    }
    let mut cuts = vec![x0, x1];
    for y in [y0, y1] {
        if y.abs() < 1.0 {
            cuts.extend([-s(y), s(y)].into_iter().filter(|&x| x > x0 && x < x1));
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.windows(2)
        .map(|w| {
            let (a, b) = (w[0], w[1]);
            let mid = 0.5 * (a + b);
            let circle = g(b) - g(a);
            if y1.min(s(mid)) <= y0.max(-s(mid)) {
                return 0.0;
            }
            let top = if y1 < s(mid) { y1 * (b - a) } else { circle };
            let bottom = if y0 > -s(mid) { y0 * (b - a) } else { -circle };
            top - bottom
        })
        .sum()
}

fn clip(a: Rect, b: Rect) -> Rect {
    [
        a[0].max(b[0]),
        a[1].max(b[1]),
        a[2].min(b[2]),
        a[3].min(b[3]),
    ]
}

fn area(r: Rect) -> f64 {
    (r[2] - r[0]).max(0.0) * (r[3] - r[1]).max(0.0)
}

impl Body {
    /// The body's volume inside `rect × z`, or the volume of the bore it leaves in the block.
    fn within(&self, bore: bool, rect: Rect, z: [f64; 2]) -> f64 {
        SLABS
            .iter()
            .zip(self.kept)
            .map(|(slab, kept)| {
                let dz = (z[1].min(slab[1]) - z[0].max(slab[0])).max(0.0);
                let rod: f64 = kept.iter().map(|&k| disk_in_rect(clip(k, rect))).sum();
                dz * if bore {
                    area(clip(WHOLE, rect)) - rod
                } else {
                    rod
                }
            })
            .sum()
    }

    /// The body as built, turned, or reflected — `place` `0`, `1..=3` (a half turn about `x`, a half
    /// turn about `y`, a quarter about `z`), `4..=6` (in `x = 0`, `y = 0`, `z = 5/2`). `Mirror` does
    /// not take a curved solid, so a reflected body is built reflected. `down` states the rod's axis
    /// along `−z`.
    fn build(&self, m: &mut Model, bore: bool, down: bool, place: usize) -> Handle<Solid> {
        let (base, axis) = if down {
            ([0.0, 0.0, 5.0], [0.0, 0.0, -1.0])
        } else {
            ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        };
        let mut x = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array(base),
            nacre_math::Vector3::from_array(axis),
            nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
            1.0,
            5.0,
        )
        .solid;
        for &([x0, y0, x1, y1], z, h) in self.cuts {
            let [x0, y0, x1, y1] = match place {
                4 => [-x1 + 0.0, y0, -x0 + 0.0, y1],
                5 => [x0, -y1 + 0.0, x1, -y0 + 0.0],
                _ => [x0, y0, x1, y1],
            };
            let z = if place == 6 { 5.0 - (z + h) } else { z };
            let t = prism(m, &[[x0, y0], [x1, y0], [x1, y1], [x0, y1]], z, h);
            let out = crate::boolean(m, BoolKind::Cut, x, t).expect("the cut");
            assert_eq!(out.len(), 1, "the rod stays one body");
            x = out[0];
        }
        if bore {
            let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 0.0], [2.0, 2.0, 5.0]);
            let out = crate::boolean(m, BoolKind::Cut, block, x).expect("the bore");
            assert_eq!(out.len(), 1, "the block stays one body");
            x = out[0];
        }
        let turn = |axis, deg| Operation::Transform {
            solid: x,
            isometry: Isometry::rotation(Rotation {
                axis,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
            }),
        };
        let op = match place {
            1 => turn(Axis::X, 180),
            2 => turn(Axis::Y, 180),
            3 => turn(Axis::Z, 90),
            _ => return x,
        };
        match apply(m, &op).expect("the placement") {
            OpOutput::Transform { solid } => solid,
            other => panic!("{other:?}"),
        }
    }
}

/// `rect × z`: a box far off, and boxes through the hole, its walls and the plane of the run. No wall
/// of these is tangent to the rod.
const TOOLS: [(Rect, [f64; 2]); 8] = [
    ([5.0, 5.0, 6.0, 6.0], [1.5, 2.5]),
    ([-2.0, -0.25, 2.0, 0.25], [1.5, 2.5]),
    ([-0.75, -2.0, 0.75, -0.75], [0.5, 3.5]),
    ([0.5, -0.25, 1.5, 0.25], [1.25, 1.75]),
    ([-0.5, -0.5, 0.5, 0.0], [1.0, 2.0]),
    (WHOLE, [1.25, 1.75]),
    (WHOLE, [2.0, 3.0]),
    ([-0.75, -0.25, 0.75, 0.25], [-1.0, 6.0]),
];

/// **Every tool, every boolean.** Each builds, to its closed-form volume, `validate` and the mesh.
/// `bodies` is the slab model's count for each cell: four digits a tool, in the order of [`OPS`].
fn against_the_tools(body: &Body, bore: bool, bodies: &str) {
    let va = body.within(bore, [-9.0, -9.0, 9.0, 9.0], [-9.0, 9.0]);
    let mut got = Vec::with_capacity(TOOLS.len());
    for (k, (rect, z)) in TOOLS.into_iter().enumerate() {
        let common = body.within(bore, rect, z);
        let wants = four(va, area(rect) * (z[1] - z[0]), common);
        let mut cells = String::new();
        for (op, want) in OPS.into_iter().zip(wants) {
            let mut m = Model::new();
            let a = body.build(&mut m, bore, false, 0);
            let lo = [rect[0], rect[1], z[0]];
            let b = crate::fixtures::cuboid(&mut m, lo, [rect[2], rect[3], z[1]]);
            let label = format!("bore = {bore}, tool {k}, {op:?}");
            let n = scored(&label, &mut m, op, a, b, want)
                .unwrap_or_else(|| panic!("{label}: refused"));
            cells.push_str(&n.to_string());
        }
        got.push(cells);
    }
    assert_eq!(got.join(" "), bodies);
}

#[test]
fn a_stepped_hole_takes_a_boolean() {
    against_the_tools(&STEPPED, false, "2110 1111 1112 1111 1110 1211 1211 1131");
}

#[test]
fn a_stepped_bore_takes_a_boolean() {
    against_the_tools(&STEPPED, true, "2110 1111 1121 1111 1101 1211 1211 1111");
}

#[test]
fn a_half_notch_at_a_panels_edge_takes_a_boolean() {
    against_the_tools(
        &PANEL_HALF,
        false,
        "2110 1111 1112 1111 1110 1211 1211 1111",
    );
}

#[test]
fn the_bore_of_a_half_notch_at_a_panels_edge_takes_a_boolean() {
    against_the_tools(&PANEL_HALF, true, "2110 1111 1121 1111 1101 1211 1211 1111");
}

#[test]
fn a_notch_off_the_axis_at_a_panels_edge_takes_a_boolean() {
    against_the_tools(&PANEL_OFF, false, "2110 1121 1112 1111 1101 1211 1211 1111");
}

#[test]
fn the_bore_of_a_notch_off_the_axis_at_a_panels_edge_takes_a_boolean() {
    against_the_tools(&PANEL_OFF, true, "2110 1112 1121 1111 1110 1211 1211 1111");
}

#[test]
fn a_cape_standing_into_a_hole_takes_a_boolean() {
    against_the_tools(&CAPE, false, "2110 1111 1112 1111 1121 1211 1211 1131");
}

#[test]
fn the_bore_of_a_cape_takes_a_boolean() {
    against_the_tools(&CAPE, true, "2110 1111 1121 1111 1112 1211 1211 1111");
}

/// The rod less a notch whose wall meets a slanted roof `z = 2 + x`, so an arc on the roof's section
/// runs into a corner and the next arc leaves it and crosses that section again; `mirrored` reflects
/// it in `x = 0`, which walks that arc from the other end of the run. `bore` is what it leaves in the
/// block.
fn slanted_run(m: &mut Model, bore: bool, mirrored: bool) -> Handle<Solid> {
    let roof = [(0.5, -0.5), (0.5, 0.0), (2.0, 0.0), (1.5, -0.5)];
    let wedge = [
        [0.0, 0.0],
        [0.0, 3.0],
        [3.0, 3.0],
        [3.0, -3.0],
        [-1.5, -3.0],
    ];
    let (roof, wedge): (Vec<_>, Vec<_>) = if mirrored {
        (
            roof.iter().rev().map(|&(z, x)| (z, -x + 0.0)).collect(),
            wedge.iter().rev().map(|&[x, y]| [-x + 0.0, y]).collect(),
        )
    } else {
        (roof.to_vec(), wedge.to_vec())
    };
    let rod = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([0.0; 3]),
        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
        nacre_math::Vector3::from_array([1.0, 0.0, 0.0]),
        1.0,
        5.0,
    )
    .solid;
    let t = super::oblique_cut_rod::zx_prism(m, &roof);
    let w = prism(m, &wedge, 0.5, 1.5);
    let cutter = crate::boolean(m, BoolKind::Fuse, t, w).expect("the cutter");
    let out = crate::boolean(m, BoolKind::Cut, rod, cutter[0]).expect("the notch");
    assert_eq!(out.len(), 1, "the rod stays one body");
    if !bore {
        return out[0];
    }
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 0.0], [2.0, 2.0, 5.0]);
    let out = crate::boolean(m, BoolKind::Cut, block, out[0]).expect("the bore");
    assert_eq!(out.len(), 1, "the block stays one body");
    out[0]
}

/// What [`slanted_run`] holds in `rect × z`, unmirrored: the rod (or the bore in the block) less the
/// wedge `x > 0 ∨ y < 2x` over `z ∈ [1/2, 2]` and the strip `y > 0, −1/2 < x < 0` under the roof
/// `z = 2 + x`. A column's height is constant in `y` between the breakpoints of its `x`, so the inner
/// integral is an exact sum; the outer one is Gauss–Legendre between the `x` where those breakpoints
/// change order, in `x = sin φ` on the disk, where the rim's square root is smooth.
fn slanted_within(bore: bool, [x0, y0, x1, y1]: Rect, [lo, hi]: [f64; 2]) -> f64 {
    const NODES: [(f64, f64); 5] = [
        (0.0, 0.568_888_888_888_888_9),
        (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_1),
        (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    ];
    let over = |a: f64, b: f64| (b.min(hi) - a.max(lo)).max(0.0);
    let rim = |x: f64| (1.0 - x * x).max(0.0).sqrt();
    let (y0, y1) = if bore {
        (y0.max(-2.0), y1.min(2.0))
    } else {
        (y0, y1)
    };
    let height = |x: f64, y: f64| {
        let whole = over(0.0, 5.0);
        if x * x + y * y >= 1.0 {
            return if bore { whole } else { 0.0 };
        }
        let gone = if x > 0.0 || y < 2.0 * x {
            over(0.5, 2.0)
        } else if y > 0.0 && x > -0.5 {
            over(0.5, 2.0 + x)
        } else {
            0.0
        };
        if bore { gone } else { whole - gone }
    };
    let column = |x: f64| {
        let mut ys = vec![y0, y1];
        ys.extend(
            [-rim(x), rim(x), 0.0, 2.0 * x]
                .into_iter()
                .filter(|&y| y > y0 && y < y1),
        );
        ys.sort_by(f64::total_cmp);
        ys.windows(2)
            .map(|w| (w[1] - w[0]) * height(x, 0.5 * (w[0] + w[1])))
            .sum::<f64>()
    };
    let reach = if bore { 2.0 } else { 1.0 };
    let (x0, x1) = (x0.max(-reach), x1.min(reach));
    if x0 >= x1 || y0 >= y1 {
        return 0.0;
    }
    let mut xs = vec![
        x0,
        x1,
        -1.0,
        1.0,
        -0.5,
        0.0,
        5f64.sqrt() / 5.0,
        -(5f64.sqrt()) / 5.0,
    ];
    xs.extend([lo - 2.0, hi - 2.0, 0.5 * y0, 0.5 * y1]);
    for y in [y0, y1] {
        if y.abs() < 1.0 {
            xs.extend([-rim(y), rim(y)]);
        }
    }
    xs.retain(|&x| x >= x0 && x <= x1);
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    let panel = |a: f64, b: f64, f: &dyn Fn(f64) -> f64| {
        const PANELS: usize = 12;
        let h = (b - a) / PANELS as f64;
        (0..PANELS)
            .map(|k| {
                let mid = a + (k as f64 + 0.5) * h;
                NODES
                    .iter()
                    .map(|&(t, w)| 0.5 * h * w * f(mid + 0.5 * h * t))
                    .sum::<f64>()
            })
            .sum::<f64>()
    };
    xs.windows(2)
        .map(|w| {
            let (a, b) = (w[0], w[1]);
            if a >= -1.0 && b <= 1.0 {
                panel(a.asin(), b.asin(), &|p: f64| column(p.sin()) * p.cos())
            } else {
                panel(a, b, &column)
            }
        })
        .sum()
}

/// **A run read from either end.** The far box's volumes are closed-form:
/// `15π/4 + (3/4)·atan 2 − 5√3/16 + 1/3`, the rod less the wedge and the strip under the roof; the
/// boxes through the body are held to [`slanted_within`], which the closed form checks, and the
/// reflected body to the same numbers.
#[test]
fn a_run_on_a_slanted_section_reads_alike_from_either_end() {
    let rod = 15.0 * std::f64::consts::FRAC_PI_4 + 0.75 * 2.0f64.atan()
        - 5.0 * 3.0f64.sqrt() / 16.0
        + 1.0 / 3.0;
    for bore in [false, true] {
        let va = if bore { 80.0 - rod } else { rod };
        let everywhere = slanted_within(bore, [-9.0, -9.0, 9.0, 9.0], [-9.0, 9.0]);
        assert!(
            (everywhere - va).abs() < 1e-9,
            "the reference: {everywhere} against {va}"
        );
        for k in [0, 1, 2, 7] {
            let (rect, z) = TOOLS[k];
            let vt = area(rect) * (z[1] - z[0]);
            let box_of = |m: &mut Model, mirrored: bool| {
                let [x0, y0, x1, y1] = rect;
                let (x0, x1) = if mirrored {
                    (-x1 + 0.0, -x0 + 0.0)
                } else {
                    (x0, x1)
                };
                crate::fixtures::cuboid(m, [x0, y0, z[0]], [x1, y1, z[1]])
            };
            let common = slanted_within(bore, rect, z);
            for mirrored in [false, true] {
                for (op, want) in OPS.into_iter().zip(four(va, vt, common)) {
                    let mut m = Model::new();
                    let a = slanted_run(&mut m, bore, mirrored);
                    let b = box_of(&mut m, mirrored);
                    let label = format!("bore = {bore}, mirrored = {mirrored}, tool {k}, {op:?}");
                    scored(&label, &mut m, op, a, b, want)
                        .unwrap_or_else(|| panic!("{label}: refused"));
                }
            }
        }
    }
}

/// **Every frame the run can be read in.** Rod and bore, the axis stated either way, each placement,
/// less a far box: the end a run's stretch of absence starts from is a product of the lateral's
/// facing, the arc's sense, the class's stored normal and the walk's side, and a wrong one is refused
/// here whatever the tool.
fn in_every_frame(body: &Body) {
    let va = body.within(false, [-9.0, -9.0, 9.0, 9.0], [-9.0, 9.0]);
    let vb = body.within(true, [-9.0, -9.0, 9.0, 9.0], [-9.0, 9.0]);
    for bore in [false, true] {
        for down in [false, true] {
            for place in 0..7 {
                let mut m = Model::new();
                let a = body.build(&mut m, bore, down, place);
                let far = crate::fixtures::cuboid(&mut m, [5.0; 3], [6.0; 3]);
                let label = format!("bore = {bore}, down = {down}, place {place}");
                let out = crate::boolean(&mut m, BoolKind::Cut, a, far)
                    .unwrap_or_else(|e| panic!("{label}: {e:?}"));
                assert_eq!(out.len(), 1, "{label}");
                assert!(nacre_validate::validate(&m).is_empty(), "{label}");
                let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                let want = if bore { vb } else { va };
                assert!((v - want).abs() < 1e-9, "{label}: {v}, want {want}");
            }
        }
    }
}

#[test]
fn a_stepped_hole_reads_alike_in_every_frame() {
    in_every_frame(&STEPPED);
}

#[test]
fn a_half_notch_at_a_panels_edge_reads_alike_in_every_frame() {
    in_every_frame(&PANEL_HALF);
}

#[test]
fn a_notch_off_the_axis_at_a_panels_edge_reads_alike_in_every_frame() {
    in_every_frame(&PANEL_OFF);
}

#[test]
fn a_cape_reads_alike_in_every_frame() {
    in_every_frame(&CAPE);
}
