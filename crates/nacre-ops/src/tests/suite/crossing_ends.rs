//! **Laterals whose boundary sections cross where the face ends** — the scoreboard.
//!
//! The rod of [`super::oblique_cut_rod`] (radius 2 about `z` through `(0, 3)`, `z ∈ [0, 10]`) cut so
//! that its top is two slanted planes, or two slanted planes and the cap, whose sections cross on
//! the cylinder where the lateral ends. No boundary section then stands above every other, so the
//! face has no single top end; the face is the region its boundary loops close, and only their arcs —
//! each section over the angles it bounds — say where it is. Every top is a function of `x` alone:
//!
//! - **roof** — the rod cut by `prism(1/2)` and then `prism(−1/2)`: top `min(5 + x/2, 5 − x/2)`, the
//!   two arcs meeting at `(0, 1, 5)` and `(0, 5, 5)`; volume `20π − 16/3`.
//! - **valley** — the rod less the part above both planes: top `max(5 + x/2, 5 − x/2)`, the same
//!   boundary planes as the roof and another face; volume `20π + 16/3`.
//! - **steep valley** — one cut by a prism whose floor is `z = 5 + 3|x|`: top `min(10, 5 + 3|x|)`; each slanted section runs
//!   over `z ∈ [−1, 11]` and so crosses the bottom rim's plane outside its own arc; volume
//!   `40π + 32 − 40·asin(5/6) − 194√11/27`.
//! - **ridge off the axis** — planes `z = 11/2 − x/2` and `z = 9/2 + x/2`: top `5 − |x − 1|/2`;
//!   volume `58π/3 − 3√3`.
//! - **chamfered** — one cut by a prism whose floor is `z = 11 − |x|`: top `min(10, 11 − |x|)`, the
//!   ridge above the rod, so the cap's arc and the two slanted arcs bound the top; volume
//!   `128π/3 − 6√3`. Cut by two prisms in turn, the second cut already meets the first one's crossing
//!   sections, so this body is built in one.
//!
//! The tools: a far box; half-space boxes above `z = 2`, `4.5` (cutting the roof's arcs), `5`
//! (touching the roof at its two corners) and `5.5` (crossing both roof sections, missing the roof);
//! half-space boxes beyond `x = 0` (the ruling through the roof's corners) and `x = 1`; the coaxial
//! pin of radius ½ over `z ∈ [−1, 11]`; and the prism above `z = 4 + x/4`.
//!
//! **Derived, not copied.** Each common part is `∬ h dA` over the part of a disk the tool covers, with
//! `h(x)` the height the tool leaves of the body's top — piecewise linear in `x`, so on each piece
//! `∫ (p + q·x)·2√(R² − x²) dx` has a closed form ([`slab`]). The bodies' volumes are checked against
//! the closed forms above.

use super::oblique_cut_rod::{mesh_volume, prism, rod, zx_prism};
use super::*;

const PI: f64 = std::f64::consts::PI;

/// A line `z = a + b·x` in the world `zx` plane.
#[derive(Clone, Copy, Debug)]
struct Line {
    a: f64,
    b: f64,
}

impl Line {
    fn at(self, x: f64) -> f64 {
        self.a + self.b * x
    }
}

/// One of the scoreboard's bodies: how it is built, its top over the rod's disk, the lines that top
/// is made of (where it may bend), and its volume in closed form.
struct Body {
    name: &'static str,
    build: fn(&mut Model) -> Handle<Solid>,
    top: fn(f64) -> f64,
    lines: &'static [Line],
    volume: f64,
}

#[derive(Clone, Copy, Debug)]
enum Tool {
    /// `[20, 21]³`.
    Far,
    /// `[−6, 6]² × [c, 12]`.
    Above(f64),
    /// `[a, 6] × [−6, 6] × [−2, 12]`.
    Beyond(f64),
    /// The coaxial rod of radius ½ over `z ∈ [−1, 11]`.
    Pin,
    /// The prism above the plane `z = line(x)` up to `z = 12`, over `x ∈ [−3, 3]`, `y ∈ [0, 6]`.
    Over(Line),
}

impl Tool {
    fn build(self, m: &mut Model) -> Handle<Solid> {
        match self {
            Tool::Far => crate::fixtures::cuboid(m, [20.0, 20.0, 20.0], [21.0, 21.0, 21.0]),
            Tool::Above(c) => crate::fixtures::cuboid(m, [-6.0, -6.0, c], [6.0, 6.0, 12.0]),
            Tool::Beyond(a) => crate::fixtures::cuboid(m, [a, -6.0, -2.0], [6.0, 6.0, 12.0]),
            Tool::Pin => {
                crate::fixtures::cylinder_with_seam(
                    m,
                    Point3::from_array([0.0, 3.0, -1.0]),
                    Vector3::from_array([0.0, 0.0, 1.0]),
                    Vector3::from_array([1.0, 0.0, 0.0]),
                    0.5,
                    12.0,
                )
                .solid
            }
            Tool::Over(l) => zx_prism(
                m,
                &[
                    (l.at(-3.0), -3.0),
                    (l.at(3.0), 3.0),
                    (12.0, 3.0),
                    (12.0, -3.0),
                ],
            ),
        }
    }

    fn volume(self) -> f64 {
        match self {
            Tool::Far => 1.0,
            Tool::Above(c) => 144.0 * (12.0 - c),
            Tool::Beyond(a) => (6.0 - a) * 12.0 * 14.0,
            Tool::Pin => 3.0 * PI,
            // The section's area over `x ∈ [−3, 3]` is `6·(12 − a)` (the slope's part is odd).
            Tool::Over(l) => 6.0 * 6.0 * (12.0 - l.a),
        }
    }

    /// `∬ h dA` of the body's part the tool covers.
    fn common(self, body: &Body) -> f64 {
        let top = body.top;
        match self {
            Tool::Far => 0.0,
            Tool::Above(c) => slab(
                |x| (top(x) - c).max(0.0),
                2.0,
                -2.0,
                &kinks(body.lines, &[Line { a: c, b: 0.0 }]),
            ),
            Tool::Beyond(a) => slab(top, 2.0, a, &kinks(body.lines, &[])),
            Tool::Pin => slab(top, 0.5, -0.5, &kinks(body.lines, &[])),
            Tool::Over(l) => slab(
                |x| (top(x) - l.at(x)).max(0.0),
                2.0,
                -2.0,
                &kinks(body.lines, &[l]),
            ),
        }
    }
}

/// Every `x` where two of the lines meet — the only places a height made of them by `min`, `max`
/// and clamping at another of them can bend.
fn kinks(lines: &[Line], more: &[Line]) -> Vec<f64> {
    let all: Vec<Line> = lines
        .iter()
        .chain(more)
        .copied()
        .chain([Line { a: 0.0, b: 0.0 }])
        .collect();
    let mut out = Vec::new();
    for (i, p) in all.iter().enumerate() {
        for q in &all[i + 1..] {
            if p.b != q.b {
                out.push((q.a - p.a) / (p.b - q.b));
            }
        }
    }
    out
}

/// `∫_lo^R h(x)·2√(R² − x²) dx` — the volume under the height `h` over the part `x ≥ lo` of the disk
/// of radius `R`, for `h` linear between consecutive `kinks`: on each piece `h = p + q·x`, read off
/// its ends, and `∫ (p + q·x)·2√(R² − x²) dx = p·(x√(R² − x²) + R²·asin(x/R)) − q·(2/3)(R² − x²)^{3/2}`.
fn slab(h: impl Fn(f64) -> f64, r: f64, lo: f64, kinks: &[f64]) -> f64 {
    let mut cuts: Vec<f64> = kinks.iter().copied().filter(|&k| k > lo && k < r).collect();
    cuts.extend([lo, r]);
    cuts.sort_by(f64::total_cmp);
    let prim = |p: f64, q: f64, x: f64| {
        let w = (r * r - x * x).max(0.0);
        p * (x * w.sqrt() + r * r * (x / r).clamp(-1.0, 1.0).asin()) - q * 2.0 / 3.0 * w * w.sqrt()
    };
    cuts.windows(2)
        .filter(|w| w[1] > w[0])
        .map(|w| {
            let (x0, x1) = (w[0], w[1]);
            let q = (h(x1) - h(x0)) / (x1 - x0);
            let p = h(x0) - q * x0;
            prim(p, q, x1) - prim(p, q, x0)
        })
        .sum()
}

const ROOF: &[Line] = &[Line { a: 5.0, b: 0.5 }, Line { a: 5.0, b: -0.5 }];
const STEEP: &[Line] = &[
    Line { a: 5.0, b: 3.0 },
    Line { a: 5.0, b: -3.0 },
    Line { a: 10.0, b: 0.0 },
];
const OFF: &[Line] = &[Line { a: 5.5, b: -0.5 }, Line { a: 4.5, b: 0.5 }];
const CHAMFER: &[Line] = &[
    Line { a: 11.0, b: -1.0 },
    Line { a: 11.0, b: 1.0 },
    Line { a: 10.0, b: 0.0 },
];

fn bodies() -> Vec<Body> {
    let (s3, s11) = (3f64.sqrt(), 11f64.sqrt());
    vec![
        Body {
            name: "roof",
            build: |m| {
                let (r, p) = (rod(m), prism(m, 0.5));
                let a = crate::boolean(m, BoolKind::Cut, r, p).expect("the slant")[0];
                let q = prism(m, -0.5);
                crate::boolean(m, BoolKind::Cut, a, q).expect("the second slant")[0]
            },
            top: |x| (5.0 + x / 2.0).min(5.0 - x / 2.0),
            lines: ROOF,
            volume: 20.0 * PI - 16.0 / 3.0,
        },
        Body {
            name: "valley",
            build: |m| {
                let r = rod(m);
                let (p, q) = (prism(m, 0.5), prism(m, -0.5));
                let v = crate::boolean(m, BoolKind::Common, p, q).expect("the part above both")[0];
                crate::boolean(m, BoolKind::Cut, r, v).expect("the valley")[0]
            },
            top: |x| (5.0 + x / 2.0).max(5.0 - x / 2.0),
            lines: ROOF,
            volume: 20.0 * PI + 16.0 / 3.0,
        },
        Body {
            name: "steep valley",
            build: |m| {
                let r = rod(m);
                let v = zx_prism(
                    m,
                    &[
                        (14.0, -3.0),
                        (5.0, 0.0),
                        (14.0, 3.0),
                        (15.0, 3.0),
                        (15.0, -3.0),
                    ],
                );
                crate::boolean(m, BoolKind::Cut, r, v).expect("the steep valley")[0]
            },
            top: |x| (5.0 + 3.0 * x.abs()).min(10.0),
            lines: STEEP,
            volume: 40.0 * PI + 32.0 - 40.0 * (5.0f64 / 6.0).asin() - 194.0 * s11 / 27.0,
        },
        Body {
            name: "ridge off the axis",
            build: |m| {
                let r = rod(m);
                let p = zx_prism(m, &[(7.0, -3.0), (4.0, 3.0), (12.0, 3.0), (12.0, -3.0)]);
                let a = crate::boolean(m, BoolKind::Cut, r, p).expect("the first slant")[0];
                let q = zx_prism(m, &[(3.0, -3.0), (6.0, 3.0), (12.0, 3.0), (12.0, -3.0)]);
                crate::boolean(m, BoolKind::Cut, a, q).expect("the second slant")[0]
            },
            top: |x| 5.0 - (x - 1.0).abs() / 2.0,
            lines: OFF,
            volume: 58.0 * PI / 3.0 - 3.0 * s3,
        },
        Body {
            name: "chamfered",
            build: |m| {
                let r = rod(m);
                let v = zx_prism(
                    m,
                    &[
                        (8.0, -3.0),
                        (11.0, 0.0),
                        (8.0, 3.0),
                        (15.0, 3.0),
                        (15.0, -3.0),
                    ],
                );
                crate::boolean(m, BoolKind::Cut, r, v).expect("the chamfer")[0]
            },
            top: |x| (11.0 - x.abs()).min(10.0),
            lines: CHAMFER,
            volume: 128.0 * PI / 3.0 - 6.0 * s3,
        },
    ]
}

/// What a body against a tool answers in the three cells of one order: built with these bodies
/// (`[Fuse, Cut, Common]`), or refused by one name.
type Want = Result<[usize; 3], RejectReason>;

fn tools() -> [Tool; 9] {
    [
        Tool::Far,
        Tool::Above(2.0),
        Tool::Above(4.5),
        Tool::Above(5.0),
        Tool::Above(5.5),
        Tool::Beyond(0.0),
        Tool::Beyond(1.0),
        Tool::Pin,
        Tool::Over(Line { a: 4.0, b: 0.25 }),
    ]
}

/// The answer per body, per tool in [`tools`]' order, per order (the body first, then the tool):
/// built with these bodies, or refused by one name. A cell reports the lowest-numbered class's
/// refusal (`arrangement::result`'s `settle_rounds`) and the operands' order numbers the classes, so
/// the two orders can name two.
///
/// What still refuses, and why (each name's face measured):
/// - the tool's plane holds the straight edge two of the body's top faces meet on — the ridge
///   `x = 0, z = 5` (`z ≥ 5`, `x ≥ 0`; `x = 1` for the off-axis ridge) and the chamfer's
///   `x = 1, z = 10` — and the trace of one of those **planar** faces declines `NoPinOnLine`, the
///   wedge groove's refusal in `oblique_cut_rod.rs`;
/// - `z = 4.5` against the roof's slanted faces, and `z = 5.5` against the valley's, cross their
///   elliptical edges: the planar face's scan cannot name a point on an edge riding the cylinder
///   (`CurvedRingWall`);
/// - the steep valley: a slant's trace runs along its own arc from the valley's bottom corner,
///   where it goes on inside the face, to the cap's corner, where it goes on outside — a run whose
///   flanks differ, which `cycle_on_class` declines (measured at that arm; `ObliqueCylinderCut`);
///   with the tool first,
///   `z ≥ 5`, `z ≥ 5.5` and `x ≥ 0` name their own refusals before it;
/// - the gate: the off-axis ridge's `z = 4.5`, `z = 5.5` and the prism's floor touch its
///   slanted sections, and the steep valley's pin.
fn wants() -> [[[Want; 2]; 9]; 5] {
    let decline = |kind| Err(RejectReason::TraceDeclined { kind, face: None });
    let ring_wall: Want = decline(crate::DeclineKind::CurvedRingWall);
    let no_pin: Want = decline(crate::DeclineKind::NoPinOnLine);
    // The gate's refusal and an oblique trace's decline share this name.
    let oblique: Want = Err(RejectReason::ObliqueCylinderCut);
    let both = |w: Want| [w, w];
    // Far apart: two bodies fused, the body itself cut, nothing in common.
    let apart = both(Ok([2, 1, 0]));
    let one = both(Ok([1, 1, 1]));
    // The pin pokes out above the top and below the floor: the pin less the body is two.
    let pin: [Want; 2] = [Ok([1, 1, 1]), Ok([1, 2, 1])];
    [
        // roof
        [
            apart,
            one,
            both(ring_wall),
            both(no_pin),
            apart,
            both(no_pin),
            one,
            pin,
            one,
        ],
        // valley
        [
            apart,
            one,
            one,
            both(no_pin),
            both(ring_wall),
            both(no_pin),
            one,
            pin,
            one,
        ],
        // steep valley
        [
            both(oblique),
            both(oblique),
            both(oblique),
            [oblique, no_pin],
            [oblique, ring_wall],
            [oblique, no_pin],
            both(oblique),
            both(oblique),
            both(oblique),
        ],
        // ridge off the axis
        [
            apart,
            one,
            both(oblique),
            both(no_pin),
            both(oblique),
            one,
            both(no_pin),
            pin,
            both(oblique),
        ],
        // chamfered
        [apart, one, one, one, one, one, both(no_pin), pin, one],
    ]
}

/// **Every body against every tool, every kind, both orders** — each cell one boolean on a fresh
/// model: refused ones by their name, built ones held to their body count, `validate`, the mesh, and
/// the boolean identities with the derived common part. Mismatches are reported together.
#[test]
fn laterals_whose_boundary_sections_cross_where_they_end() {
    let same_name = |a: RejectReason, b: RejectReason| match (a, b) {
        (
            RejectReason::TraceDeclined { kind: x, .. },
            RejectReason::TraceDeclined { kind: y, .. },
        ) => x == y,
        _ => a == b,
    };
    let mut bad: Vec<String> = Vec::new();
    for (body, wants) in bodies().iter().zip(wants()) {
        let mut m = Model::new();
        let b = (body.build)(&mut m);
        let got = nacre_props::mass_props(&m, b).expect("props").volume;
        let derived = slab(body.top, 2.0, -2.0, &kinks(body.lines, &[]));
        assert!(
            (derived - body.volume).abs() <= 1e-12 * body.volume,
            "{}: the integral {derived} vs the closed form {}",
            body.name,
            body.volume
        );
        assert!(
            (got - body.volume).abs() <= 1e-12 * body.volume,
            "{}: {got} vs {}",
            body.name,
            body.volume
        );
        for (tool, want) in tools().into_iter().zip(wants) {
            let c = tool.common(body);
            for swap in [false, true] {
                for (k, kind) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
                    .into_iter()
                    .enumerate()
                {
                    let at = format!("{} × {tool:?} — {kind:?} swap={swap}", body.name);
                    let mut m = Model::new();
                    let (b, t) = ((body.build)(&mut m), tool.build(&mut m));
                    let (x, y) = if swap { (t, b) } else { (b, t) };
                    let out = crate::boolean(&mut m, kind, x, y);
                    match (want[usize::from(swap)], out) {
                        (Err(name), Err(BoolError::Rejected { reason, .. }))
                            if same_name(reason, name) => {}
                        (Ok(counts), Ok(out)) => {
                            let n = counts[k];
                            let got = built_volume(&m, &out, n, &at);
                            let derived = match (kind, swap) {
                                (BoolKind::Fuse, _) => body.volume + tool.volume() - c,
                                (BoolKind::Cut, false) => body.volume - c,
                                (BoolKind::Cut, true) => tool.volume() - c,
                                (BoolKind::Common, _) => c,
                            };
                            if (got - derived).abs() > 1e-9 * derived.max(1.0) {
                                bad.push(format!("{at}: {got} vs {derived}"));
                            }
                            let mesh = mesh_volume(&m, &out);
                            if (mesh - derived).abs() > 2e-3 * derived.max(1.0) {
                                bad.push(format!("{at}: mesh {mesh} vs {derived}"));
                            }
                        }
                        (want, got) => bad.push(format!("{at}: want {want:?}, got {got:?}")),
                    }
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} cells:\n{}", bad.len(), bad.join("\n"));
}
