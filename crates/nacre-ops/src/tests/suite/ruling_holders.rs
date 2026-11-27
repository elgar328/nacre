//! **A lateral's boundary along one ruling whose holder changes** — the unit rod (seam on `+x̂`)
//! cut by prisms stacked in height, so that one ruling of the hole's edge is held by one plane below
//! a height and another above it. Along that ruling the cycle runs on through the node where the
//! planes change: a class across the axis there crosses the boundary at a node between two ruling
//! edges (`cycle_on_class`, `ruling_hole`), and each wall class holding a stretch of the run reads
//! its own stretch and the other plane's as `quiet` (`cycle_on_wall`).
//!
//! The bodies (every cut a prism over a convex polygon, from a height, for a height):
//!
//! - **half, then slope** — `y < 0` over `z ∈ [1, 2]`, then `x − 2y > 1` over `[2, 3]`: the ruling
//!   `(1, 0)`, on the seam, is held by `y = 0` and then by `x − 2y = 1`; its mirror in `x = 0` holds
//!   the change on `(−1, 0)`, off the seam;
//! - **slope, then half** — the same two the other way up;
//! - **kinked wedges** — two wedges over `[1, 2]` and `[2, 3]` whose kinked chords both run from
//!   `(1, 0)` to `(0, −1)`, so both rulings change hands;
//! - **three times** — half, slope, half over `[1, 2]`, `[2, 3]`, `[3, 4]`: one run in three
//!   stretches;
//! - **to the top** — half over `[1, 2]`, slope from `2` past the rod's top: the run ends on the top
//!   rim, past which the face holds neither side;
//! - **rim to rim** — slope from below the rod's foot to `2`, half from `2` past its top: the run
//!   ends on both rims;
//! - **tangent box** — the rod joined to the box `[0, 1] × [0, 2] × [1, 2]`, whose wall `x = 1`
//!   touches the rod along `(1, 0)`, then `y > 0` cut over `[2, 3]`: the tangent class holds the
//!   run's lower stretch;
//! - **under a roof** — `y > 0` cut from `1` up to the oblique roof `z = 2 + x/4`, and `x + 2y > 1`
//!   from the roof up to `7/2`: the planes change at a node on an oblique class.
//!
//! Each is also the bore it leaves in the block `[−2, 2]² × [0, 5]`.
//!
//! **What is derived, and what is read.** A prism body's section on each stretch of height is the
//! unit disk less one convex polygon (or joined to one), so every common part with a box is a sum of
//! areas of the disk inside convex polygons, in closed form ([`disk_in`]). The body counts come from
//! a slab model of the same sections outside the kernel — but where two bodies meet only along a
//! line (a box's wall tangent to the rod, or its corner on a wall), which the model's polygon of the
//! circle merges and the kernel answers as two — and so do the cells refused as self-contact (`s`,
//! there or a box's corner standing on a wall). The roof's sections move with height: its cells are
//! held instead to the box (`B − A` with `A ∩ B`), to the three parts of `A ∪ B` and to the same
//! cells turned, and its body counts are the kernel's; an `o` or `g` cell is the refusal its cuts
//! make without the change. A `u` cell is the name the kernel gives, which no model outside it
//! supplies: the box's wall holds the ruling with the two wedge walls, a pencil the class table has
//! not learned, so the class's order around a node of that line sees two edges leave along one ray
//! (`UnorderedEdges` — todo «한 선을 품는 두 벽»).

use super::one_sided_notch::{four, prism};
use super::scoreboard::{BoolResult, OPS, run};
use super::*;

type Pts = &'static [[f64; 2]];
/// Polygons, each a list of corners.
type Polygons = Vec<Vec<[f64; 2]>>;

const BELOW: Pts = &[[-3.0, -3.0], [3.0, -3.0], [3.0, 0.0], [-3.0, 0.0]];
/// `x − 2y > 1`: its wall meets the circle at `(1, 0)` and `(−3/5, −4/5)`.
const SLOPE: Pts = &[[-3.0, -3.0], [3.0, -3.0], [3.0, 1.0], [-3.0, -2.0]];
const KINK_P: Pts = &[[1.0, 0.0], [0.25, -0.25], [0.0, -1.0], [3.0, -3.0]];
const KINK_Q: Pts = &[[1.0, 0.0], [0.1, -0.5], [0.0, -1.0], [3.0, -3.0]];
const ABOVE: Pts = &[[-3.0, 0.0], [3.0, 0.0], [3.0, 3.0], [-3.0, 3.0]];
/// `x + 2y > 1`, which within the disk lies in `y > 0`.
const SLOPE_UP: Pts = &[[-3.0, 2.0], [3.0, -1.0], [3.0, 3.0], [-3.0, 3.0]];
/// The block a bore is cut in.
const BLOCK: Pts = &[[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]];
/// The tangent box's footprint.
const TANGENT_BOX: Pts = &[[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 2.0]];

#[derive(Clone, Copy)]
enum Body {
    /// The rod less each prism `(polygon, from, height)`.
    Cuts(&'static [(Pts, f64, f64)]),
    Tangent,
    Roof,
}

const HALF_SLOPE: Body = Body::Cuts(&[(BELOW, 1.0, 1.0), (SLOPE, 2.0, 1.0)]);
const SLOPE_HALF: Body = Body::Cuts(&[(SLOPE, 1.0, 1.0), (BELOW, 2.0, 1.0)]);
const KINKS: Body = Body::Cuts(&[(KINK_P, 1.0, 1.0), (KINK_Q, 2.0, 1.0)]);
const THRICE: Body = Body::Cuts(&[(BELOW, 1.0, 1.0), (SLOPE, 2.0, 1.0), (BELOW, 3.0, 1.0)]);
const TO_TOP: Body = Body::Cuts(&[(BELOW, 1.0, 1.0), (SLOPE, 2.0, 3.5)]);
const RIM_TO_RIM: Body = Body::Cuts(&[(SLOPE, -0.5, 2.5), (BELOW, 2.0, 3.5)]);

/// `[x0, y0, z0, x1, y1, z1]`: a box far off, and boxes through the hole across the height where the
/// walls change, around the changing ruling, across the whole rod there, with a wall `x = 0` holding
/// the ruling `(0, −1)`, outside the rod with a wall `y = 0` holding `(1, 0)`, through the axis,
/// through the walls' far side, and two across the change with a wall tangent to the rod on the
/// changing ruling `(1, 0)` and on its mirror image `(−1, 0)`.
const TOOLS: [[f64; 6]; 10] = [
    [5.0, 5.0, 1.5, 6.0, 6.0, 2.5],
    [-0.75, -1.5, 1.5, 0.75, -0.4, 2.5],
    [0.5, -0.25, 1.25, 1.5, 0.25, 2.75],
    [-2.0, -2.0, 1.75, 2.0, 2.0, 2.25],
    [0.0, -0.75, 2.5, 1.25, 1.5, 3.5],
    [1.25, -1.5, 2.5, 1.5, 0.0, 3.5],
    [-0.4, -0.4, 0.5, 0.2, 0.2, 3.5],
    [-1.5, -1.5, 0.5, 1.5, -0.75, 3.5],
    [1.0, -0.5, 1.5, 1.5, 0.5, 2.5],
    [-1.5, -0.5, 1.5, -1.0, 0.5, 2.5],
];

/// `pts` reflected in `x = 0`, kept counter-clockwise.
fn mirrored(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    pts.iter().rev().map(|&[x, y]| [-x + 0.0, y]).collect()
}

fn one(got: BoolResult, what: &str) -> Handle<Solid> {
    let out = got.unwrap_or_else(|e| panic!("{what}: {e:?}"));
    assert_eq!(out.len(), 1, "{what}: one body");
    out[0]
}

impl Body {
    /// The body, reflected in `x = 0` if `mirror`, or the bore it leaves in the block; `down` states
    /// the rod's axis along `−z`.
    fn build(self, m: &mut Model, bore: bool, mirror: bool, down: bool) -> Handle<Solid> {
        let (base, axis) = if down {
            ([0.0, 0.0, 5.0], [0.0, 0.0, -1.0])
        } else {
            ([0.0; 3], [0.0, 0.0, 1.0])
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
        let place = |pts: Pts| if mirror { mirrored(pts) } else { pts.to_vec() };
        match self {
            Body::Cuts(cuts) => {
                for &(pts, z, h) in cuts {
                    let t = prism(m, &place(pts), z, h);
                    x = one(crate::boolean(m, BoolKind::Cut, x, t), "a cut");
                }
            }
            Body::Tangent => {
                let t = prism(m, &place(TANGENT_BOX), 1.0, 1.0);
                x = one(crate::boolean(m, BoolKind::Fuse, x, t), "the tangent box");
                let t = prism(m, &place(ABOVE), 2.0, 1.0);
                x = one(crate::boolean(m, BoolKind::Cut, x, t), "the cut");
            }
            Body::Roof => {
                let below: [(f64, f64); 4] = [(1.0, -3.0), (1.0, 3.0), (2.75, 3.0), (1.25, -3.0)];
                let above: [(f64, f64); 4] = [(1.25, -3.0), (2.75, 3.0), (3.5, 3.0), (3.5, -3.0)];
                for (roof, pts) in [(below, ABOVE), (above, SLOPE_UP)] {
                    // `(z, x)` in the frame `zx_prism` extrudes along `+y`.
                    let roof: Vec<(f64, f64)> = if mirror {
                        roof.iter().rev().map(|&(z, x)| (z, -x + 0.0)).collect()
                    } else {
                        roof.to_vec()
                    };
                    let r = super::oblique_cut_rod::zx_prism(m, &roof);
                    let p = prism(m, &place(pts), 0.5, 3.5);
                    let t = one(
                        crate::boolean(m, BoolKind::Common, r, p),
                        "the roofed prism",
                    );
                    x = one(crate::boolean(m, BoolKind::Cut, x, t), "a cut");
                }
            }
        }
        if bore {
            let block = prism(m, BLOCK, 0.0, 5.0);
            x = one(crate::boolean(m, BoolKind::Cut, block, x), "the bore");
        }
        x
    }

    /// The disk's section over the height `z`: the polygons it is cut by there (`minus`) and joined
    /// to (`plus`), reflected if `mirror`.
    fn section(self, mirror: bool, z: f64) -> (Polygons, Polygons) {
        let place = |pts: Pts| if mirror { mirrored(pts) } else { pts.to_vec() };
        let within = |lo: f64, h: f64| lo < z && z < lo + h;
        match self {
            Body::Cuts(cuts) => (
                cuts.iter()
                    .filter(|c| within(c.1, c.2))
                    .map(|c| place(c.0))
                    .collect(),
                Vec::new(),
            ),
            Body::Tangent => (
                if within(2.0, 1.0) {
                    vec![place(ABOVE)]
                } else {
                    Vec::new()
                },
                if within(1.0, 1.0) {
                    vec![place(TANGENT_BOX)]
                } else {
                    Vec::new()
                },
            ),
            Body::Roof => unreachable!("the roof's sections move with height"),
        }
    }

    /// The body's volume inside the box `[x0, y0, z0, x1, y1, z1]`, or the bore's, in closed form:
    /// on each stretch of constant section, the disk inside the box less its parts inside the
    /// polygons it is cut by, plus the joined polygons' parts inside the box and outside the disk.
    /// The polygons cut from one stretch do not overlap.
    fn within(self, bore: bool, mirror: bool, b: [f64; 6]) -> f64 {
        let rect = vec![[b[0], b[1]], [b[3], b[1]], [b[3], b[4]], [b[0], b[4]]];
        let mut zs = [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            b[2].clamp(0.0, 5.0),
            b[5].clamp(0.0, 5.0),
        ];
        zs.sort_by(f64::total_cmp);
        zs.windows(2)
            .filter(|w| w[0] >= b[2] && w[1] <= b[5])
            .map(|w| {
                let (minus, plus) = self.section(mirror, 0.5 * (w[0] + w[1]));
                let in_rect = |p: &[[f64; 2]]| clip(&rect, p);
                let mut rod = disk_in(&rect);
                for p in &minus {
                    rod -= disk_in(&in_rect(p));
                }
                for p in &plus {
                    let q = in_rect(p);
                    rod += area(&q) - disk_in(&q);
                }
                let here = if bore {
                    area(&in_rect(BLOCK)) - rod
                } else {
                    rod
                };
                (w[1] - w[0]) * here
            })
            .sum()
    }
}

/// The signed area of a polygon (counter-clockwise positive).
fn area(p: &[[f64; 2]]) -> f64 {
    (0..p.len())
        .map(|i| {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            a[0] * b[1] - a[1] * b[0]
        })
        .sum::<f64>()
        / 2.0
}

/// `p` clipped to the convex counter-clockwise polygon `clip_to`, one of its edges' half-planes at a
/// time — right for convex `p` and `clip_to`, which is all this file asks.
fn clip(clip_to: &[[f64; 2]], p: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut out = p.to_vec();
    for i in 0..clip_to.len() {
        let (a, b) = (clip_to[i], clip_to[(i + 1) % clip_to.len()]);
        let side = |q: [f64; 2]| (b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0]);
        let input = std::mem::take(&mut out);
        for j in 0..input.len() {
            let (u, v) = (input[j], input[(j + 1) % input.len()]);
            let (su, sv) = (side(u), side(v));
            if su >= 0.0 {
                out.push(u);
            }
            if (su >= 0.0) != (sv >= 0.0) {
                let t = su / (su - sv);
                out.push([u[0] + t * (v[0] - u[0]), u[1] + t * (v[1] - u[1])]);
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}

/// **The area of the unit disk inside the polygon `p`** (counter-clockwise; empty is 0), by Green's
/// theorem: each edge contributes the part of the triangle it spans with the centre that lies in
/// the disk — split where the edge crosses the circle, a stretch inside counts its triangle, one
/// outside the circular sector over it.
fn disk_in(p: &[[f64; 2]]) -> f64 {
    let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    let dot = |a: [f64; 2], b: [f64; 2]| a[0] * b[0] + a[1] * b[1];
    let mut total = 0.0;
    for i in 0..p.len() {
        let (a, b) = (p[i], p[(i + 1) % p.len()]);
        let d = [b[0] - a[0], b[1] - a[1]];
        // |a + t·d|² = 1: t² dd + 2t a·d + aa − 1 = 0.
        let (dd, ad, aa) = (dot(d, d), dot(a, d), dot(a, a));
        let mut ts = vec![0.0, 1.0];
        let disc = ad * ad - dd * (aa - 1.0);
        if disc > 0.0 {
            for t in [(-ad - disc.sqrt()) / dd, (-ad + disc.sqrt()) / dd] {
                if t > 0.0 && t < 1.0 {
                    ts.push(t);
                }
            }
        }
        ts.sort_by(f64::total_cmp);
        for w in ts.windows(2) {
            let at = |t: f64| [a[0] + t * d[0], a[1] + t * d[1]];
            let (u, v) = (at(w[0]), at(w[1]));
            let m = at(0.5 * (w[0] + w[1]));
            // A stretch inside the disk has its midpoint strictly inside; one touching the circle
            // at its midpoint (an edge tangent to it) lies outside but for that point.
            total += if dot(m, m) < 1.0 {
                cross(u, v) / 2.0
            } else {
                cross(u, v).atan2(dot(u, v)) / 2.0
            };
        }
    }
    total
}

/// One cell: the body count, or `s` for a result touching itself (`SelfTouchingResult`,
/// `NonManifoldResultEdge`, `NonManifoldVertex`), `r` for `NoClearRay` (the corner of a box on a
/// wall, where the result touches itself and no ray finds a clear way), `u` for `UnorderedEdges`,
/// `o` for `ObliqueCylinderCut`, `g` for `CylinderGateUndecided`. A built cell is valid, covered by
/// its mesh and, where `want` is given, of that volume.
fn cell(label: &str, m: &Model, got: BoolResult, want: Option<f64>) -> (char, f64) {
    match got {
        Ok(out) => {
            let bad = nacre_validate::validate(m);
            assert!(bad.is_empty(), "{label}: {bad:?}");
            let v: f64 = out
                .iter()
                .map(|&s| nacre_props::mass_props(m, s).expect("props").volume)
                .sum();
            if let Some(want) = want {
                assert!((v - want).abs() < 1e-9, "{label}: {v}, want {want}");
            }
            mesh_covers_faces(label, m, &out);
            let n = u32::try_from(out.len()).expect("a few bodies");
            (char::from_digit(n, 10).expect("under ten bodies"), v)
        }
        Err(BoolError::Rejected {
            reason:
                RejectReason::SelfTouchingResult
                | RejectReason::NonManifoldResultEdge
                | RejectReason::NonManifoldVertex,
            ..
        }) => ('s', f64::NAN),
        Err(BoolError::Rejected {
            reason: RejectReason::NoClearRay,
            ..
        }) => ('r', f64::NAN),
        Err(BoolError::Rejected {
            reason: RejectReason::UnorderedEdges,
            ..
        }) => ('u', f64::NAN),
        Err(BoolError::Rejected {
            reason: RejectReason::ObliqueCylinderCut,
            ..
        }) => ('o', f64::NAN),
        Err(BoolError::Rejected {
            reason: RejectReason::CylinderGateUndecided,
            ..
        }) => ('g', f64::NAN),
        Err(e) => panic!("{label}: {e:?}"),
    }
}

/// The volume of `[x0, y0, z0, x1, y1, z1]`.
fn box_volume(b: [f64; 6]) -> f64 {
    (b[3] - b[0]) * (b[4] - b[1]) * (b[5] - b[2])
}

/// **Every tool, every boolean**, as four characters a tool in the order of [`OPS`]. Each built cell
/// is valid and covered by its mesh, at its closed-form volume — or, under the roof, at the volume
/// its four operations agree on. `turn` places body and tool alike (a half turn about `x` or `y`, a
/// quarter about `z`; `0` leaves them).
fn cells(body: Body, bore: bool, mirror: bool, down: bool, turn: usize) -> String {
    let derived = !matches!(body, Body::Roof);
    let va = derived.then(|| body.within(bore, mirror, [-9.0, -9.0, -9.0, 9.0, 9.0, 9.0]));
    let mut out = Vec::with_capacity(TOOLS.len());
    for (k, b) in TOOLS.into_iter().enumerate() {
        let wants = va.map(|va| four(va, box_volume(b), body.within(bore, mirror, b)));
        let mut got = String::new();
        let mut volumes = [f64::NAN; 4];
        for (j, op) in OPS.into_iter().enumerate() {
            let mut m = Model::new();
            let a = body.build(&mut m, bore, mirror, down);
            let a = turned(&mut m, a, turn);
            let t = crate::fixtures::cuboid(&mut m, [b[0], b[1], b[2]], [b[3], b[4], b[5]]);
            let t = turned(&mut m, t, turn);
            let label = format!("bore = {bore}, mirror = {mirror}, turn {turn}, tool {k}, {op:?}");
            let result = run(&mut m, op, a, t);
            let (c, v) = cell(&label, &m, result, wants.map(|w| w[j]));
            got.push(c);
            volumes[j] = v;
        }
        if !derived && volumes.iter().all(|v| v.is_finite()) {
            // `B − A` and `A ∩ B` make up the box, and `A ∪ B` is the three parts.
            let [fuse, a_less_b, b_less_a, common] = volumes;
            assert!(
                (b_less_a + common - box_volume(b)).abs() < 1e-9,
                "tool {k}: B − A and A ∩ B are not the box"
            );
            assert!(
                (fuse - a_less_b - b_less_a - common).abs() < 1e-9,
                "tool {k}: A ∪ B is not its three parts"
            );
        }
        out.push(got);
    }
    out.join(" ")
}

/// `s` turned by `turn` (1: a half turn about `x`, 2: about `y`, 3: a quarter about `z`).
fn turned(m: &mut Model, s: Handle<Solid>, turn: usize) -> Handle<Solid> {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let (axis, deg) = match turn {
        0 => return s,
        1 => (Axis::X, 180),
        2 => (Axis::Y, 180),
        3 => (Axis::Z, 90),
        _ => unreachable!("three turns"),
    };
    let op = Operation::Transform {
        solid: s,
        isometry: Isometry::rotation(Rotation {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).expect("an angle"),
        }),
    };
    match apply(m, &op).expect("the turn") {
        OpOutput::Transform { solid } => solid,
        other => panic!("{other:?}"),
    }
}

/// The body's cells as built, and the same with body and tools turned or the rod's axis stated
/// the other way: all alike.
fn alike_placed(body: Body, bore: bool, mirror: bool, want: &str) {
    assert_eq!(cells(body, bore, mirror, false, 0), want, "as built");
    for turn in 1..=3 {
        assert_eq!(cells(body, bore, mirror, false, turn), want, "turn {turn}");
    }
    assert_eq!(
        cells(body, bore, mirror, true, 0),
        want,
        "the axis along −z"
    );
}

/// The cells of a body as built.
fn as_built(body: Body, bore: bool, mirror: bool, want: &str) {
    assert_eq!(cells(body, bore, mirror, false, 0), want);
}

#[test]
fn a_wall_changing_on_the_seams_ruling_is_read_wherever_the_rod_is_placed() {
    alike_placed(
        HALF_SLOPE,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1s11 1112 2110 2110",
    );
}

#[test]
fn a_wall_changing_on_the_seams_ruling_takes_a_boolean_as_a_bore() {
    as_built(
        HALF_SLOPE,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 s111 1121 1s01 1s01",
    );
}

#[test]
fn a_wall_changing_off_the_seam_takes_a_boolean() {
    as_built(
        HALF_SLOPE,
        false,
        true,
        "2110 1111 1111 1211 1121 2110 1111 1112 2110 2110",
    );
    as_built(
        HALF_SLOPE,
        true,
        true,
        "2110 1111 1111 1211 1112 1101 1111 1121 1s01 1s01",
    );
}

#[test]
fn a_wall_changing_the_other_way_up_takes_a_boolean() {
    as_built(
        SLOPE_HALF,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1s11 1112 2110 2110",
    );
    as_built(
        SLOPE_HALF,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 s111 1121 1s01 1s01",
    );
}

/// Tools 4 and 5 hold a ruling with both wedge walls (`u`, the unlearned pencil).
#[test]
fn kinked_wedges_changing_both_rulings_are_read_wherever_the_rod_is_placed() {
    alike_placed(
        KINKS,
        false,
        false,
        "2110 1111 1111 1211 uu1u uu1u 1s01 1111 2110 2110",
    );
}

#[test]
fn kinked_wedges_changing_both_rulings_take_a_boolean_as_a_bore() {
    as_built(
        KINKS,
        true,
        false,
        "2110 1111 1111 1211 uu1u uu0u 2110 1111 1s01 1s01",
    );
}

#[test]
fn a_run_held_by_three_stretches_takes_a_boolean() {
    as_built(
        THRICE,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1s21 1112 2110 2110",
    );
    as_built(
        THRICE,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 s112 1121 1s01 1s01",
    );
}

#[test]
fn a_run_ending_on_the_top_rim_takes_a_boolean() {
    as_built(
        TO_TOP,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1s11 1112 2110 2110",
    );
    as_built(
        TO_TOP,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 s111 1121 1s01 1s01",
    );
}

#[test]
fn a_run_from_rim_to_rim_takes_a_boolean() {
    as_built(
        RIM_TO_RIM,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1s11 1111 2110 2110",
    );
    as_built(
        RIM_TO_RIM,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 s111 1111 1s01 1s01",
    );
}

#[test]
fn a_tangent_class_holding_a_stretch_is_read_wherever_the_rod_is_placed() {
    alike_placed(
        Body::Tangent,
        false,
        false,
        "2110 1111 1111 1211 1111 2110 1111 1111 s110 2110",
    );
}

#[test]
fn a_tangent_class_holding_a_stretch_takes_a_boolean_as_a_bore() {
    as_built(
        Body::Tangent,
        true,
        false,
        "2110 1111 1111 1211 1111 1101 1111 1111 1s01 1s01",
    );
}

/// Tool 3, a slab across the axis, and tools 8 and 9, tangent to the rod, refuse as the roof's lower
/// cut does without the change (`ObliqueCylinderCut`, `CylinderGateUndecided` — todo «비스듬한 평면 ×
/// 원통», «세계 진술 없는 클래스가 낀 원통 불리언»).
#[test]
fn a_change_on_an_oblique_roof_is_read_wherever_the_rod_is_placed() {
    alike_placed(
        Body::Roof,
        false,
        false,
        "2110 1111 1111 oooo 1111 2110 1111 1111 gggg gggg",
    );
}

#[test]
fn a_change_on_an_oblique_roof_takes_a_boolean_as_a_bore() {
    as_built(
        Body::Roof,
        true,
        false,
        "2110 1111 1111 oooo 1111 1101 1111 1111 gggg gggg",
    );
}
