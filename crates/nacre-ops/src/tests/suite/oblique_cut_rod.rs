//! **A rod's end cut by a plane oblique to its axis** — the one shape whose rim is a whole ellipse,
//! so every consumer of an elliptical edge (mesh, measures, `validate`, STEP) is locked on it.
//!
//! The rod: radius 2 about `z` through `(0, 3)`, seam `+x`, `z ∈ [0, 10]`. The cut: a prism over
//! the plane `z = 5 + k·x` — its section in `zx` the region between that plane and `z = Z` over
//! `x ∈ [−3, 3]`, extruded over `y ∈ [0, 6]` — for `k = 1/2`, and `k = 2` (the cap 63° off square,
//! past where a winding read along the axis would call it misoriented).
//!
//! Derived, not copied: the rod's top rim is `z(θ) = 5 + 2k cos θ`, so the part below the plane has
//! volume `πr²·5 = 20π` (the part above too), its lateral area `2πr·5 = 20π`, the elliptical cap's
//! `πr·a` with `a = r√(1 + k²)`, and centroid `x̄ = ∬x·z_top / V = k/5`, `z̄ = ∬z_top²/2 / V =
//! 2.5 + k²/10`; the prism is `6·6·(Z − 5)`.

use super::*;

const PI: f64 = std::f64::consts::PI;
const SLOPES: [f64; 2] = [0.5, 2.0];

/// The prism's roof: clear of the rod's top and of the plane's highest point over `x ∈ [−3, 3]`.
fn roof(k: f64) -> f64 {
    12.0f64.max(5.0 + 3.0 * k + 1.0)
}

pub(crate) fn rod(m: &mut Model) -> Handle<Solid> {
    crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([0.0, 3.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        2.0,
        10.0,
    )
    .solid
}

pub(crate) fn prism(m: &mut Model, k: f64) -> Handle<Solid> {
    let top = roof(k);
    // The world `zx` frame: `u = z`, `v = x`, extruded along `+y`.
    let pts = [
        (5.0 - 3.0 * k, -3.0),
        (5.0 + 3.0 * k, 3.0),
        (top, 3.0),
        (top, -3.0),
    ]
    .map(|(u, v)| nacre_math::Point2::from_array([u, v]));
    let ring = crate::Ring2d::polygon_decimal(pts.to_vec()).expect("a quadrilateral");
    let profile = crate::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(m, Axis::Y);
    let Ok(OpOutput::Extrude { solid, .. }) = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 6.0,
        },
    ) else {
        panic!("the prism extrudes")
    };
    solid
}

/// The rod less the prism: the rod cut obliquely at `z = 5 + k·x`.
fn cut_rod(k: f64) -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let (r, p) = (rod(&mut m), prism(&mut m, k));
    let out = crate::boolean(&mut m, BoolKind::Cut, r, p).expect("the rod cuts");
    assert_eq!(out.len(), 1);
    (m, out[0])
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs().max(1.0)
}

/// Every kind in both orders builds, at the derived volumes.
#[test]
fn a_rods_end_cut_obliquely_builds_every_kind() {
    for k in SLOPES {
        let p = 36.0 * (roof(k) - 5.0);
        every_kind(
            |m: &mut Model| (rod(m), prism(m, k)),
            [
                [
                    Want::Built(1, 20.0 * PI + p),
                    Want::Built(1, 20.0 * PI),
                    Want::Built(1, 20.0 * PI),
                ],
                [
                    Want::Built(1, 20.0 * PI + p),
                    Want::Built(1, p - 20.0 * PI),
                    Want::Built(1, 20.0 * PI),
                ],
            ],
        );
    }
}

/// The cut rod is a valid solid, its top rim one whole ellipse on the oblique plane, and that
/// rim's seam vertex where the generator through `(2, 3)` meets the plane, `(2, 3, 5 + 2k)`.
#[test]
fn a_rods_oblique_end_is_one_whole_ellipse() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let violations = nacre_validate::validate(&m);
        assert!(violations.is_empty(), "{k}: {violations:?}");
        let ellipses: Vec<_> = m
            .reachable_from(&[s])
            .edges
            .iter()
            .filter_map(|&e| match m.edge_curve(e) {
                nacre_geom::Curve::Ellipse(c) if !c.is_circle() => Some((e, *c)),
                _ => None,
            })
            .collect();
        assert_eq!(ellipses.len(), 1, "{k}: one elliptical rim");
        let [v0, v1] = m.edge(ellipses[0].0).vertices;
        assert_eq!(v0, v1, "{k}: a whole rim");
        assert!(
            matches!(m.vertex(v0), nacre_topo::Vertex::OnSeam(_)),
            "{k}: its vertex is the seam's"
        );
        let p = m.vertex_point(v0).as_array();
        assert!(
            close(p[0], 2.0, 1e-15) && close(p[1], 3.0, 1e-15) && close(p[2], 5.0 + 2.0 * k, 1e-15),
            "{k}: {p:?}"
        );
        assert!((ellipses[0].1.semi_major() - 2.0 * (1.0 + k * k).sqrt()).abs() < 1e-12);
    }
}

/// Volume, area and centroid against the derived values, and a box that holds the shape.
#[test]
fn a_rods_oblique_end_measures_as_derived() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let mp = nacre_props::mass_props(&m, s).expect("props");
        assert!(
            close(mp.volume, 20.0 * PI, 1e-12),
            "{k}: volume {}",
            mp.volume
        );
        let area = 4.0 * PI + 20.0 * PI + 4.0 * (1.0 + k * k).sqrt() * PI;
        assert!(
            close(mp.area, area, 1e-12),
            "{k}: area {} vs {area}",
            mp.area
        );
        let c = nacre_props::centroid(&m, s).expect("centroid").as_array();
        let (x, z) = (k / 5.0, 2.5 + k * k / 10.0);
        assert!(
            close(c[0], x, 1e-12) && close(c[1], 3.0, 1e-12) && close(c[2], z, 1e-12),
            "{k}: {c:?}"
        );
        let (lo, hi) = nacre_props::bounds(&m, s).expect("bounds");
        let (lo, hi) = (lo.as_array(), hi.as_array());
        let top = 5.0 + 2.0 * k;
        assert!(
            lo[0] <= -2.0 && hi[0] >= 2.0 && lo[2] <= 0.0 && hi[2] >= top,
            "{lo:?} {hi:?}"
        );
        assert!(
            hi[2] <= top + 1e-12 && lo[1] >= 1.0 - 1e-12,
            "{lo:?} {hi:?}"
        );
    }
}

/// The mesh: no lateral triangle leans off the cylinder, the elliptical rim's polyline turns within
/// the budget at its major ends too (an ellipse turns `a/b` times as fast as its circle there), and
/// the triangles enclose the solid's volume to the mesh's resolution.
#[test]
fn a_rods_oblique_end_meshes() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let cfg = nacre_tess::TessConfig::default();
        let t = nacre_tess::tessellate(&m, &cfg).expect("the mesh");
        let leaning = nacre_tess::leaning_laterals(&m, &t, cfg.max_angle_deg);
        assert!(leaning.is_empty(), "{leaning:?}");
        let kinked = nacre_tess::kinked_edges(&m, &t, cfg.max_angle_deg);
        assert!(kinked.is_empty(), "{k}: {kinked:?}");
        let unshared = nacre_tess::unshared_columns(&m, &t);
        assert!(unshared.is_empty(), "{k}: {unshared:?}");
        let mut volume = 0.0;
        for f in &m.reachable_from(&[s]).faces {
            for &tri in &t.by_face[f] {
                let [a, b, c] = t.triangles.get(tri).vertices.map(|v| t.vertices.get(v).pos);
                let o = Point3::origin();
                volume += (a - o).dot((b - o).cross(c - o)) / 6.0;
            }
        }
        assert!(close(volume, 20.0 * PI, 2e-3), "{k}: mesh volume {volume}");
    }
}

/// STEP writes the rim as an `ELLIPSE` and proves an uncertainty for it.
#[test]
fn a_rods_oblique_end_writes_an_ellipse() {
    for k in SLOPES {
        let (m, s) = cut_rod(k);
        let text = nacre_step::to_step(&m, "2000-01-01T00:00:00").expect("STEP");
        assert!(text.contains("ELLIPSE("), "{k}: no ELLIPSE written");
        match nacre_step::uncertainty(&m, &[s]) {
            nacre_step::Uncertainty::Proven(d) => assert!(d < 1e-12, "{k}: uncertainty {d}"),
            other => panic!("{k}: not proven: {other:?}"),
        }
    }
}

/// **A tool fed against the cut rod** (`k = 1/2`) — the re-entry scoreboard's rows. `volume` is the
/// tool's and `common` its common part with the rod, both derived; with `R = 20π` the rod's, a cell
/// that builds has `Fuse = R + T − C`, `rod − tool = R − C`, `tool − rod = T − C`.
struct Tool {
    name: &'static str,
    build: fn(&mut Model) -> Handle<Solid>,
    volume: f64,
    common: f64,
    /// A plane of the tool's faces touches the rod's slanted rim at one point — its highest,
    /// `(2, 3, 6)`, or its lowest, `(−2, 3, 4)` — whether or not the face reaches it.
    touches_rim: bool,
    /// A plane of the tool's faces across the axis crosses the slanted rim (`4 < z < 6`).
    crosses_rim: bool,
}

fn rod_on(m: &mut Model, at: [f64; 2], r: f64, z: [f64; 2]) -> Handle<Solid> {
    crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([at[0], at[1], z[0]]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        r,
        z[1] - z[0],
    )
    .solid
}

/// The tools, each against the rod `x² + (y − 3)² ≤ 4`, `0 ≤ z ≤ 5 + x/2`. The commons integrate the
/// rod's height over the part of its disk the tool covers: above `z = c` the rod reaches only where
/// `x > 2c − 10`, a segment of the disk.
fn tools() -> Vec<Tool> {
    let s3 = 3.0f64.sqrt();
    vec![
        Tool {
            name: "a coaxial bore of radius 1",
            build: |m| rod_on(m, [0.0, 3.0], 1.0, [-1.0, 11.0]),
            volume: 12.0 * PI,
            common: 5.0 * PI,
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a box over the rod's foot",
            build: |m| crate::fixtures::cuboid(m, [1.0, 1.0, -1.0], [4.0, 5.0, 2.0]),
            volume: 36.0,
            common: 2.0 * (4.0 * PI / 3.0 - s3),
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a second slant, z = 5 − x/2 and above",
            build: |m| prism(m, -0.5),
            volume: 252.0,
            common: 16.0 / 3.0,
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a box above z = 5.5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 5.5], [3.0, 6.0, 8.0]),
            volume: 90.0,
            common: 1.5 * s3 - 2.0 * PI / 3.0,
            touches_rim: false,
            crosses_rim: true,
        },
        Tool {
            name: "a box above z = 4.5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 4.5], [3.0, 6.0, 8.0]),
            volume: 126.0,
            common: 4.0 * PI / 3.0 + 1.5 * s3,
            touches_rim: false,
            crosses_rim: true,
        },
        Tool {
            name: "a box above z = 5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 5.0], [3.0, 6.0, 8.0]),
            volume: 108.0,
            common: 8.0 / 3.0,
            touches_rim: false,
            crosses_rim: true,
        },
        Tool {
            name: "a pin of radius 1/2 through (1, 3)",
            build: |m| rod_on(m, [1.0, 3.0], 0.5, [-1.0, 11.0]),
            volume: 3.0 * PI,
            common: 11.0 * PI / 8.0,
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a pin of radius 1/2 sunk on the axis",
            build: |m| rod_on(m, [0.0, 3.0], 0.5, [4.0, 8.0]),
            volume: PI,
            common: PI / 4.0,
            touches_rim: true,
            crosses_rim: false,
        },
        Tool {
            name: "a wall box from x = 1.5 up to z = 5.5",
            build: |m| crate::fixtures::cuboid(m, [1.5, 0.0, -1.0], [4.0, 6.0, 5.5]),
            volume: 97.5,
            common: 22.0 * 0.75f64.acos() - 33.0 * 7.0f64.sqrt() / 8.0,
            touches_rim: false,
            crosses_rim: true,
        },
        Tool {
            name: "the slant's own prism",
            build: |m| prism(m, 0.5),
            volume: 252.0,
            common: 0.0,
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a far box",
            build: |m| crate::fixtures::cuboid(m, [10.0, 10.0, 10.0], [11.0, 11.0, 11.0]),
            volume: 1.0,
            common: 0.0,
            touches_rim: false,
            crosses_rim: false,
        },
        Tool {
            name: "a box touching the rim's top at z = 6",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 6.0], [3.0, 6.0, 8.0]),
            volume: 72.0,
            common: 0.0,
            touches_rim: true,
            crosses_rim: false,
        },
    ]
}

/// A boolean's entry: the front door, or the stages behind it.
type Entry =
    fn(&mut Model, BoolKind, Handle<Solid>, Handle<Solid>) -> Result<Vec<Handle<Solid>>, BoolError>;

/// A refusal a cell is held to: by its name, or — a tracer's decline — by the decline's kind (the
/// face it names is the model's handle, which the cell does not state).
#[derive(Clone, Copy, Debug)]
enum Refusal {
    Name(RejectReason),
    Decline(crate::DeclineKind),
}

impl Refusal {
    fn is(self, reason: RejectReason) -> bool {
        match (self, reason) {
            (Self::Name(name), r) => r == name,
            (Self::Decline(kind), RejectReason::TraceDeclined { kind: k, .. }) => k == kind,
            (Self::Decline(_), _) => false,
        }
    }
}

/// Every cell of every tool against `operand` — each kind with the operand first, then swapped —
/// through `run`, refused as `want(tool)` says, raised in a file ending as it says.
fn every_cell_refused(
    operand: fn() -> (Model, Handle<Solid>),
    run: Entry,
    want: impl Fn(&Tool) -> (Refusal, &'static str) + Sync,
) {
    use crate::error::raise_probe::RAISED;
    crate::ledger::owned(|| {
        for tool in tools() {
            for swap in [false, true] {
                for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                    let at = format!("{} — {kind:?} swap={swap}", tool.name);
                    let (name, file) = want(&tool);
                    let (mut m, s) = operand();
                    let t = (tool.build)(&mut m);
                    let (x, y) = if swap { (t, s) } else { (s, t) };
                    let r0 = RAISED.len();
                    match run(&mut m, kind, x, y) {
                        Err(BoolError::Rejected { reason, .. }) if name.is(reason) => {
                            let sites: Vec<&str> =
                                RAISED.mine_since(r0).into_iter().map(|(_, f)| f).collect();
                            assert!(
                                sites.iter().any(|f| f.ends_with(file)),
                                "{at}: raised at {sites:?}, not {file}"
                            );
                        }
                        out => panic!("{at}: want {name:?}, got {out:?}"),
                    }
                }
            }
        }
    });
}

/// **The scoreboard's commons are what the forward road builds** — every operand on it has no
/// elliptical edge, so the door admits it: a flat tool is cut by the slant first and then met with
/// the rod; a tool with a cylinder is met with the rod first (the slant would leave an ellipse on
/// it) and the piece then cut by the slant. Either way the volume is the derived `common`, so a cell
/// that builds behind the door is held to a number a road already answers, not to arithmetic alone.
#[test]
fn the_scoreboards_commons_are_what_the_forward_road_builds() {
    let solids = |m: &Model, out: Vec<Handle<Solid>>| -> f64 {
        out.into_iter()
            .map(|s| nacre_props::mass_props(m, s).expect("props").volume)
            .sum()
    };
    for tool in tools() {
        let mut m = Model::new();
        let r = rod(&mut m);
        let t = (tool.build)(&mut m);
        let tv = nacre_props::mass_props(&m, t).expect("props").volume;
        assert!(close(tv, tool.volume, 1e-12), "{}: tool {tv}", tool.name);
        let flat = m.reachable_from(&[t]).faces.iter().all(|&f| {
            matches!(
                m.surface(m.face(f).surface),
                nacre_topo::Surface::Plane { .. }
            )
        });
        let p = prism(&mut m, 0.5);
        let got = if flat {
            let below = crate::boolean(&mut m, BoolKind::Cut, t, p)
                .unwrap_or_else(|e| panic!("{}: tool − slant {e:?}", tool.name));
            let mut v = 0.0;
            let mut r = Some(r);
            for s in below {
                // Each piece meets a fresh rod: a boolean consumes its operands.
                let rod_here = r.take().unwrap_or_else(|| rod(&mut m));
                let out = crate::boolean(&mut m, BoolKind::Common, rod_here, s)
                    .unwrap_or_else(|e| panic!("{}: rod ∩ piece {e:?}", tool.name));
                v += solids(&m, out);
            }
            v
        } else {
            let piece = crate::boolean(&mut m, BoolKind::Common, r, t)
                .unwrap_or_else(|e| panic!("{}: rod ∩ tool {e:?}", tool.name));
            let mut v = 0.0;
            let mut p = Some(p);
            for s in piece {
                let slant = p.take().unwrap_or_else(|| prism(&mut m, 0.5));
                let out = crate::boolean(&mut m, BoolKind::Cut, s, slant)
                    .unwrap_or_else(|e| panic!("{}: piece − slant {e:?}", tool.name));
                v += solids(&m, out);
            }
            v
        };
        assert!(
            close(got, tool.common, 1e-9),
            "{}: {got} vs {}",
            tool.name,
            tool.common
        );
    }
}

/// **A result with an elliptical edge is refused as an operand, at the door** — every tool, every
/// kind, both orders: `ObliqueCylinderCut`, raised by the boolean's entrance (the stages behind it
/// read a curve on a cylinder as a circle across its axis), not by any of them.
#[test]
fn a_cut_rod_fed_back_is_refused_at_the_door() {
    every_cell_refused(
        || cut_rod(0.5),
        crate::boolean,
        |_| {
            (
                Refusal::Name(RejectReason::ObliqueCylinderCut),
                "nacre-ops/src/boolean.rs",
            )
        },
    );
}

/// **What the stages behind the door answer for the cut rod** — the scoreboard the re-entry work
/// flips row by row, for the rod and the rod with a window in its side. A tool with a plane touching
/// the slanted rim at a point stops at the cylinder gate (`ClassMeet::Undrawable`); every other cell
/// passes it — the rod's own cap reaches its lateral and is recorded as its section — and the
/// tracer. A tool with a plane across the axis crossing the slanted rim stops at a straight angle
/// past the arrangement (`StraightAngle`); the rest stop at the lateral's band rows, which read the face's ends by
/// their heights, and a slanted rim has none (`DeclineKind::CylSpan`).
#[test]
fn a_cut_rod_behind_the_door_stops_at_the_tracer() {
    for operand in [|| cut_rod(0.5), windowed_cut_rod] {
        every_cell_refused(
            operand,
            crate::boolean::behind_the_elliptical_door,
            |tool| {
                if tool.touches_rim {
                    (
                        Refusal::Name(RejectReason::ObliqueCylinderCut),
                        "nacre-ops/src/planes/cyl_gate.rs",
                    )
                } else if tool.crosses_rim {
                    // By name alone: the classes are read in parallel past a first refusal and one
                    // at a time in the serial build, so the guard that rings first differs.
                    (Refusal::Name(RejectReason::StraightAngle), "")
                } else {
                    (
                        Refusal::Decline(crate::DeclineKind::CylSpan),
                        "nacre-ops/src/bands.rs",
                    )
                }
            },
        );
    }
}

/// The cut rod with a window through its side first (`x ≥ 1`, `y ∈ [2, 4]`, `z ∈ [2, 3]`): its
/// lateral carries two circles at `z = 2, 3` besides its rims, so the circles alone give a range
/// (`[0, 3]`) — which stops short of the slanted rim.
fn windowed_cut_rod() -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let r = rod(&mut m);
    let w = crate::fixtures::cuboid(&mut m, [1.0, 2.0, 2.0], [3.0, 4.0, 3.0]);
    let r = crate::boolean(&mut m, BoolKind::Cut, r, w).expect("the window")[0];
    let p = prism(&mut m, 0.5);
    let out = crate::boolean(&mut m, BoolKind::Cut, r, p).expect("the slant cuts");
    assert_eq!(out.len(), 1);
    (m, out[0])
}

/// **A moved cut rod keeps its edges' kinds** — translated by a third and turned a quarter about
/// `z`, it validates, keeps one elliptical rim (re-derived from its moved carriers, the kind theirs),
/// and its volume.
#[test]
fn a_moved_cut_rod_keeps_its_elliptical_rim() {
    use nacre_exact::{Isometry, Rat};
    let third = Rat::new(1, 3).expect("a third");
    let zero = Rat::from_int(0);
    let (mut m, mut s) = cut_rod(0.5);
    for iso in [
        Isometry::translation([third, zero, third]),
        super::fixtures::rot_iso(Axis::Z, 90),
    ] {
        s = super::fixtures::xf(&mut m, s, iso);
    }
    let violations = nacre_validate::validate(&m);
    assert!(violations.is_empty(), "{violations:?}");
    let ellipses = m
        .reachable_from(&[s])
        .edges
        .iter()
        .filter(|&&e| matches!(m.edge_curve(e), nacre_geom::Curve::Ellipse(c) if !c.is_circle()))
        .count();
    assert_eq!(ellipses, 1);
    let v = nacre_props::mass_props(&m, s).expect("props").volume;
    assert!(close(v, 20.0 * PI, 1e-12), "{v}");
}
