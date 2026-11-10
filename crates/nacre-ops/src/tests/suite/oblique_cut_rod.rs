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
    zx_prism(
        m,
        &[
            (5.0 - 3.0 * k, -3.0),
            (5.0 + 3.0 * k, 3.0),
            (top, 3.0),
            (top, -3.0),
        ],
    )
}

/// A polygon `(z, x)` in the world `zx` frame (`u = z`, `v = x`) extruded along `+y` over `[0, 6]`.
pub(super) fn zx_prism(m: &mut Model, pts: &[(f64, f64)]) -> Handle<Solid> {
    let pts = pts
        .iter()
        .map(|&(u, v)| nacre_math::Point2::from_array([u, v]))
        .collect();
    let ring = crate::Ring2d::polygon_decimal(pts).expect("a polygon");
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
    /// The bodies a cell that builds leaves — `[Fuse, Cut, Common]`, the rod first, then the tool.
    bodies: [[usize; 3]; 2],
    /// The tool's volume inside the windowed rod's window ([`windowed_cut_rod`]) — the part of
    /// `common` that rod lacks.
    in_window: f64,
    /// That part is a body of its own in `tool − rod` (rod material above and below it).
    window_piece: bool,
    /// The tool's surface touches the window's face `x = 1` along a line: the windowed rod less the
    /// tool touches itself there.
    touches_window: bool,
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
            bodies: [[1, 1, 1], [1, 2, 1]],
            in_window: 0.0,
            window_piece: false,
            touches_window: true,
        },
        Tool {
            name: "a box over the rod's foot",
            build: |m| crate::fixtures::cuboid(m, [1.0, 1.0, -1.0], [4.0, 5.0, 2.0]),
            volume: 36.0,
            common: 2.0 * (4.0 * PI / 3.0 - s3),
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a second slant, z = 5 − x/2 and above",
            build: |m| prism(m, -0.5),
            volume: 252.0,
            common: 16.0 / 3.0,
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a box above z = 5.5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 5.5], [3.0, 6.0, 8.0]),
            volume: 90.0,
            common: 1.5 * s3 - 2.0 * PI / 3.0,
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a box above z = 4.5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 4.5], [3.0, 6.0, 8.0]),
            volume: 126.0,
            common: 4.0 * PI / 3.0 + 1.5 * s3,
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a box above z = 5",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 5.0], [3.0, 6.0, 8.0]),
            volume: 108.0,
            common: 8.0 / 3.0,
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a pin of radius 1/2 through (1, 3)",
            build: |m| rod_on(m, [1.0, 3.0], 0.5, [-1.0, 11.0]),
            volume: 3.0 * PI,
            common: 11.0 * PI / 8.0,
            touches_rim: false,
            bodies: [[1, 1, 1], [1, 2, 1]],
            in_window: PI / 8.0,
            window_piece: true,
            touches_window: false,
        },
        Tool {
            name: "a pin of radius 1/2 sunk on the axis",
            build: |m| rod_on(m, [0.0, 3.0], 0.5, [4.0, 8.0]),
            volume: PI,
            common: PI / 4.0,
            touches_rim: true,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a wall box from x = 1.5 up to z = 5.5",
            build: |m| crate::fixtures::cuboid(m, [1.5, 0.0, -1.0], [4.0, 6.0, 5.5]),
            volume: 97.5,
            common: 22.0 * 0.75f64.acos() - 33.0 * 7.0f64.sqrt() / 8.0,
            touches_rim: false,
            bodies: [[1; 3]; 2],
            in_window: 2.0 * PI / 3.0 + s3 - 3.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "the slant's own prism",
            build: |m| prism(m, 0.5),
            volume: 252.0,
            common: 0.0,
            touches_rim: false,
            bodies: [[1, 1, 0], [1, 1, 0]],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a far box",
            build: |m| crate::fixtures::cuboid(m, [10.0, 10.0, 10.0], [11.0, 11.0, 11.0]),
            volume: 1.0,
            common: 0.0,
            touches_rim: false,
            bodies: [[2, 1, 0], [2, 1, 0]],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
        Tool {
            name: "a box touching the rim's top at z = 6",
            build: |m| crate::fixtures::cuboid(m, [-3.0, 0.0, 6.0], [3.0, 6.0, 8.0]),
            volume: 72.0,
            common: 0.0,
            touches_rim: true,
            bodies: [[1; 3]; 2],
            in_window: 0.0,
            window_piece: false,
            touches_window: false,
        },
    ]
}

/// One cell's answer: refused by this name, raised in a file ending as the `str` says (empty: any),
/// or built into this many bodies.
#[derive(Clone, Copy, Debug)]
enum Cell {
    Refused(RejectReason, &'static str),
    Built(usize),
}

/// The cut rod, plain or with its window ([`windowed_cut_rod`]), and its volume.
fn operand(windowed: bool) -> (Model, Handle<Solid>, f64) {
    let rod = 20.0 * PI;
    if windowed {
        let (m, s) = windowed_cut_rod();
        (m, s, rod - (3f64.sqrt() + 2.0 * PI / 3.0 - 2.0))
    } else {
        let (m, s) = cut_rod(0.5);
        (m, s, rod)
    }
}

/// Every cell of every tool against the rod — each kind with the rod first, then swapped — held to
/// `want(tool)`: a refusal to its name and raise file; a build to its bodies,
/// `validate`, the mesh, its derived volume (relative `1e-9`) and the mesh's enclosed volume (`2e-3`).
fn scoreboard(windowed: bool, want: impl Fn(&Tool) -> [[Cell; 3]; 2] + Sync) {
    use crate::error::raise_probe::RAISED;
    crate::ledger::owned(|| {
        for tool in tools() {
            for (swap, row) in [false, true].into_iter().zip(want(&tool)) {
                for (kind, cell) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
                    .into_iter()
                    .zip(row)
                {
                    let at = format!("{} — {kind:?} swap={swap} windowed={windowed}", tool.name);
                    let (mut m, s, rod) = operand(windowed);
                    let t = (tool.build)(&mut m);
                    let (x, y) = if swap { (t, s) } else { (s, t) };
                    let r0 = RAISED.len();
                    match (cell, crate::boolean(&mut m, kind, x, y)) {
                        (Cell::Refused(name, file), Err(BoolError::Rejected { reason, .. }))
                            if reason == name =>
                        {
                            let sites: Vec<&str> =
                                RAISED.mine_since(r0).into_iter().map(|(_, f)| f).collect();
                            assert!(
                                sites.iter().any(|f| f.ends_with(file)),
                                "{at}: raised at {sites:?}, not {file}"
                            );
                        }
                        (Cell::Built(n), Ok(out)) => {
                            let got = built_volume(&m, &out, n, &at);
                            let c = tool.common - if windowed { tool.in_window } else { 0.0 };
                            let derived = match (kind, swap) {
                                (BoolKind::Fuse, _) => rod + tool.volume - c,
                                (BoolKind::Cut, false) => rod - c,
                                (BoolKind::Cut, true) => tool.volume - c,
                                (BoolKind::Common, _) => c,
                            };
                            assert!(close(got, derived, 1e-9), "{at}: {got} vs {derived}");
                            let mesh = mesh_volume(&m, &out);
                            assert!(close(mesh, derived, 2e-3), "{at}: mesh {mesh} vs {derived}");
                        }
                        (cell, out) => panic!("{at}: want {cell:?}, got {out:?}"),
                    }
                }
            }
        }
    });
}

/// The volume the mesh's triangles enclose over `solids`' faces.
pub(super) fn mesh_volume(m: &Model, solids: &[Handle<Solid>]) -> f64 {
    let t = nacre_tess::tessellate(m, &nacre_tess::TessConfig::default()).expect("the mesh");
    let o = Point3::origin();
    m.reachable_from(solids)
        .faces
        .iter()
        .flat_map(|f| &t.by_face[f])
        .map(|&tri| {
            let [a, b, c] = t.triangles.get(tri).vertices.map(|v| t.vertices.get(v).pos);
            (a - o).dot((b - o).cross(c - o)) / 6.0
        })
        .sum()
}

/// **The scoreboard's commons are what a second road builds** — one on which no operand carries a
/// slanted rim: a flat tool is cut by the slant first and then met with the rod; a tool with a
/// cylinder is met with the rod first (the slant would leave an ellipse on it) and the piece then
/// cut by the slant. Either way the volume is the derived `common` — less the part inside the window,
/// for the rod with its window — so a cell that builds is held to a number a road that never reads a
/// slanted rim answers, not to arithmetic alone.
#[test]
fn the_scoreboards_commons_are_what_the_forward_road_builds() {
    let solids = |m: &Model, out: Vec<Handle<Solid>>| -> f64 {
        out.into_iter()
            .map(|s| nacre_props::mass_props(m, s).expect("props").volume)
            .sum()
    };
    for windowed in [false, true] {
        // The rod as the slant finds it: with its window cut first, or plain.
        let body = |m: &mut Model| {
            let r = rod(m);
            if !windowed {
                return r;
            }
            let w = crate::fixtures::cuboid(m, [1.0, 2.0, 2.0], [3.0, 4.0, 3.0]);
            crate::boolean(m, BoolKind::Cut, r, w).expect("the window")[0]
        };
        for tool in tools() {
            let mut m = Model::new();
            let r = body(&mut m);
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
                    let rod_here = r.take().unwrap_or_else(|| body(&mut m));
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
            let want = tool.common - if windowed { tool.in_window } else { 0.0 };
            assert!(
                close(got, want, 1e-9),
                "{} (windowed {windowed}): {got} vs {want}",
                tool.name
            );
        }
    }
}

/// **A cut rod fed back to a boolean** — the rod and the rod with a window in its side, every tool,
/// every kind, both orders. A tool with a plane touching the slanted rim at a point stops at the
/// cylinder gate (`ClassMeet::Undrawable`). Every other cell builds, at its derived volume — a plane
/// across the axis crossing the slanted rim included, and the wall box that seals the window's rest
/// into a cavity, in either order — but the windowed rod less a tool touching its window along a
/// line, which touches itself (`SelfTouchingResult`).
#[test]
fn a_cut_rod_fed_back_builds_but_where_a_plane_touches_its_slanted_rim() {
    for windowed in [false, true] {
        scoreboard(windowed, |tool| {
            if tool.touches_rim {
                let gate = Cell::Refused(
                    RejectReason::ObliqueCylinderCut,
                    "nacre-ops/src/planes/cyl_gate.rs",
                );
                return [[gate; 3]; 2];
            }
            let mut cells = tool.bodies.map(|row| row.map(Cell::Built));
            if windowed && tool.window_piece {
                cells[1][1] = Cell::Built(tool.bodies[1][1] + 1);
            }
            if windowed && tool.touches_window {
                cells[0][1] = Cell::Refused(
                    RejectReason::SelfTouchingResult,
                    "nacre-ops/src/assembly/self_touch.rs",
                );
            }
            cells
        });
    }
}

/// **Other bodies with slanted rims fed back** — five more shapes against four tools, every kind,
/// both orders, each cell at its derived volume (relative `1e-9`), validated and meshed:
/// - a plate `[−5, 5]² × [0, 4]` drilled along `(0, 3, 4)` with radius 1 — both rims slanted, the
///   material outside the bore — of volume `400 − 5π` (the bore crosses the slab `5/4` times its
///   thickness);
/// - the rod halved at `y = 3` then cut by the slant — a lateral whose slanted rim is an arc ending on
///   two rulings, the cap a half-ellipse — of volume `10π`;
/// - the rod notched on both sides (`zx_prism`): `x ≤ −3/2` over `1 ≤ z ≤ 5 + x`, and `x ≥ 3/2` over
///   `11/2 − x ≤ z ≤ 9` — two slanted notch floors that cross on the cylinder at `x = 1/4`, between
///   the lateral's whole rims `z = 0` and `z = 10`, each still read as below or above the face by the
///   rest of its boundary (`arrangement`'s `rim_side`). With `S = 4·acos(3/4) − (3/4)√7` the area of the disk past
///   `|x| = 3/2` and `7√7/12` its first moment, the notches take `(4S − 7√7/12) + (7S/2 + 7√7/12)`,
///   so the body is `40π − 15S/2`. The pin is refused in every cell by the nesting of a class's rings,
///   which finds no witness (`RingHasNoWitness`, `nesting`'s `inside`): the notch floors' planes
///   section the pin's cylinder clear of the floors — the first notch's prism against the pin builds
///   in every cell;
/// - the rod notched by the wedge between `z = 4 + x` and `z = 6 − x` over `x ≤ 1` — two slanted
///   notch faces meeting on an edge `x = 1, z = 5` that crosses the lateral — of volume
///   `104π/3 − 6√3` (the wedge's height `2 − 2x` over the disk short of `x = 1`). Past the far box
///   every cell is refused: `z = 2` touches the wedge floor's section at its lowest point (the gate,
///   `ObliqueCylinderCut`), the box `x ≥ 1` — whose ruling runs through the lateral's corner where
///   the notch's two sections meet — stops at the ruling sweep (`DeclineKind::Ruling`, as the
///   crossing-ends ridges in `crossing_ends.rs` do), and the pin as for the two notches;
/// - the rod with a window through `x ≤ −3/2` between `z = 4 + x` and `z = 6 − x` — a floor and a
///   ceiling crossing on the cylinder at `x = 1`, the window's sides rulings — of volume
///   `40π − 2S − 7√7/6`; its rulings run between the two sections, which order only on a generator,
///   so the operand's boundary is read as curved (`CurvedOperandBoundary`) and the trace declines
///   a face of the window (`DeclineKind::OuterRing`) in every cell but the touch at `z = 2`.
///
/// The tools: a far box, the half-space box `z ≥ 2`, the half-space box `x ≥ 1`, and the pin of
/// radius ½ about `(0, 3)`. The plate's bore touches `x = 1` along a line, so the plate less that box
/// touches itself; its bore and the pin are two cylinders of crossing axes, which the gate refuses.
#[test]
fn other_bodies_with_slanted_rims_fed_back() {
    use RejectReason::{
        CylinderPairContact, ObliqueCylinderCut, RingHasNoWitness, SelfTouchingResult,
    };
    let (s3, pi) = (3f64.sqrt(), PI);
    let plate = |m: &mut Model| {
        let p = crate::fixtures::cuboid(m, [-5.0, -5.0, 0.0], [5.0, 5.0, 4.0]);
        let d = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([0.0, -2.25, -1.0]),
            Vector3::from_array([0.0, 0.6, 0.8]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            1.0,
            7.5,
        )
        .solid;
        crate::boolean(m, BoolKind::Cut, p, d).expect("the drilled plate")[0]
    };
    let half = |m: &mut Model| {
        let r = rod(m);
        let b = crate::fixtures::cuboid(m, [-3.0, -1.0, -1.0], [3.0, 3.0, 11.0]);
        let h = crate::boolean(m, BoolKind::Cut, r, b).expect("the half")[0];
        let p = prism(m, 0.5);
        crate::boolean(m, BoolKind::Cut, h, p).expect("the slanted half")[0]
    };
    let notched = |m: &mut Model| {
        let mut r = rod(m);
        for pts in [
            &[(1.0, -3.0), (1.0, -1.5), (3.5, -1.5), (2.0, -3.0)][..],
            &[(4.0, 1.5), (2.5, 3.0), (9.0, 3.0), (9.0, 1.5)][..],
        ] {
            let n = zx_prism(m, pts);
            r = crate::boolean(m, BoolKind::Cut, r, n).expect("a notch")[0];
        }
        r
    };
    let wedged = |m: &mut Model| {
        let r = rod(m);
        let w = zx_prism(m, &[(1.0, -3.0), (5.0, 1.0), (9.0, -3.0)]);
        crate::boolean(m, BoolKind::Cut, r, w).expect("the wedge")[0]
    };
    let ruled = |m: &mut Model| {
        let r = rod(m);
        let w = zx_prism(m, &[(1.0, -3.0), (2.5, -1.5), (7.5, -1.5), (9.0, -3.0)]);
        crate::boolean(m, BoolKind::Cut, r, w).expect("the window")[0]
    };
    type Build = fn(&mut Model) -> Handle<Solid>;
    let tools: [(Build, f64); 4] = [
        (
            |m| crate::fixtures::cuboid(m, [20.0, 20.0, 20.0], [21.0, 21.0, 21.0]),
            1.0,
        ),
        (
            |m| crate::fixtures::cuboid(m, [-6.0, -6.0, 2.0], [6.0, 6.0, 12.0]),
            1440.0,
        ),
        (
            |m| crate::fixtures::cuboid(m, [1.0, -6.0, -2.0], [6.0, 6.0, 12.0]),
            840.0,
        ),
        (|m| rod_on(m, [0.0, 3.0], 0.5, [-1.0, 11.0]), 3.0 * pi),
    ];
    // Per body: its volume, and per tool the common part — `Err` where every cell refuses by that
    // name — with the cells' bodies `[Fuse, Cut, Common]`, the body first then the tool, and the one
    // cell (kind, swap) refused by a name of its own.
    type Row = (
        Result<f64, RejectReason>,
        [[usize; 3]; 2],
        Option<(BoolKind, bool, RejectReason)>,
    );
    let far: Row = (Ok(0.0), [[2, 1, 0], [2, 1, 0]], None);
    let ones = [[1; 3]; 2];
    let outer_ring = RejectReason::TraceDeclined {
        kind: crate::DeclineKind::OuterRing,
        face: None,
    };
    // The disk past `|x| = 3/2`: its area and first moment.
    let (seg, moment) = (
        4.0 * 0.75f64.acos() - 0.75 * 7f64.sqrt(),
        7.0 * 7f64.sqrt() / 12.0,
    );
    let bodies: [(Build, f64, [Row; 4]); 5] = [
        (
            plate,
            400.0 - 5.0 * pi,
            [
                far,
                (Ok(200.0 - 2.5 * pi), ones, None),
                (
                    Ok(160.0),
                    ones,
                    Some((BoolKind::Cut, false, SelfTouchingResult)),
                ),
                (Err(ObliqueCylinderCut), ones, None),
            ],
        ),
        (
            half,
            10.0 * pi,
            [
                far,
                (Ok(6.0 * pi), ones, None),
                (Ok(5.0 * (2.0 * pi / 3.0 - s3 / 2.0) + s3 / 2.0), ones, None),
                (Ok(5.0 * pi / 8.0), ones, None),
            ],
        ),
        (
            notched,
            40.0 * pi - 7.5 * seg,
            [
                far,
                (Ok(32.0 * pi - 6.5 * seg), ones, None),
                (
                    Ok(10.0 * (4.0 * pi / 3.0 - s3) - 3.5 * seg - moment),
                    ones,
                    None,
                ),
                (Err(RingHasNoWitness), ones, None),
            ],
        ),
        (
            wedged,
            104.0 * pi / 3.0 - 6.0 * s3,
            [
                far,
                (Err(ObliqueCylinderCut), ones, None),
                (
                    Err(RejectReason::TraceDeclined {
                        kind: crate::DeclineKind::Ruling,
                        face: None,
                    }),
                    ones,
                    None,
                ),
                (Err(RingHasNoWitness), ones, None),
            ],
        ),
        (
            ruled,
            40.0 * pi - 2.0 * seg - 2.0 * moment,
            [
                (Err(outer_ring), far.1, None),
                (Err(ObliqueCylinderCut), ones, None),
                (Err(outer_ring), ones, None),
                (Err(outer_ring), ones, None),
            ],
        ),
    ];
    // A trace's decline is named by its kind; the face it names is the body's own handle.
    let same_name = |a: RejectReason, b: RejectReason| match (a, b) {
        (
            RejectReason::TraceDeclined { kind: x, .. },
            RejectReason::TraceDeclined { kind: y, .. },
        ) => x == y,
        _ => a == b,
    };
    for (body, volume, rows) in bodies {
        let mut m = Model::new();
        let b = body(&mut m);
        let got = nacre_props::mass_props(&m, b).expect("props").volume;
        assert!(close(got, volume, 1e-12), "the body: {got} vs {volume}");
        for ((tool, tool_volume), (common, counts, refused)) in tools.iter().zip(rows) {
            for (swap, row) in [false, true].into_iter().zip(counts) {
                for (kind, n) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
                    .into_iter()
                    .zip(row)
                {
                    let at = format!(
                        "body of {volume:.3}, tool of {tool_volume:.3} — {kind:?} swap={swap}"
                    );
                    let mut m = Model::new();
                    let (b, t) = (body(&mut m), tool(&mut m));
                    let (x, y) = if swap { (t, b) } else { (b, t) };
                    let out = crate::boolean(&mut m, kind, x, y);
                    let own = refused.filter(|&(k, s, _)| k == kind && s == swap);
                    match (common, out) {
                        (Err(name), Err(BoolError::Rejected { reason, .. })) => {
                            assert!(
                                same_name(reason, name) || reason == CylinderPairContact,
                                "{at}: {reason:?}"
                            );
                        }
                        (Ok(_), Err(BoolError::Rejected { reason, .. })) if own.is_some() => {
                            assert_eq!(Some(reason), own.map(|(_, _, r)| r), "{at}");
                        }
                        (Ok(c), Ok(out)) => {
                            let got = built_volume(&m, &out, n, &at);
                            let derived = match (kind, swap) {
                                (BoolKind::Fuse, _) => volume + tool_volume - c,
                                (BoolKind::Cut, false) => volume - c,
                                (BoolKind::Cut, true) => tool_volume - c,
                                (BoolKind::Common, _) => c,
                            };
                            assert!(close(got, derived, 1e-9), "{at}: {got} vs {derived}");
                            let mesh = mesh_volume(&m, &out);
                            assert!(close(mesh, derived, 2e-3), "{at}: mesh {mesh} vs {derived}");
                        }
                        (want, got) => panic!("{at}: want {want:?}, got {got:?}"),
                    }
                }
            }
        }
    }
    // The first notch's prism itself against the pin — its floor's plane sections the pin's cylinder
    // clear of the floor, as on the notched rod, and the two bodies stand apart: every cell builds.
    // The prism is `6 · 3/2 · (1 + 5/2)/2`.
    let (notch, prism_volume) = (
        |m: &mut Model| zx_prism(m, &[(1.0, -3.0), (1.0, -1.5), (3.5, -1.5), (2.0, -3.0)]),
        15.75,
    );
    let (pin, pin_volume) = tools[3];
    for swap in [false, true] {
        let cut = if swap { pin_volume } else { prism_volume };
        for (kind, n, derived) in [
            (BoolKind::Fuse, 2, prism_volume + pin_volume),
            (BoolKind::Cut, 1, cut),
            (BoolKind::Common, 0, 0.0),
        ] {
            let at = format!("the notch prism and the pin — {kind:?} swap={swap}");
            let mut m = Model::new();
            let (a, b) = (notch(&mut m), pin(&mut m));
            let (x, y) = if swap { (b, a) } else { (a, b) };
            let out = crate::boolean(&mut m, kind, x, y).unwrap_or_else(|e| panic!("{at}: {e:?}"));
            let got = built_volume(&m, &out, n, &at);
            assert!(
                close(got, derived, 1e-9) || n == 0,
                "{at}: {got} vs {derived}"
            );
        }
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

/// **A box turned off the quarters beside the cut rod is refused by the gate's name** — the box
/// `[8, 10] × [−1, 1] × [1, 3]` turned 30° about `z`, clear of the rod. Its turned walls have no world
/// statement; one cuts the rod's infinite cylinder with its faces apart, so its class draws the
/// rulings on its own arrangement (`Silence::Sectioned`), and where such a ruling ends on the slanted
/// rim the point is a nameless plane's meet with a section, stated only by intervals, whose height
/// the ruling split cannot place — every cell refuses `CylinderGateUndecided`. The same box beside
/// the rod before its cut builds (the common part is empty), so the refusal is the slanted rim's.
#[test]
fn a_turned_box_beside_the_cut_rod_is_refused_by_the_gates_name() {
    let turned = |m: &mut Model| {
        let b = crate::fixtures::cuboid(m, [8.0, -1.0, 1.0], [10.0, 1.0, 3.0]);
        super::fixtures::xf(m, b, super::fixtures::rot_iso(Axis::Z, 30))
    };
    let mut m = Model::new();
    let (r, t) = (rod(&mut m), turned(&mut m));
    let out = crate::boolean(&mut m, BoolKind::Common, r, t).expect("the uncut rod builds");
    assert!(out.is_empty(), "nothing in common: {out:?}");
    for swap in [false, true] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let (mut m, s) = cut_rod(0.5);
            let t = turned(&mut m);
            let (x, y) = if swap { (t, s) } else { (s, t) };
            match crate::boolean(&mut m, kind, x, y) {
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderGateUndecided,
                    ..
                }) => {}
                other => panic!("{kind:?} swap={swap}: {other:?}"),
            }
        }
    }
}
