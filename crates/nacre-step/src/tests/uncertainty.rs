use super::*;
use crate::to_step_solids;
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::{Point3, Vector3};
use nacre_ops::{Edge2d, OpOutput, Operation, Ring2d, SketchFrame};
use num_bigint::BigInt;

const STAMP: &str = "2026-10-26T00:00:00Z";

fn turn(m: &mut Model, s: Handle<Solid>, axis: Axis, deg: i128) -> Handle<Solid> {
    let isometry = Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).expect("a whole-degree angle"),
    });
    moved(m, s, isometry)
}

fn moved(m: &mut Model, solid: Handle<Solid>, isometry: Isometry) -> Handle<Solid> {
    let Ok(OpOutput::Transform { solid }) =
        nacre_ops::apply(m, &Operation::Transform { solid, isometry })
    else {
        panic!("the motion applies")
    };
    m.rebuild_adjacency();
    solid
}

fn cuboid(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
    nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi))
}

fn cylinder(m: &mut Model, base: [f64; 3], r: f64, h: f64) -> Handle<Solid> {
    nacre_ops::fixtures::cylinder(
        m,
        Point3::from_array(base),
        Vector3::from_array([0.0, 0.0, 1.0]),
        r,
        h,
    )
    .solid
}

/// A prism of `ring` on the world `xy` plane, `h` along `+z`.
fn prism(m: &mut Model, ring: Ring2d, h: f64) -> Handle<Solid> {
    let profile = nacre_ops::from_paths(vec![ring])
        .expect("one profile")
        .remove(0);
    let frame = SketchFrame::world(m, Axis::Z);
    let Ok(OpOutput::Extrude { solid, .. }) = nacre_ops::apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: h,
        },
    ) else {
        panic!("the prism extrudes")
    };
    m.rebuild_adjacency();
    solid
}

/// The deep chain the refine door pays for: a box turned 7° about `z`, then moved `1/7` along
/// `x` 300 times.
fn deep_box(m: &mut Model) -> Handle<Solid> {
    let mut s = cuboid(m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    s = turn(m, s, Axis::Z, 7);
    for _ in 0..300 {
        let step = Isometry::translation([
            Rat::new(1, 7).expect("1/7"),
            Rat::from_int(0),
            Rat::from_int(0),
        ]);
        s = moved(m, s, step);
    }
    s
}

/// The populations the bound is held on, each as a model and the solids written from it.
fn fixtures() -> Vec<(&'static str, Model, Handle<Solid>)> {
    let r = Rat::from_int;
    let mut out = Vec::new();
    let mut build = |name: &'static str, f: &dyn Fn(&mut Model) -> Handle<Solid>| {
        let mut m = Model::new();
        let s = f(&mut m);
        m.rebuild_adjacency();
        out.push((name, m, s));
    };
    build("integer box", &|m| cuboid(m, [0.0; 3], [2.0, 3.0, 4.0]));
    build("turned box", &|m| {
        let s = cuboid(m, [0.1, 0.7, 0.3], [1.2, 1.5, 2.4]);
        turn(m, s, Axis::Y, 37)
    });
    build("cylinder", &|m| cylinder(m, [1.0, 0.5, 0.0], 0.75, 2.0));
    build("bored plate", &|m| {
        let plate = cuboid(m, [-3.0, -3.0, 0.0], [3.0, 3.0, 2.0]);
        let bore = cylinder(m, [0.3, 0.2, -1.0], 1.25, 4.0);
        let s = nacre_ops::boolean(m, nacre_ops::BoolKind::Cut, plate, bore).expect("a bore")[0];
        turn(m, s, Axis::X, 23)
    });
    build("refined deep box", &|m| {
        let s = deep_box(m);
        nacre_ops::refine_caches(m);
        s
    });
    build("turned slot", &|m| {
        let arc = |c: [i128; 2]| Edge2d::Arc {
            center: c.map(r),
            r2: r(4),
            ccw: true,
        };
        let ring = Ring2d::new(
            vec![[r(2), r(0)], [r(2), r(10)], [r(-2), r(10)], [r(-2), r(0)]],
            vec![Edge2d::Line, arc([0, 10]), Edge2d::Line, arc([0, 0])],
        )
        .expect("a slot");
        let s = prism(m, ring, 3.0);
        turn(m, s, Axis::X, 30)
    });
    build("leaf", &|m| {
        let arc = |c: [i128; 2]| Edge2d::Arc {
            center: c.map(r),
            r2: r(25),
            ccw: false,
        };
        let ring = Ring2d::new(
            vec![[r(0), r(0)], [r(5), r(5)]],
            vec![arc([5, 0]), arc([0, 5])],
        )
        .expect("a leaf");
        prism(m, ring, 2.0)
    });
    build("turned hollow box", &|m| {
        let a = cuboid(m, [0.0; 3], [4.0; 3]);
        let b = cuboid(m, [1.0; 3], [3.0; 3]);
        let void = m.reversed_shell(m.solid(b).outer);
        let outer = m.solid(a).outer;
        let hollow = m.push_solid(Solid {
            outer,
            cavities: vec![void],
        });
        m.restore_live(vec![hollow]);
        m.rebuild_adjacency();
        turn(m, hollow, Axis::X, 23)
    });
    // The cavity holds the larger gaps: an exact box around the unrefined deep box's shell,
    // reversed — a file whose bound comes from its cavity alone.
    build("deep cavity", &|m| {
        let deep = deep_box(m);
        let void = m.reversed_shell(m.solid(deep).outer);
        let outer = cuboid(m, [40.0, -2.0, -1.0], [50.0, 6.0, 5.0]);
        let outer = m.solid(outer).outer;
        let hollow = m.push_solid(Solid {
            outer,
            cavities: vec![void],
        });
        m.restore_live(vec![hollow]);
        hollow
    });
    out
}

/// A finite `f64` as an integer on the common grid `2⁻¹⁰⁷⁴` — every `f64` lands on it whole.
fn big(v: f64) -> BigInt {
    let bits = v.to_bits();
    let raw = ((bits >> 52) & 0x7ff) as i64;
    let frac = bits & ((1u64 << 52) - 1);
    let (mantissa, exp) = match raw {
        0 => (frac, -1074),
        _ => (frac | (1u64 << 52), raw - 1075),
    };
    let m = BigInt::from(mantissa) << (exp + 1074) as usize;
    if v.is_sign_negative() { -m } else { m }
}

type B3 = [BigInt; 3];

fn b3(p: [f64; 3]) -> B3 {
    p.map(big)
}

fn bsub(a: &B3, b: &B3) -> B3 {
    core::array::from_fn(|k| &a[k] - &b[k])
}

fn bdot(a: &B3, b: &B3) -> BigInt {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn bcross(a: &B3, b: &B3) -> B3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

/// `(r − u)²·aa ≤ ww ≤ (r + u)²·aa` — a point whose offset from an axis, crossed with the axis
/// direction `a` (`aa = a·a`), has square norm `ww` stands within `u` of radius `r`.
fn within_round(ww: &BigInt, aa: &BigInt, r: &BigInt, u: &BigInt) -> bool {
    let hi = (r + u) * (r + u) * aa;
    let lo = if r > u {
        (r - u) * (r - u) * aa
    } else {
        BigInt::from(0)
    };
    lo <= *ww && *ww <= hi
}

/// **The written vertices stand within `u` of every curve and surface they lie on, and a line's
/// far end within `u` of its faces' surfaces**, checked in
/// exact integers on the numbers the file holds — squared, so no root is taken. A vertex off a
/// plane by `h·|n|` breaks `(n·(v−o))² ≤ u²(n·n)`; off a round surface, the squared radius of its
/// offset leaves `[(r−u)², (r+u)²]`; off a line, `|(v−p)×d|² ≤ u²|d|²`; off a circle, its height
/// over the circle's plane and its radius each. Returns what broke.
fn vertex_gaps_beyond(m: &Model, s: Handle<Solid>, u: f64) -> Vec<String> {
    let ub = big(u);
    let mut bad = Vec::new();
    let solid = m.solid(s);
    for shell in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
        for &fh in &m.shell(shell).faces {
            let face = m.face(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    let [v0, v1] = m.edge(he.edge).vertices;
                    let p = b3(m.vertex_point(v0).as_array());
                    // The far end of a line's stretch: the foot of `v1`, `p + t·d` with
                    // `t = (v1 − p)·d / d·d` — kept whole as `f·dd` (`dd = d·d`), so the
                    // edge's own gap from the face's surface is checked where it peaks for a
                    // plane, and at the end for a round face.
                    if let nacre_geom::Curve::Line(l) = m.edge_curve(he.edge) {
                        let d = b3(l.direction().as_array());
                        let dd = bdot(&d, &d);
                        let t = bdot(&bsub(&b3(m.vertex_point(v1).as_array()), &p), &d);
                        let foot: B3 = core::array::from_fn(|k| &p[k] * &dd + &t * &d[k]);
                        let scaled = |x: &B3| -> B3 { core::array::from_fn(|k| &x[k] * &dd) };
                        match m.surface_cache(face.surface) {
                            SurfaceCache::Plane { plane, .. } => {
                                let (o, n) =
                                    (b3(plane.origin().as_array()), b3(plane.normal().as_array()));
                                let h = bdot(&n, &bsub(&foot, &scaled(&o)));
                                if &h * &h > &ub * &ub * bdot(&n, &n) * &dd * &dd {
                                    bad.push(format!(
                                        "e{}'s far end off its plane",
                                        he.edge.index()
                                    ));
                                }
                            }
                            SurfaceCache::Cylinder { cylinder, .. } => {
                                let (a0, a) = (
                                    b3(cylinder.axis().origin().as_array()),
                                    b3(cylinder.axis().direction().as_array()),
                                );
                                let w = bcross(&bsub(&foot, &scaled(&a0)), &a);
                                let aa = bdot(&a, &a) * &dd * &dd;
                                if !within_round(&bdot(&w, &w), &aa, &big(cylinder.radius()), &ub) {
                                    bad.push(format!(
                                        "e{}'s far end off its cylinder",
                                        he.edge.index()
                                    ));
                                }
                            }
                        }
                    }
                    for vh in [v0, v1] {
                        let v = b3(m.vertex_point(vh).as_array());
                        match m.surface_cache(face.surface) {
                            SurfaceCache::Plane { plane, .. } => {
                                let (o, n) =
                                    (b3(plane.origin().as_array()), b3(plane.normal().as_array()));
                                let h = bdot(&n, &bsub(&v, &o));
                                if &h * &h > &ub * &ub * bdot(&n, &n) {
                                    bad.push(format!("v{} off its plane", vh.index()));
                                }
                            }
                            SurfaceCache::Cylinder { cylinder, .. } => {
                                let (a0, a) = (
                                    b3(cylinder.axis().origin().as_array()),
                                    b3(cylinder.axis().direction().as_array()),
                                );
                                let w = bcross(&bsub(&v, &a0), &a);
                                let r = big(cylinder.radius());
                                if !within_round(&bdot(&w, &w), &bdot(&a, &a), &r, &ub) {
                                    bad.push(format!("v{} off its cylinder", vh.index()));
                                }
                            }
                        }
                        match m.edge_curve(he.edge) {
                            nacre_geom::Curve::Line(l) => {
                                let d = b3(l.direction().as_array());
                                let x = bcross(&bsub(&v, &p), &d);
                                if bdot(&x, &x) > &ub * &ub * bdot(&d, &d) {
                                    bad.push(format!("v{} off its line", vh.index()));
                                }
                            }
                            nacre_geom::Curve::Circle(c) => {
                                let (cc, n) =
                                    (b3(c.center().as_array()), b3(c.normal().as_array()));
                                let d = bsub(&v, &cc);
                                let h = bdot(&n, &d);
                                let nn = bdot(&n, &n);
                                if &h * &h > &ub * &ub * &nn {
                                    bad.push(format!("v{} off its circle's plane", vh.index()));
                                }
                                let y = bcross(&d, &n);
                                if !within_round(&bdot(&y, &y), &nn, &big(c.radius()), &ub) {
                                    bad.push(format!("v{} off its circle's radius", vh.index()));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    bad
}

/// The `u` a fixture's file states, as a positive proven distance.
fn proven(m: &Model, s: Handle<Solid>) -> f64 {
    match uncertainty(m, &[s]) {
        Uncertainty::Proven(u) => {
            assert!(u > 0.0 && u.is_finite(), "{u}");
            u
        }
        other => panic!("{other:?}"),
    }
}

/// **The stated uncertainty bounds every gap among the written numbers**, on populations that
/// reach every kind of meeting the writer emits: a plane and a round face, straight and circular
/// edges, turned (the numbers rounded) and not, a pierce corner on a turned slot, two arcs of
/// different circles, a deep chain the refine door paid for, and a cavity shell. The check is in
/// exact integers on the written `f64`s and shares no arithmetic with the bound.
#[test]
fn the_uncertainty_bounds_every_gap_the_file_holds() {
    let mut bad = Vec::new();
    for (name, m, s) in fixtures() {
        let u = proven(&m, s);
        for b in vertex_gaps_beyond(&m, s, u) {
            bad.push(format!("{name}: {b}"));
        }
    }
    assert!(bad.is_empty(), "{bad:#?}");
}

/// **The file states that distance**, read back: step-io's `precision` is the uncertainty's own
/// bits.
#[test]
fn the_file_states_the_uncertainty_it_proves() {
    for (name, m, s) in fixtures() {
        let u = proven(&m, s);
        let text = to_step_solids(&m, &[s], STAMP).expect("export");
        let (model, _) = step_io::read(text.as_bytes()).expect("re-read");
        let read = model.scene().units().precision;
        assert_eq!(
            read.map(f64::to_bits),
            Some(u.to_bits()),
            "{name}: {read:?} vs {u:e}"
        );
    }
}

/// **A model the refine door has not paid for is measured too** — its figures stand where their
/// construction put them, and the file says how far apart they are. The deep box before the door
/// is proven and bounded like any other, and here the door, settling every figure to the nearest
/// `f64`, tightens it.
#[test]
fn an_unrefined_model_states_what_its_figures_hold() {
    let mut m = Model::new();
    let s = deep_box(&mut m);
    m.rebuild_adjacency();
    let before = proven(&m, s);
    let bad = vertex_gaps_beyond(&m, s, before);
    assert!(bad.is_empty(), "{bad:#?}");
    nacre_ops::refine_caches(&mut m);
    let after = proven(&m, s);
    assert!(after < before, "refined {after:e}, figures {before:e}");
}

/// **Where every gap is exactly zero, the file states its coordinate spacing** — half the `f64`
/// spacing at its largest coordinate, positive as the standard requires. An integer box's
/// corners lie exactly on its planes and lines; its largest coordinate is 4.
#[test]
fn an_exact_model_states_its_coordinate_spacing() {
    let mut m = Model::new();
    let s = cuboid(&mut m, [0.0; 3], [2.0, 3.0, 4.0]);
    m.rebuild_adjacency();
    assert_eq!(measure(&m, &[s]).expect("measured").max, 0.0);
    assert_eq!(
        proven(&m, s),
        (f64::from_bits(4f64.to_bits() + 1) - 4.0) / 2.0
    );
}

/// **Nothing to write proves nothing**, and the file states the writer's convention.
#[test]
fn an_empty_export_states_the_writers_convention() {
    let m = Model::new();
    assert_eq!(
        uncertainty(&m, &[]),
        Uncertainty::NotProven(NotProven::NothingExported)
    );
    let text = to_step_solids(&m, &[], STAMP).expect("an empty file");
    let (model, _) = step_io::read(text.as_bytes()).expect("re-read");
    assert_eq!(model.scene().units().precision, Some(1e-7));
}
