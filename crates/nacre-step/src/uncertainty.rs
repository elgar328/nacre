//! The uncertainty a file states: how far apart the numbers it writes may stand where they meet.

use nacre_exact::HpBounded;
use nacre_geom::Curve;
use nacre_store::Handle;
use nacre_topo::{Edge, Model, Shell, Solid, Surface, SurfaceCache, Vertex};
use std::collections::HashMap;

/// What [`uncertainty`] found for a set of solids.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum Uncertainty {
    /// A distance, in the file's length unit, that no gap among the written numbers exceeds:
    /// between a vertex and each curve and surface it lies on, and between an edge's curve — over
    /// the stretch its two vertices bound, a whole circle for an arc — and the surfaces of its
    /// faces. Never below half the spacing of `f64`s at the largest coordinate written — the
    /// resolution the file's coordinates are written at — so always positive.
    Proven(f64),
    /// No such distance was proven, and the file states the writer's convention instead.
    NotProven(NotProven),
}

/// Why [`uncertainty`] proved no distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum NotProven {
    /// There is no solid to write.
    NothingExported,
    /// A gap could not be bounded: a written direction or radius too close to zero to divide
    /// by, or a number past the range of `f64`.
    Undecided,
}

/// **The uncertainty a STEP file of `solids` states** — computed on the very numbers the file
/// writes (the model's caches as they stand), so it says how well the file holds together, not
/// how close it stands to the exact geometry. How close it stands is what
/// `nacre_ops::refine_caches` (`nacre::ops::refine_caches`) settles and reports before an export.
///
/// [`to_step`](crate::to_step), [`to_step_solids`](crate::to_step_solids) and
/// [`write_step`](crate::write_step) write [`Uncertainty::Proven`]'s distance as the file's
/// `UNCERTAINTY_MEASURE_WITH_UNIT`, and the writer's convention — `1e-7` mm — where it is
/// [`Uncertainty::NotProven`].
pub fn uncertainty(model: &Model, solids: &[Handle<Solid>]) -> Uncertainty {
    if solids.is_empty() {
        return Uncertainty::NotProven(NotProven::NothingExported);
    }
    match measure(model, solids) {
        Some(gaps) => Uncertainty::Proven(gaps.max.max(half_ulp(gaps.largest))),
        None => Uncertainty::NotProven(NotProven::Undecided),
    }
}

/// Every gap of `solids` bounded — `None` where one cannot be.
pub(crate) fn measure(model: &Model, solids: &[Handle<Solid>]) -> Option<Gaps> {
    let mut gaps = Gaps::default();
    let mut met = Met::default();
    for &s in solids {
        let solid = model.solid(s);
        for shell in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
            gaps.shell(model, shell, &mut met)?;
        }
    }
    Some(gaps)
}

/// The bits every distance is computed at: the written `f64`s enter exactly, and what rounding
/// the products and quotients add is carried in each value's radius.
const PREC: usize = 256;

type V3 = [HpBounded; 3];

/// The distances already bounded, so a pair met again is not measured again: a vertex from a
/// face's surface (two edges of a loop meet at each vertex), and an edge's second vertex from its
/// curve (each edge bounds two faces).
#[derive(Default)]
struct Met {
    off_surface: HashMap<(Handle<Vertex>, Handle<Surface>), f64>,
    off_curve: HashMap<Handle<Edge>, f64>,
}

/// The largest gap bounded so far and the largest coordinate written.
#[derive(Default, Debug)]
pub(crate) struct Gaps {
    pub(crate) max: f64,
    pub(crate) largest: f64,
}

impl Gaps {
    fn gap(&mut self, d: f64) {
        self.max = self.max.max(d);
    }

    fn point(&mut self, p: [f64; 3]) -> Option<V3> {
        for c in p {
            self.largest = self.largest.max(c.abs());
        }
        v3(p)
    }

    /// Every face of `shell`: each edge of its loops and the edge's two vertices against the
    /// face's surface, and the two vertices against the edge's curve.
    fn shell(&mut self, model: &Model, shell: Handle<Shell>, met: &mut Met) -> Option<()> {
        for &fh in &model.shell(shell).faces {
            let face = model.face(fh);
            let surface = Carrier::of(self, model.surface_cache(face.surface))?;
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    let ends = model.edge(he.edge).vertices;
                    let points = ends.map(|v| model.vertex_point(v).as_array());
                    let [v0, v1] = [self.point(points[0])?, self.point(points[1])?];
                    let curve = Path::of(self, model.edge_curve(he.edge))?;
                    let far = match met.off_curve.get(&he.edge) {
                        Some(&d) => d,
                        None => {
                            self.gap(curve.to_from(&v0, &v0)?);
                            let d = curve.to_from(&v1, &v0)?;
                            self.gap(d);
                            met.off_curve.insert(he.edge, d);
                            d
                        }
                    };
                    let mut off = [0.0; 2];
                    for (k, (vh, v)) in ends.into_iter().zip([&v0, &v1]).enumerate() {
                        off[k] = match met.off_surface.get(&(vh, face.surface)) {
                            Some(&d) => d,
                            None => {
                                let d = surface.to(v)?;
                                self.gap(d);
                                met.off_surface.insert((vh, face.surface), d);
                                d
                            }
                        };
                    }
                    self.gap(curve.against(&surface, &v0, &v1, off, far)?);
                }
            }
        }
        Some(())
    }
}

/// A written surface.
enum Carrier {
    /// Through `origin`, normal `normal` (not taken as a unit).
    Plane { origin: V3, normal: V3 },
    /// About the axis through `origin` along `axis` (not taken as a unit), of `radius`.
    Cylinder {
        origin: V3,
        axis: V3,
        radius: HpBounded,
    },
}

impl Carrier {
    fn of(gaps: &mut Gaps, cache: &SurfaceCache) -> Option<Carrier> {
        Some(match cache {
            SurfaceCache::Plane { plane, .. } => Carrier::Plane {
                origin: gaps.point(plane.origin().as_array())?,
                normal: v3(plane.normal().as_array())?,
            },
            SurfaceCache::Cylinder { cylinder, .. } => Carrier::Cylinder {
                origin: gaps.point(cylinder.axis().origin().as_array())?,
                axis: v3(cylinder.axis().direction().as_array())?,
                radius: hp(cylinder.radius())?,
            },
        })
    }

    /// How far `v` stands from the surface, at most.
    fn to(&self, v: &V3) -> Option<f64> {
        match self {
            Carrier::Plane { origin, normal } => off_plane(v, origin, normal),
            Carrier::Cylinder {
                origin,
                axis,
                radius,
            } => off_round(&cross(&sub(v, origin), axis), axis, radius),
        }
    }
}

/// A written edge curve. One lives on the stack for the edge it is read from, so the line's
/// smaller variant costs nothing worth boxing the circle for.
#[allow(clippy::large_enum_variant)]
enum Path {
    /// Through the edge's first vertex — the point the file gives the line — along `dir`.
    Line { dir: V3 },
    /// About `centre` in the plane normal to `normal` (not taken as a unit), of `radius`.
    Circle {
        centre: V3,
        normal: V3,
        radius: HpBounded,
    },
    /// The file's `ELLIPSE`: `centre + a cos u·x + b sin u·y`, `x` the written reference direction
    /// made unit and normal to the written normal, `y` that unit normal crossed with `x` — the
    /// placement as a reader builds it from the written numbers.
    Ellipse {
        centre: V3,
        x: V3,
        y: V3,
        a: HpBounded,
        b: HpBounded,
    },
}

impl Path {
    fn of(gaps: &mut Gaps, curve: &Curve) -> Option<Path> {
        Some(match curve {
            Curve::Line(l) => Path::Line {
                dir: v3(l.direction().as_array())?,
            },
            // The numbers `build_edge` writes, read from the same `written`: a circle's centre,
            // axis and radius, or an ellipse's placement and semi-axes.
            Curve::Ellipse(c) if c.is_circle() => {
                let w = c.written();
                Path::Circle {
                    centre: gaps.point(w.center.as_array())?,
                    normal: v3(w.normal.as_array())?,
                    radius: hp(w.semi_1)?,
                }
            }
            Curve::Ellipse(c) => {
                let w = c.written();
                let n = unit(&v3(w.normal.as_array())?)?;
                let r = v3(w.ref_dir.as_array())?;
                let along = dot(&r, &n);
                let x = unit(&core::array::from_fn(|k| {
                    r[k].sub(&along.mul(&n[k], PREC), PREC)
                }))?;
                Path::Ellipse {
                    centre: gaps.point(w.center.as_array())?,
                    y: cross(&n, &x),
                    x,
                    a: hp(w.semi_1)?,
                    b: hp(w.semi_2)?,
                }
            }
        })
    }

    /// How far `v` stands from the curve, at most — `start` being the line's point.
    fn to_from(&self, v: &V3, start: &V3) -> Option<f64> {
        match self {
            Path::Line { dir } => {
                let x = cross(&sub(v, start), dir);
                root_up(&ratio(&dot(&x, &x), &dot(dir, dir))?)
            }
            Path::Circle {
                centre,
                normal,
                radius,
            } => {
                let d = sub(v, centre);
                let axial = off_plane(v, centre, normal)?;
                let radial = off_round(&cross(&d, normal), normal, radius)?;
                Some(sum_up(axial, radial))
            }
            // No closed form reaches the nearest point (a quartic), so the bound is the distance to
            // **a** point of the curve: the angle `v` itself points to in the ellipse's own
            // coordinates, `(cos u, sin u) ∝ ((v − c)·x / a, (v − c)·y / b)` — normalized in the
            // intervals, whose true values lie on the unit circle, so the true point lies on the
            // curve and inside the interval point.
            Path::Ellipse { centre, x, y, a, b } => {
                let d = sub(v, centre);
                let (p, q) = (ratio(&dot(&d, x), a)?, ratio(&dot(&d, y), b)?);
                let k = p.mul(&p, PREC).add(&q.mul(&q, PREC), PREC).inv_sqrt(PREC)?;
                let (cu, su) = (p.mul(&k, PREC), q.mul(&k, PREC));
                let (ac, bs) = (a.mul(&cu, PREC), b.mul(&su, PREC));
                let gap: V3 = core::array::from_fn(|i| {
                    d[i].sub(&ac.mul(&x[i], PREC), PREC)
                        .sub(&bs.mul(&y[i], PREC), PREC)
                });
                root_up(&dot(&gap, &gap))
            }
        }
    }

    /// How far the curve's stretch between `v0` and `v1` stands from `surface`, at most — `off`
    /// being the two vertices' distances from the surface and `far` `v1`'s from the curve.
    fn against(&self, surface: &Carrier, v0: &V3, v1: &V3, off: [f64; 2], far: f64) -> Option<f64> {
        match (self, surface) {
            // The distance to a plane is affine along a line, so it peaks at an end of the
            // stretch: at the line's point (`v0` itself) or at the foot of `v1`, which stands no
            // further from the plane than `v1` plus `v1`'s distance from the line.
            (Path::Line { .. }, Carrier::Plane { .. }) => Some(off[0].max(sum_up(off[1], far))),
            // The distance from the axis changes along the line by at most the length walked
            // times the sine of the angle between the two — the length walked being the foot of
            // `v1`'s parameter, `(v1 − v0)·d / |d|`.
            (Path::Line { dir }, Carrier::Cylinder { axis, .. }) => {
                let along = dot(&sub(v1, v0), dir);
                let turn = cross(dir, axis);
                let num = along.mul(&along, PREC).mul(&dot(&turn, &turn), PREC);
                let dd = dot(dir, dir);
                let den = dd.mul(&dd, PREC).mul(&dot(axis, axis), PREC);
                Some(sum_up(off[0], root_up(&ratio(&num, &den)?)?))
            }
            // The centre's height over the plane, plus the radius times the sine of the tilt
            // between the two normals.
            (
                Path::Circle {
                    centre,
                    normal: m,
                    radius,
                },
                Carrier::Plane { origin, normal: n },
            ) => {
                let height = off_plane(centre, origin, n)?;
                let tilt = cross(n, m);
                let sine = root_up(&ratio(
                    &dot(&tilt, &tilt),
                    &dot(n, n).mul(&dot(m, m), PREC),
                )?)?;
                Some(sum_up(height, mul_up(upper(radius)?, sine)))
            }
            // A point of the circle stands from the axis within `[r_c|cos φ| − e, r_c + e]`, `e`
            // the centre's distance from the axis and `φ` the tilt between the circle's normal
            // and the axis; `1 − |cos φ| ≤ sin² φ`, so the gap is at most
            // `e + |r_c − r| + r_c·sin² φ`.
            (
                Path::Circle {
                    centre,
                    normal: m,
                    radius: rc,
                },
                Carrier::Cylinder {
                    origin,
                    axis,
                    radius,
                },
            ) => {
                let off = cross(&sub(centre, origin), axis);
                let e = root_up(&ratio(&dot(&off, &off), &dot(axis, axis))?)?;
                let tilt = cross(m, axis);
                let sine2 = upper(&ratio(
                    &dot(&tilt, &tilt),
                    &dot(m, m).mul(&dot(axis, axis), PREC),
                )?)?;
                let radii = upper(&rc.sub(radius, PREC))?;
                Some(sum_up(sum_up(e, radii), mul_up(upper(rc)?, sine2)))
            }
            // The whole ellipse: the centre's height over the plane, plus the most
            // `n·(a cos u·x + b sin u·y)` reaches, `√((a n·x)² + (b n·y)²)`, over `|n|`.
            (Path::Ellipse { centre, x, y, a, b }, Carrier::Plane { origin, normal: n }) => {
                let height = off_plane(centre, origin, n)?;
                let (p, q) = (a.mul(&dot(n, x), PREC), b.mul(&dot(n, y), PREC));
                let reach = root_up(&ratio(
                    &p.mul(&p, PREC).add(&q.mul(&q, PREC), PREC),
                    &dot(n, n),
                )?)?;
                Some(sum_up(height, reach))
            }
            // A point of the ellipse stands from the axis at `|W|/|A|`, `W = D + cos u·P + sin u·Q`
            // with `D`, `P`, `Q` the centre's offset and the two semi-diameters crossed with the
            // axis `A`. `| |W|/|A| − R | ≤ | |W|² − R²|A|² | / (|A|² R)`, and
            // `|W|² = D·D + (P·P + Q·Q)/2 + cos 2u·(P·P − Q·Q)/2 + sin 2u·P·Q + 2(cos u·D·P + sin u·D·Q)`,
            // so the gap is at most the constant part's size, plus the two amplitudes
            // `√(((P·P − Q·Q)/2)² + (P·Q)²)` and `2√((D·P)² + (D·Q)²)`, over `|A|² R`.
            (
                Path::Ellipse { centre, x, y, a, b },
                Carrier::Cylinder {
                    origin,
                    axis,
                    radius,
                },
            ) => {
                let scaled =
                    |s: &HpBounded, v: &V3| -> V3 { core::array::from_fn(|i| s.mul(&v[i], PREC)) };
                let dd = cross(&sub(centre, origin), axis);
                let pp = cross(&scaled(a, x), axis);
                let qq = cross(&scaled(b, y), axis);
                let aa = dot(axis, axis);
                let den = aa.mul(radius, PREC);
                let (p2, q2, pq) = (dot(&pp, &pp), dot(&qq, &qq), dot(&pp, &qq));
                let (dp, dq) = (dot(&dd, &pp), dot(&dd, &qq));
                let half = hp(0.5)?;
                let constant = dot(&dd, &dd)
                    .add(&p2.add(&q2, PREC).mul(&half, PREC), PREC)
                    .sub(&radius.mul(radius, PREC).mul(&aa, PREC), PREC);
                let swing = p2.sub(&q2, PREC).mul(&half, PREC);
                let den2 = den.mul(&den, PREC);
                let level = upper(&ratio(&constant, &den)?)?;
                let wobble = root_up(&ratio(
                    &swing.mul(&swing, PREC).add(&pq.mul(&pq, PREC), PREC),
                    &den2,
                )?)?;
                let shift = root_up(&ratio(
                    &dp.mul(&dp, PREC).add(&dq.mul(&dq, PREC), PREC),
                    &den2,
                )?)?;
                Some(sum_up(sum_up(level, wobble), mul_up(2.0, shift)))
            }
        }
    }
}

/// `|n·(v − o)| / |n|`, at most.
fn off_plane(v: &V3, o: &V3, n: &V3) -> Option<f64> {
    let h = dot(&sub(v, o), n);
    root_up(&ratio(&h.mul(&h, PREC), &dot(n, n))?)
}

/// `| |w| / |a| − r |` at most — the distance from a round surface of radius `r` of a point whose
/// offset from the axis, crossed with the axis direction `a`, is `w`. Without a root:
/// `| |w|²/|a|² − r² | / (|w|/|a| + r) ≤ | |w|² − r²|a|² | / (|a|² r)`.
fn off_round(w: &V3, a: &V3, r: &HpBounded) -> Option<f64> {
    let aa = dot(a, a);
    let excess = dot(w, w).sub(&r.mul(r, PREC).mul(&aa, PREC), PREC);
    upper(&ratio(&excess, &aa.mul(r, PREC))?)
}

/// `v / |v|`, `None` where `|v|` may be zero.
fn unit(v: &V3) -> Option<V3> {
    let k = dot(v, v).inv_sqrt(PREC)?;
    Some(core::array::from_fn(|i| v[i].mul(&k, PREC)))
}

fn hp(v: f64) -> Option<HpBounded> {
    HpBounded::of_f64(v, PREC)
}

fn v3(p: [f64; 3]) -> Option<V3> {
    Some([hp(p[0])?, hp(p[1])?, hp(p[2])?])
}

fn sub(a: &V3, b: &V3) -> V3 {
    core::array::from_fn(|k| a[k].sub(&b[k], PREC))
}

fn dot(a: &V3, b: &V3) -> HpBounded {
    a[0].mul(&b[0], PREC)
        .add(&a[1].mul(&b[1], PREC), PREC)
        .add(&a[2].mul(&b[2], PREC), PREC)
}

fn cross(a: &V3, b: &V3) -> V3 {
    let x = |i: usize, j: usize| a[i].mul(&b[j], PREC).sub(&a[j].mul(&b[i], PREC), PREC);
    [x(1, 2), x(2, 0), x(0, 1)]
}

/// `num / den`, `None` where `den` may be zero.
fn ratio(num: &HpBounded, den: &HpBounded) -> Option<HpBounded> {
    num.div(den, PREC)
}

/// An `f64` at or above every magnitude `x` holds.
fn upper(x: &HpBounded) -> Option<f64> {
    x.abs_upper_f64()
}

/// An `f64` at or above `√|x|`: the square root of the upper bound, one step up.
fn root_up(x: &HpBounded) -> Option<f64> {
    Some(step_up(upper(x)?.sqrt()))
}

fn sum_up(a: f64, b: f64) -> f64 {
    step_up(a + b)
}

fn mul_up(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return 0.0;
    }
    step_up(a * b).max(f64::from_bits(1))
}

/// The next `f64` above a non-negative one — so a rounded-to-nearest result still bounds the
/// exact one. Zero stays zero: every caller's zero is exact (a sum of zeros, a root of zero).
fn step_up(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        f64::from_bits(x.to_bits() + 1)
    }
}

/// Half the spacing of `f64`s at `x` — at least the smallest normal `f64`.
fn half_ulp(x: f64) -> f64 {
    ((f64::from_bits(x.to_bits() + 1) - x) / 2.0).max(f64::MIN_POSITIVE)
}

#[cfg(test)]
#[path = "tests/uncertainty.rs"]
mod tests;
