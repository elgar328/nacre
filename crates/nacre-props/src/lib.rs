//! Exact derived properties of a nacre solid — volume, surface area, centroid,
//! bounding box, and the per-face facts ([`FaceProps`]) a caller needs to *name*
//! a face from its geometry rather than by an index that shifts under it.
//!
//! Computed **analytically on the exact geometry** via the divergence theorem —
//! not from a tessellation, so the numbers match an analytic oracle (OCCT) to
//! near machine precision rather than to mesh resolution. This crate is a
//! read-only analysis over the topology [`Model`]; it lives above `nacre-topo`
//! (not inside it) so the topology crate stays a truth-only aggregate, mirroring
//! `nacre-validate` / `nacre-tess` / `nacre-step`.
//!
//! Consumers: user "part volume / surface area" queries, the boolean
//! volume-conservation invariants (proptest), and the OCCT oracle diff
//! (`nacre-oracle`).
//!
//! Precondition: the solid is a valid closed, outward-oriented b-rep (what the
//! producers emit and `validate` accepts). Coordinates are the cache side of the
//! truth/cache split, so f64 arithmetic here is appropriate.

#![cfg_attr(not(test), deny(clippy::print_stdout, clippy::print_stderr))]
use nacre_geom::{Curve, Surface};
use nacre_math::{Point3, Vector3};
use nacre_store::Handle;
use nacre_topo::{Face, Loop, Model, Solid};

use std::f64::consts::PI;

/// Volume and surface area of one solid. The centroid is [`centroid`], a
/// separate function: it needs a planar boundary, and folding it in here would
/// take volume and area down with it on a curved solid. (Inertia is still
/// deferred — no consumer has asked.)
///
/// No `PartialEq`: the fields are `f64`, and production geometry never compares
/// coordinates with a bare `==` (tests use a tolerance).
#[derive(Clone, Copy, Debug)]
pub struct MassProps {
    pub volume: f64,
    pub area: f64,
}

/// A shape this analysis does not yet cover.
///
/// Named for the crate, not for mass: [`bounds`] and [`face_props`] return it too,
/// and an error type called `MassError` coming back from a bounding-box query
/// would be a name that lies.
///
/// There is deliberately **no** `UnsupportedSurface`: the surface `match` is
/// exhaustive over [`Surface`], so a new variant (NURBS, sphere…) is a compile
/// error here until it is handled — the same "new variant forces handling"
/// idiom `Surface::distance` and `to_step` use.
#[derive(Debug)]
pub enum PropsError {
    /// A planar face bounded by something other than a straight polygon or a
    /// single full circle (e.g. a future line+arc mix). None occur today.
    UnsupportedBoundary,
    /// Not returned: a curved face's inner loops — holes, and a lateral's second rim — are
    /// integrated like every other boundary (`lateral_moments`).
    UnsupportedInnerLoop,
    /// [`centroid`] met a **curved** face. Volume and area handle those
    /// analytically, but the cone argument the centroid rests on needs a planar
    /// base, so it is refused rather than approximated. Volume and area on the
    /// same solid are unaffected — which is why the centroid is its own
    /// function and not a [`MassProps`] field.
    CentroidOfCurvedFace,
}

/// Exact mass properties of one solid via the divergence theorem.
///
/// `V = (1/3) ∮ (r − R)·n̂ dA`, area `= Σ face area`, summed face by face. The
/// local reference point `R` (a vertex of the solid) is subtracted before the
/// flux dot product: the closed-surface identity `∮ n̂ dA = 0` makes `V`
/// independent of `R`, while choosing `R` near the solid removes the
/// catastrophic cancellation a far-from-origin placement would otherwise cause.
pub fn mass_props(model: &Model, solid: Handle<Solid>) -> Result<MassProps, PropsError> {
    let solid = model.solid(solid);
    let outer = model.shell(solid.outer);

    // R = the first vertex of the first face of the outer shell; any vertex on
    // the solid works. The closed-surface identity ∮ n̂ dA = 0 holds over the
    // *full* boundary — outer shell plus every cavity shell — so V is
    // R-independent; keeping R on the outer shell keeps the numbers small.
    let first_face = model.face(outer.faces[0]);
    let reference = model.vertex_point(he_start(model, first_face.outer.half_edges[0])?);

    let mut volume_flux = 0.0;
    let mut area = 0.0;
    // Outer boundary, then each inner cavity shell. A cavity's faces are wound
    // with their outward normals pointing into the void
    // (containment), so its flux is negative and subtracts the void's volume;
    // its (unsigned) area adds — both surfaces bound material.
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &face in &model.shell(sh).faces {
            let (a, flux) = face_contribution(model, model.face(face), reference)?;
            area += a;
            volume_flux += flux;
        }
    }
    Ok(MassProps {
        volume: volume_flux / 3.0,
        area,
    })
}

/// Area, centroid and outward normal of one face — the facts a caller filters on
/// when it has to *name* a face ("the one facing +Z, highest up") rather than
/// click it. Re-deciding it every run from geometry is what keeps a script from
/// depending on face indices that shift when the model upstream changes.
#[derive(Clone, Copy, Debug)]
pub struct FaceProps {
    /// The face's area, **holes subtracted**.
    pub area: f64,
    /// The area-weighted centroid of the face's region, holes subtracted.
    pub centroid: Point3,
    /// The outward normal — `None` for a **curved** face, which has no single
    /// one. A filter over every face of a solid can then skip those instead of
    /// failing, which is what "pick the face facing +Z" wants.
    pub normal: Option<Vector3>,
}

/// **The face's outward normal at one point of it** — [`FaceProps::normal`]'s answer for the
/// faces that have no single one.
///
/// A planar face faces the same way everywhere, so `p` changes nothing there; a cylindrical one
/// faces radially, so it changes everything. Both come back as the **outward** normal — the
/// surface's natural direction turned by the face's `orientation` — which is the sense a caller
/// means when it asks which way a face looks:
///
/// * a boss's wall points **away** from its axis, and
/// * a bore's wall points **toward** it, because the material is outside the wall.
///
/// That flip is the kernel's convention (the boolean's `flip` is exactly it), so it lives here
/// rather than in each caller that draws or measures a face.
///
/// `p` is taken on trust: it should lie on the face's surface, and a mesh corner or a sampled
/// point does. `None` only where the direction genuinely has no name — a point on a cylinder's
/// own axis.
pub fn face_normal_at(model: &Model, face: Handle<Face>, p: Point3) -> Option<Vector3> {
    let f = model.face(face);
    let sign = f64::from(f.orientation.sign());
    Some(model.surface_cache(f.surface).normal_at(p)? * sign)
}

/// [`FaceProps`] of one face.
pub fn face_props(model: &Model, face: Handle<Face>) -> Result<FaceProps, PropsError> {
    let face = model.face(face);
    let sign = f64::from(face.orientation.sign());
    // ★ The cache, because every answer below is an `f64` area or centroid: the integral is
    // evaluated numerically, and an exact normal would be rounded into the same product a line
    // later. Kind questions go to the truth; this one's answer is a number.
    match model.surface_cache(face.surface) {
        Surface::Plane(plane) => {
            let (area, centroid) = planar_region(model, face)?;
            Ok(FaceProps {
                area,
                centroid,
                normal: Some(plane.normal() * sign),
            })
        }
        Surface::Cylinder(cyl) => {
            let m = lateral_moments(model, face, cyl)?;
            // The mean of `p = a₀ + z·axis + r·n̂(θ)` over the region: the axis point at the
            // mean axial station, plus `r·∬n̂ / ∬1` radially (exactly zero for a full band, whose
            // mean is on the axis). Every `sign` cancels in a ratio of two moments.
            Ok(FaceProps {
                area: cyl.radius() * sign * m.j1,
                centroid: cyl.axis().origin()
                    + cyl.axis().direction() * (m.jz / m.j1 + m.z0)
                    + m.jn * (cyl.radius() / m.j1),
                normal: None,
            })
        }
    }
}

/// The centroid of one solid.
///
/// Reuses the divergence-theorem decomposition [`mass_props`] already performs:
/// the solid is the signed sum of **cones** from a reference point `R` over each
/// planar face, and a cone's centroid sits three quarters of the way from its
/// apex to the base's centroid — for *any* base outline, because the
/// cross-sections are scaled copies. So `C = R + Σ Vᵢ·¾(cᵢ − R) / Σ Vᵢ` with
/// `Vᵢ = Aᵢ·n̂ᵢ·(cᵢ − R)/3`, all of which [`face_props`] already supplies. Holes
/// need no special case: `Aᵢ` and `cᵢ` are the region's, holes already removed.
///
/// A **curved** face breaks that argument (the cone over it is not one), so a
/// solid carrying any is refused rather than approximated — volume and area on
/// the same solid keep working, which is why this is not a [`MassProps`] field.
pub fn centroid(model: &Model, solid: Handle<Solid>) -> Result<Point3, PropsError> {
    let solid = model.solid(solid);
    let outer = model.shell(solid.outer);
    let first_face = model.face(outer.faces[0]);
    let reference = model.vertex_point(he_start(model, first_face.outer.half_edges[0])?);

    let mut volume = 0.0;
    let mut moment = Vector3::zero();
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &face in &model.shell(sh).faces {
            let f = face_props(model, face)?;
            let normal = f.normal.ok_or(PropsError::CentroidOfCurvedFace)?;
            let arm = f.centroid - reference;
            let v = f.area * normal.dot(arm) / 3.0;
            volume += v;
            moment += arm * (0.75 * v);
        }
    }
    if volume == 0.0 {
        return Err(PropsError::UnsupportedBoundary);
    }
    Ok(reference + moment * (1.0 / volume))
}

/// Axis-aligned bounding box of one solid, as `(min, max)`.
///
/// **Curve-aware, not a hull of the vertices.** A cylinder's lateral face bulges
/// past its two seam vertices, so a vertex sweep would report a box that is
/// silently too small. A full circle of radius `r` and normal `n̂` extends
/// `±r·√(1 − (n̂·e)²)` along each axis `e`, which is exact; a cylinder band is
/// bounded by its two rim circles, so walking every edge covers it.
pub fn bounds(model: &Model, solid: Handle<Solid>) -> Result<(Point3, Point3), PropsError> {
    let solid = model.solid(solid);
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut grow = |p: Point3, pad: [f64; 3]| {
        let c = p.as_array();
        for i in 0..3 {
            lo[i] = f64::min(lo[i], c[i] - pad[i]);
            hi[i] = f64::max(hi[i], c[i] + pad[i]);
        }
    };
    // The outer shell bounds the solid; a cavity lies inside it by construction.
    for &face in &model.shell(solid.outer).faces {
        let f = model.face(face);
        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
            for &he in &lp.half_edges {
                match edge_curve(model, he) {
                    Curve::Line(_) => grow(model.vertex_point(he_start(model, he)?), [0.0; 3]),
                    Curve::Circle(circle) => {
                        // ★ The **whole** circle's extent, whatever the edge's endpoints say — so
                        // an arc grows the box too much, never too little. A bound that
                        // errs outward stays a bound, which is why this one needs no change when
                        // arcs arrive; `sample_edge` in `nacre-tess` is the one that would not.
                        let n = circle.normal().as_array();
                        let r = circle.radius();
                        let pad = std::array::from_fn(|i| r * (1.0 - n[i] * n[i]).max(0.0).sqrt());
                        grow(circle.center(), pad);
                    }
                }
            }
        }
    }
    if lo[0] > hi[0] {
        return Err(PropsError::UnsupportedBoundary);
    }
    Ok((Point3::from_array(lo), Point3::from_array(hi)))
}

/// `(area, raw volume flux)` of one face — flux is `∮ (r − R)·n̂ dA` without the
/// global `1/3`. Outward normal = surface natural normal × orientation sign.
fn face_contribution(
    model: &Model,
    face: &Face,
    reference: Point3,
) -> Result<(f64, f64), PropsError> {
    let sign = f64::from(face.orientation.sign());
    // ★ The cache, for the same reason as `face_props`: the divergence theorem is being
    // evaluated in `f64`, so the realization is the right description to read.
    match model.surface_cache(face.surface) {
        Surface::Plane(plane) => {
            // Outer boundary, minus each inner loop (a hole): area and first
            // moment are additive, so both the area and the flux subtract the
            // hole's contribution (divergence theorem, n̂ constant on a plane).
            let normal = plane.normal() * sign;
            let (area, centroid) = planar_region(model, face)?;
            Ok((area, normal.dot(centroid - reference) * area))
        }
        Surface::Cylinder(cyl) => {
            // With `p = a₀ + z·axis + r·n̂(θ)` and outward `sign·n̂`, `(p−R)·n̂ = r + (a₀−R)·n̂`
            // (the axial term drops: `axis·n̂ = 0`), so over the region
            //   ∮(p−R)·n̂_out dA = sign·r·( r·∬1 + (a₀−R)·∬n̂ ).
            // ★ The moments below are already `sign`-carrying, so the two signs **cancel** and
            // none is written here; the area, which must come out positive, keeps one.
            let m = lateral_moments(model, face, cyl)?;
            let radius = cyl.radius();
            let a0 = cyl.axis().origin();
            Ok((
                radius * sign * m.j1,
                radius * (radius * m.j1 + (a0 - reference).dot(m.jn)),
            ))
        }
    }
}

/// Area and area-weighted centroid of a planar face's **region** — its outer
/// loop with every hole removed. Area and first moment are both additive, so a
/// hole subtracts each.
fn planar_region(model: &Model, face: &Face) -> Result<(f64, Point3), PropsError> {
    let (mut area, c_out) = planar_face(model, &face.outer)?;
    let mut moment = Vector3::zero();
    for hole in &face.inner {
        let (a_in, c_in) = planar_face(model, hole)?;
        area -= a_in;
        moment -= (c_in - c_out) * a_in;
    }
    if area == 0.0 {
        return Err(PropsError::UnsupportedBoundary);
    }
    Ok((area, c_out + moment * (1.0 / area)))
}

/// Area and area-weighted centroid of one planar loop: a straight polygon, or a
/// single full circle (a cylinder cap).
fn planar_face(model: &Model, outer: &Loop) -> Result<(f64, Point3), PropsError> {
    let all_lines = outer
        .half_edges
        .iter()
        .all(|he| matches!(edge_curve(model, *he), Curve::Line(_)));

    if all_lines {
        let points = loop_points(model, outer)?;
        // The signed-area-weighted centroid lives in `nacre-geom` (pure geometry, and
        // `nacre-ops` needs the same computation for its sketch frames).
        nacre_geom::planar_region_area_centroid(&points, &[]).ok_or(PropsError::UnsupportedBoundary)
    } else if outer.half_edges.len() == 1 {
        // ★ **`len() == 1` is what makes the whole-circle formula safe when arcs arrive**, and the
        // reason is structural rather than a check: an arc has two *different* endpoints, so it
        // cannot close a loop alone — the smallest arc-bounded loop is an arc plus a chord, two
        // half-edges, which lands in the `else` and declines. Written down so the next reader does
        // not "fix" a formula that is not wrong. (A loop *mixing* arcs with lines is the same
        // decline, and `nacre-validate::loop_winding` names it as its own landing site.)
        match edge_curve(model, outer.half_edges[0]) {
            Curve::Circle(circle) => {
                let area = PI * circle.radius() * circle.radius();
                Ok((area, circle.center()))
            }
            _ => Err(PropsError::UnsupportedBoundary),
        }
    } else {
        // ★ **The mixed arm**: a loop of chords and circular arcs. The chord polygon is
        // a signed fan from the first vertex; each arc half-edge then adds or removes its
        // **circular segment** — area `r²(Δθ − sin Δθ)/2`, centroid `4r·sin³(Δθ/2) /
        // (3(Δθ − sin Δθ))` out along the bisector — with Δθ read from the edge's stored
        // `[from, to]` order (**CCW about the axis**, the arc convention;
        // [`Circle::angle_of`] is the one spelling of the angle; a closed edge `[v, v]` is the
        // whole turn) and the sign from which way this traversal walks it. A digon (chord + arc) is the degenerate case the fan
        // contributes nothing to: one segment is the whole answer.
        //
        // Every contribution is scalarized on one reference normal — the first arc's circle
        // normal, which is the cylinder axis and so collinear with the cap plane's normal
        // either way up — so segments and fan combine by plain signs, and the signs cancel out
        // of `moment / area` (the centroid is winding-independent, like the polygon road's).
        let points = loop_points(model, outer)?;
        let n_ref = outer
            .half_edges
            .iter()
            .find_map(|&he| match edge_curve(model, he) {
                Curve::Circle(c) => Some(c.normal()),
                _ => None,
            })
            .ok_or(PropsError::UnsupportedBoundary)?;
        let base = points[0];
        let mut area = 0.0;
        let mut moment = Vector3::zero(); // Σ aᵢ·(cᵢ − base)
        for i in 1..points.len().saturating_sub(1) {
            let (u, v) = (points[i] - base, points[i + 1] - base);
            let a = 0.5 * u.cross(v).dot(n_ref);
            area += a;
            // (base + pᵢ + pᵢ₊₁)/3 − base = (u + v)/3.
            moment += (u + v) * (a / 3.0);
        }
        for he in &outer.half_edges {
            let Curve::Circle(c) = edge_curve(model, *he) else {
                continue;
            };
            let [va, vb] = model.edge(he.edge).vertices;
            let t0 = c.angle_of(model.vertex_point(va));
            let t1 = c.angle_of(model.vertex_point(vb));
            // A closed edge `[v, v]` is the whole turn — its ends cannot say so.
            let dt = if va == vb {
                std::f64::consts::TAU
            } else {
                (t1 - t0).rem_euclid(std::f64::consts::TAU)
            };
            let seg = c.segment_area(dt);
            if seg <= 0.0 {
                continue; // a zero span contributes nothing
            }
            let sign = if he.forward { 1.0 } else { -1.0 } * c.normal().dot(n_ref).signum();
            let c_seg = c.segment_centroid(t0, dt);
            area += sign * seg;
            moment += (c_seg - base) * (sign * seg);
        }
        if area == 0.0 {
            return Err(PropsError::UnsupportedBoundary);
        }
        let centroid = base + moment * (1.0 / area);
        Ok((area.abs(), centroid))
    }
}

/// The `(θ, z)` moments of a lateral face's **region**, read off its boundary.
///
/// Each is carried with the sign of the walk's sense in the chart, which is the face's
/// `Orientation::sign()` (see [`lateral_moments`]); consumers either multiply it back in (the
/// area) or take a ratio of two moments (the centroid) or have it cancel against the outward
/// normal's own sign (the flux).
struct LateralMoments {
    /// `∬ dθ dz` — the chart area.
    j1: f64,
    /// `∬ n̂(θ) dθ dz`, in world coordinates.
    jn: Vector3,
    /// `∬ (z − z0) dθ dz`.
    jz: f64,
    /// The axial station, from the axis origin, that `jz` measures `z` from.
    z0: f64,
}

/// A lateral face's moments, by **one line integral over every loop it has**.
///
/// ★★★★★ **This is where the rectangle assumption died.** The old spelling found the loop's two
/// rims, read one `Δθ` and one height off them, and multiplied — true only for a `(θ, z)`
/// *rectangle*, so a face with a hole, or one whose outer walk detours around a notch, had to be
/// refused (`UnsupportedInnerLoop`, `UnsupportedBoundary`). Green's theorem asks the boundary
/// instead: for any region `R` in the chart and any `g` with `∂G/∂z = g`,
///
/// ```text
///     ∬_R g dθ dz  =  ∮_∂R −G dθ        (∂R walked CCW in the chart)
/// ```
///
/// and the four `g` this crate needs — `1`, `cos θ`, `sin θ` (area and flux) and `z` (the
/// centroid) — all have closed-form `G`. Three consequences carry the curved arms:
///
/// - **Inner loops are not a special case.** The identity sums over *every* boundary component in
///   the sense it is wound, so a hole subtracts itself. No `UnsupportedInnerLoop`.
/// - **A straight edge contributes nothing**: a ruling sits at one θ, so `dθ = 0`.
/// - **The sign is the walk's.** `θ̂ × ẑ = r̂`, so the chart is right-handed about the *surface's*
///   normal; a loop is CCW about its *face's* outward normal, which is `sign · r̂`. Hence every
///   moment carries `sign`, and the area's `sign · j1` must come out **positive** — it is not
///   wrapped in `abs()`, because a mis-threaded walk showing up as a negative area is the only
///   cheap net for one.
///
/// Honest declines (`UnsupportedBoundary`): an empty loop, or a region whose chart area is not
/// positive.
///
/// ☑ Unchanged assumption: a rim is a **circle** cut by a plane ⊥ the axis, so `z` is constant
/// along an arc. An oblique cut gives an ellipse, and that is out of scope here.
fn lateral_moments(
    model: &Model,
    face: &Face,
    cyl: &nacre_geom::Cylinder,
) -> Result<LateralMoments, PropsError> {
    let axis = cyl.axis().direction();
    let a0 = cyl.axis().origin();
    let (u, w) = (cyl.ref_dir(), axis.cross(cyl.ref_dir()));
    let theta_of = |p: Point3| -> f64 {
        let rel = p - a0;
        let q = rel - axis * rel.dot(axis);
        q.dot(w).atan2(q.dot(u))
    };
    let axial_of = |p: Point3| (p - a0).dot(axis);
    // ★ `z` is measured from the face's first vertex, to keep the magnitudes small wherever the
    // solid sits. It is legal because `∮cos θ dθ` and `∮sin θ dθ` vanish over every closed loop
    // and `∮dθ` over the face's loops together (a rim loop alone turns ±2π; the two rims turn
    // opposite ways), so `z → z + c` moves neither `j1` nor `jn`; only `jz` shifts, and `z0`
    // travels with it to the one place it is read.
    let first = *face
        .outer
        .half_edges
        .first()
        .ok_or(PropsError::UnsupportedBoundary)?;
    let z0 = axial_of(model.vertex_point(he_start(model, first)?));
    let (mut j1, mut jz, mut jc, mut js) = (0.0, 0.0, 0.0, 0.0);
    for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
        for &he in &lp.half_edges {
            let Curve::Circle(c) = edge_curve(model, he) else {
                continue; // straight: a ruling, at one θ ⇒ dθ = 0
            };
            let z = axial_of(c.center()) - z0;
            let [va, vb] = model.edge(he.edge).vertices;
            // ★★★★★ **The travel comes from the edge's own stored order, never from a wrapped
            // angle difference.** An arc edge stores `[from, to]` **CCW about the axis**
            // (`EdgeKey::Arc`), and `he.forward` says which way this face walks it. A wrapped
            // `θ_to − θ_from` is `±π` for a **half-turn** arc with no way to tell which — and a
            // boss standing on a wall cuts its own rim into exactly two half turns, so the
            // ambiguous case is not a corner of the population, it *is* the population.
            let (travel, t_start, t_end) = if va == vb {
                // A closed rim: one full turn, and its trig moments are **exactly** zero —
                // written as zero rather than differenced, which is what keeps the full-band
                // populations on their old numbers.
                let turn = std::f64::consts::TAU;
                (if he.forward { turn } else { -turn }, 0.0, 0.0)
            } else {
                let (from, to) = (
                    theta_of(model.vertex_point(va)),
                    theta_of(model.vertex_point(vb)),
                );
                let span = (to - from).rem_euclid(std::f64::consts::TAU);
                if he.forward {
                    (span, from, to)
                } else {
                    (-span, to, from)
                }
            };
            j1 -= z * travel;
            jz -= 0.5 * z * z * travel;
            if va != vb {
                // `sin`/`cos` are periodic, so these need no unwrapping — the same argument the
                // chain-end spelling rested on.
                jc -= z * (t_end.sin() - t_start.sin());
                js += z * (t_end.cos() - t_start.cos());
            }
        }
    }
    if j1 == 0.0 {
        return Err(PropsError::UnsupportedBoundary);
    }
    Ok(LateralMoments {
        j1,
        jn: u * jc + w * js,
        jz,
        z0,
    })
}

/// The ordered start vertices of a loop's half-edges (the boundary polyline).
fn loop_points(model: &Model, outer: &Loop) -> Result<Vec<Point3>, PropsError> {
    outer
        .half_edges
        .iter()
        .map(|&he| Ok(model.vertex_point(he_start(model, he)?)))
        .collect()
}

/// [`Model::he_start`], kept as a local name. (The `Result` is always `Ok`: every
/// edge is bounded by type, so the unsupported-boundary arm has nothing to catch. The wrapper
/// stays because five callers still `?` it; collapsing it is a separate change.)
fn he_start(
    model: &Model,
    he: nacre_topo::HalfEdge,
) -> Result<Handle<nacre_topo::Vertex>, PropsError> {
    Ok(model.he_start(he))
}

/// The `Curve` carried by a half-edge's edge.
fn edge_curve(model: &Model, he: nacre_topo::HalfEdge) -> &Curve {
    model.edge_curve(he.edge)
}

#[cfg(test)]
#[path = "tests/lib.rs"]
mod tests;
