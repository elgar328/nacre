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
//! Consumers: user "part volume / surface area" queries, the M5 boolean
//! volume-conservation invariants (proptest), and the OCCT oracle diff
//! (`nacre-oracle`).
//!
//! Precondition: the solid is a valid closed, outward-oriented b-rep (what the
//! producers emit and `validate` accepts). Coordinates are the cache side of the
//! truth/cache split (design.md §0), so f64 arithmetic here is appropriate.

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
    /// A **curved** face carries an inner loop (a hole). Planar holes are
    /// supported; a hole in a cylindrical face has no producer yet, so it is
    /// rejected rather than mis-integrated.
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
    let solid = model.solids.get(solid);
    let outer = model.shells.get(solid.outer);

    // R = the first vertex of the first face of the outer shell; any vertex on
    // the solid works. The closed-surface identity ∮ n̂ dA = 0 holds over the
    // *full* boundary — outer shell plus every cavity shell — so V is
    // R-independent; keeping R on the outer shell keeps the numbers small.
    let first_face = model.faces.get(outer.faces[0]);
    let reference = model.vertex_point(he_start(model, first_face.outer.half_edges[0])?);

    let mut volume_flux = 0.0;
    let mut area = 0.0;
    // Outer boundary, then each inner cavity shell. A cavity's faces are wound
    // with their outward normals pointing into the void (design §8 M5
    // containment), so its flux is negative and subtracts the void's volume;
    // its (unsigned) area adds — both surfaces bound material.
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &face in &model.shells.get(sh).faces {
            let (a, flux) = face_contribution(model, model.faces.get(face), reference)?;
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
    let f = model.faces.get(face);
    let sign = f64::from(f.orientation.sign());
    Some(model.surface(f.surface).normal_at(p)? * sign)
}

/// [`FaceProps`] of one face.
pub fn face_props(model: &Model, face: Handle<Face>) -> Result<FaceProps, PropsError> {
    let face = model.faces.get(face);
    let sign = f64::from(face.orientation.sign());
    match model.surface(face.surface) {
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
            // mean axial station, plus `r·∬n̂ / ∬1` radially (exactly zero for a full band —
            // the old on-axis answer). Every `sign` cancels in a ratio of two moments.
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
    let solid = model.solids.get(solid);
    let outer = model.shells.get(solid.outer);
    let first_face = model.faces.get(outer.faces[0]);
    let reference = model.vertex_point(he_start(model, first_face.outer.half_edges[0])?);

    let mut volume = 0.0;
    let mut moment = Vector3::zero();
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &face in &model.shells.get(sh).faces {
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
    let solid = model.solids.get(solid);
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
    for &face in &model.shells.get(solid.outer).faces {
        let f = model.faces.get(face);
        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
            for &he in &lp.half_edges {
                match edge_curve(model, he) {
                    Curve::Line(_) => grow(model.vertex_point(he_start(model, he)?), [0.0; 3]),
                    Curve::Circle(circle) => {
                        // ★ The **whole** circle's extent, whatever the edge's endpoints say — so
                        // an arc (M6-2b) grows the box too much, never too little. A bound that
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
    match model.surface(face.surface) {
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
        // ★ **The mixed arm** (M6-2b): a loop of chords and circular arcs. The chord polygon is
        // a signed fan from the first vertex; each arc half-edge then adds or removes its
        // **circular segment** — area `r²(Δθ − sin Δθ)/2`, centroid `4r·sin³(Δθ/2) /
        // (3(Δθ − sin Δθ))` out along the bisector — with Δθ read from the edge's stored
        // `[from, to]` order (**CCW about the axis**, the M6-2b convention;
        // [`Circle::angle_of`] is the one spelling of the angle) and the sign from which way
        // this traversal walks it. A digon (chord + arc) is the degenerate case the fan
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
            let [va, vb] = model.edges.get(he.edge).vertices;
            let t0 = c.angle_of(model.vertex_point(va));
            let t1 = c.angle_of(model.vertex_point(vb));
            let dt = (t1 - t0).rem_euclid(std::f64::consts::TAU);
            let seg = c.segment_area(dt);
            if seg <= 0.0 {
                continue; // a degenerate (closed or zero) span contributes nothing
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
/// - **A straight edge contributes nothing**: a ruling and the seam both sit at one θ, so
///   `dθ = 0`. An outer walk that bridges a hole along the seam therefore reads the same as if
///   the hole were a separate loop — the bridge is traversed twice, and each traversal is zero.
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
/// along an arc. An oblique cut gives an ellipse, and that is M6-3's business — the old spelling
/// grouped rim arcs by axial station on the same premise.
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
    // ★ `z` is measured from the face's first vertex, for the reason the old `axial_range` did the
    // same: keep the magnitudes small wherever the solid sits. It is legal because `∮dθ`, `∮cos θ
    // dθ` and `∮sin θ dθ` all vanish over a closed loop, so `z → z + c` moves neither `j1` nor
    // `jn`; only `jz` shifts, and `z0` travels with it to the one place it is read.
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
                continue; // straight: a ruling or a seam piece, at one θ ⇒ dθ = 0
            };
            let z = axial_of(c.center()) - z0;
            let [va, vb] = model.edges.get(he.edge).vertices;
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

/// [`Model::he_start`], kept as a local name. (Its `Result` wrapper died with S8: every edge
/// is bounded by type, so the unsupported-boundary arm had nothing left to catch.)
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
mod tests {
    use super::*;
    use nacre_math::{Point2, Vector3};
    use nacre_ops::SketchFrame;
    use nacre_ops::{OpOutput, Operation, Profile2d, apply};
    use nacre_scalar::Axis;

    /// Relative-or-absolute comparison sized for accumulated f64 error.
    fn close(got: f64, expected: f64) -> bool {
        (got - expected).abs() <= 1e-9 * expected.abs().max(1.0)
    }

    // --- golden ---

    #[test]
    fn cube_mass_matches_analytic() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, 24.0), "volume {}", props.volume); // 2·3·4
        assert!(close(props.area, 52.0), "area {}", props.area); // 2(6+8+12)
    }

    #[test]
    fn cylinder_mass_matches_analytic() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, 20.0 * PI), "volume {}", props.volume); // πr²h
        assert!(close(props.area, 28.0 * PI), "area {}", props.area); // 2πr² + 2πrh
    }

    /// ★ **Two doors, one fact.** A planar face's normal is answered twice — by
    /// [`face_props`] for the whole face and by [`face_normal_at`] at a point — and the two
    /// must say the same thing, or a caller would get a different answer depending on which
    /// it happened to ask. (The point plays no part on a plane, which this also shows.)
    #[test]
    fn the_two_normal_doors_agree_on_a_plane() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        m.rebuild_adjacency();
        let shell = m.solids.get(s).outer;
        for &fh in &m.shells.get(shell).faces {
            let whole = face_props(&m, fh)
                .unwrap()
                .normal
                .expect("a box face is planar");
            // Two different points of the same face — its own centroid and a corner.
            let corner = m.vertex_point(he_start(&m, m.faces.get(fh).outer.half_edges[0]).unwrap());
            for p in [face_props(&m, fh).unwrap().centroid, corner] {
                let at = face_normal_at(&m, fh, p).expect("a planar face has a normal anywhere");
                assert!((at - whole).norm() < 1e-12, "{at:?} vs {whole:?}");
            }
        }
    }

    /// The cylinder branch of the same door, where the point is the whole question: a free
    /// cylinder's wall faces away from its axis, everywhere on it.
    #[test]
    fn a_free_cylinders_wall_faces_away_from_its_axis() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        m.rebuild_adjacency();
        let shell = m.solids.get(s).outer;
        let wall = *m
            .shells
            .get(shell)
            .faces
            .iter()
            .find(|&&fh| matches!(m.surface(m.faces.get(fh).surface), Surface::Cylinder(_)))
            .expect("a lateral face");
        assert!(
            face_props(&m, wall).unwrap().normal.is_none(),
            "a curved face has no single normal — that is why the point door exists"
        );
        let cyl = match m.surface(m.faces.get(wall).surface) {
            Surface::Cylinder(c) => *c,
            _ => unreachable!(),
        };
        for k in 0..8 {
            let u = std::f64::consts::TAU * f64::from(k) / 8.0;
            let p = cyl.point_at(u, 2.5);
            let n = face_normal_at(&m, wall, p).expect("off the axis");
            let radial = p - cyl.axis().origin();
            let radial = radial - cyl.axis().direction() * radial.dot(cyl.axis().direction());
            assert!(
                n.dot(radial) > 0.0,
                "u={u}: the wall faces away from the axis"
            );
            assert!((n.norm() - 1.0).abs() < 1e-12, "unit");
        }
    }

    #[test]
    fn concave_extrude_mass() {
        // An L-shaped profile: a 2×2 square with the top-right 1×1 corner removed.
        // Concave, so the area-weighted centroid differs from the vertex average —
        // a vertex-average bug in the flux would fail the volume assertion.
        let pts = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let profile =
            Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).unwrap();
        let dist = 3.0;
        let mut m = Model::new();
        let op = Operation::Extrude {
            frame: SketchFrame::world(&m, Axis::Z),
            profile,
            dist,
        };
        let OpOutput::Extrude { solid: s, .. } = apply(&mut m, &op).unwrap() else {
            unreachable!("extrude yields Extrude output");
        };

        // Independent expected base area via the 2D shoelace (not mass_props).
        let a_l = shoelace(&pts);
        assert!((a_l - 3.0).abs() < 1e-12); // 4 − 1

        let props = mass_props(&m, s).unwrap();
        assert!(close(props.volume, a_l * dist), "volume {}", props.volume);
    }

    // --- bounds / face_props / centroid ---

    /// **The trap this exists for.** A cylinder's lateral face bulges past its two
    /// seam vertices, so a hull of the vertices reports a box that is too small in
    /// the radial directions — a wrong answer with nothing to tip the caller off.
    /// The axis is deliberately oblique so the analytic term is not axis-aligned.
    #[test]
    fn bounds_follow_the_curve_not_the_vertices() {
        let axis = Vector3::from_array([0.0, 3.0, 4.0]); // unit (0, .6, .8)
        let mut m = Model::new();
        let s = m.add_cylinder(Point3::origin(), axis, 2.0, 5.0);
        let (lo, hi) = bounds(&m, s).unwrap();

        // x is fully perpendicular to the axis, so the barrel spans the diameter.
        assert!(close(lo[0], -2.0) && close(hi[0], 2.0), "{lo:?} {hi:?}");
        // Along y and z the circles foreshorten by √(1 − (n̂·e)²).
        let u = axis.normalize().unwrap().as_array();
        let end = Point3::origin() + axis.normalize().unwrap() * 5.0;
        for i in [1, 2] {
            let pad = 2.0 * (1.0 - u[i] * u[i]).sqrt();
            assert!(close(lo[i], -pad), "axis {i} lo {}", lo[i]);
            assert!(
                close(hi[i], end.as_array()[i] + pad),
                "axis {i} hi {}",
                hi[i]
            );
        }

        // And the vertex hull really is smaller — otherwise this test proves nothing.
        let mut vlo = [f64::INFINITY; 3];
        let mut vhi = [f64::NEG_INFINITY; 3];
        for (vh, _) in m.vertices.iter() {
            let p = m.vertex_point(vh).as_array();
            for i in 0..3 {
                vlo[i] = vlo[i].min(p[i]);
                vhi[i] = vhi[i].max(p[i]);
            }
        }
        assert!(
            vlo[0] > lo[0] + 1.0,
            "the vertex hull must be visibly wrong here: {vlo:?} vs {lo:?}"
        );
    }

    #[test]
    fn bounds_of_a_box_are_its_corners() {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([-1.0, 2.0, 0.5]),
            Point3::from_array([3.0, 4.0, 9.0]),
        );
        let (lo, hi) = bounds(&m, s).unwrap();
        assert!(
            close(lo[0], -1.0) && close(lo[1], 2.0) && close(lo[2], 0.5),
            "{lo:?}"
        );
        assert!(
            close(hi[0], 3.0) && close(hi[1], 4.0) && close(hi[2], 9.0),
            "{hi:?}"
        );
    }

    /// The query a script needs to *name* a face: filter by normal, then by position.
    #[test]
    fn faces_can_be_picked_by_their_geometry() {
        let mut m = Model::new();
        let s = m.add_cuboid(Point3::origin(), Point3::from_array([2.0, 3.0, 4.0]));
        let up = Vector3::from_array([0.0, 0.0, 1.0]);
        let top = m
            .shells
            .get(m.solids.get(s).outer)
            .faces
            .iter()
            .map(|&f| (f, face_props(&m, f).unwrap()))
            .filter(|(_, p)| p.normal.is_some_and(|n| n.dot(up) > 0.5))
            .max_by(|a, b| a.1.centroid[2].partial_cmp(&b.1.centroid[2]).unwrap())
            .unwrap();
        assert!(close(top.1.area, 6.0), "area {}", top.1.area); // 2·3
        assert!(close(top.1.centroid[2], 4.0), "z {}", top.1.centroid[2]);
        assert!(close(top.1.centroid[0], 1.0) && close(top.1.centroid[1], 1.5));
    }

    /// A hole moves the face's centroid; using the outer ring alone would not.
    #[test]
    fn a_face_with_an_off_centre_hole_reports_the_region() {
        let sq = |a: f64, b: f64| {
            vec![
                Point2::from_array([a, a]),
                Point2::from_array([b, a]),
                Point2::from_array([b, b]),
                Point2::from_array([a, b]),
            ]
        };
        let mut m = Model::new();
        let op = Operation::Extrude {
            frame: SketchFrame::world(&m, Axis::Z),
            profile: Profile2d::with_holes(sq(0.0, 10.0), vec![sq(1.0, 3.0)]).unwrap(),
            dist: 1.0,
        };
        let OpOutput::Extrude { faces, .. } = apply(&mut m, &op).unwrap() else {
            unreachable!()
        };
        // faces[1] is the top cap — the one carrying the inner loop.
        let p = face_props(&m, faces[1]).unwrap();
        assert!(close(p.area, 96.0), "area {}", p.area); // 100 − 4
        let want = (100.0 * 5.0 - 4.0 * 2.0) / 96.0;
        assert!(
            close(p.centroid[0], want) && close(p.centroid[1], want),
            "{:?}",
            p.centroid
        );
    }

    /// Asymmetric on purpose: a symmetric solid puts the centroid at the middle
    /// whatever the arithmetic does, so it proves nothing.
    #[test]
    fn centroid_of_an_l_prism_is_not_its_bounding_centre() {
        let pts = [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let mut m = Model::new();
        let op = Operation::Extrude {
            frame: SketchFrame::world(&m, Axis::Z),
            profile: Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect())
                .unwrap(),
            dist: 3.0,
        };
        let OpOutput::Extrude { solid: s, .. } = apply(&mut m, &op).unwrap() else {
            unreachable!()
        };
        let c = centroid(&m, s).unwrap();
        // Area 3 = unit squares at (.5,.5), (1.5,.5), (.5,1.5) ⇒ centroid (5/6, 5/6).
        assert!(close(c[0], 5.0 / 6.0) && close(c[1], 5.0 / 6.0), "{c:?}");
        assert!(close(c[2], 1.5), "{c:?}");
    }

    /// A cavity carries the opposite sign, so an off-centre void must push the
    /// centroid away from it. Dropping that sign is invisible on a centred void.
    #[test]
    fn an_off_centre_void_pushes_the_centroid_away() {
        let mut m = Model::new();
        let outer = m.add_cuboid(Point3::origin(), Point3::from_array([10.0; 3]));
        let inner = m.add_cuboid(
            Point3::from_array([1.0; 3]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        let r = nacre_ops::boolean(&mut m, nacre_ops::BoolKind::Cut, outer, inner).unwrap();
        m.rebuild_adjacency();
        assert_eq!(r.len(), 1);
        let c = centroid(&m, r[0]).unwrap();
        // (1000·5 − 8·2)/992 on every axis: the void at (2,2,2) drags the rest up.
        let want = (1000.0 * 5.0 - 8.0 * 2.0) / 992.0;
        assert!(
            want > 5.0,
            "the void must be off-centre for this to prove anything"
        );
        for i in 0..3 {
            assert!(close(c[i], want), "axis {i}: {} vs {want}", c[i]);
        }
    }

    #[test]
    fn a_curved_solid_refuses_a_centroid_but_still_reports_volume() {
        let mut m = Model::new();
        let s = m.add_cylinder(
            Point3::origin(),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            5.0,
        );
        assert!(matches!(
            centroid(&m, s),
            Err(PropsError::CentroidOfCurvedFace)
        ));
        assert!(close(mass_props(&m, s).unwrap().volume, 20.0 * PI));
    }

    /// A `2·hw` square centred on a `size` cube's lid.
    ///
    /// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
    /// onto the face's plane and the axes are `u = +x_hat`, `v = +y_hat`, so a frame point `(a, b)`
    /// is world `(a, b, size)`.
    fn centred_on_the_lid(size: f64, hw: f64) -> Profile2d {
        let (cx, cy) = (0.5 * size, 0.5 * size);
        Profile2d::polygon(
            [
                [cx - hw, cy + hw],
                [cx - hw, cy - hw],
                [cx + hw, cy - hw],
                [cx + hw, cy + hw],
            ]
            .iter()
            .map(|&p| Point2::from_array(p))
            .collect(),
        )
        .unwrap()
    }

    /// Pad a `2·hw` square boss of height `dist` on a `size` cube's top face,
    /// returning the padded solid's mass.
    fn cube_then_pad(size: f64, hw: f64, dist: f64) -> MassProps {
        let sq = |s: f64| {
            Profile2d::polygon(
                [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                    .iter()
                    .map(|&p| Point2::from_array(p))
                    .collect(),
            )
            .unwrap()
        };
        let mut m = Model::new();
        let __w1 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w1,
                profile: sq(size),
                dist: size,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let boss = centred_on_the_lid(size, hw);
        let OpOutput::PadOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: faces[1],
                profile: boss,
                dist,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        mass_props(&m, solid).unwrap()
    }

    #[test]
    fn pad_boss_mass() {
        // Unit cube + a 0.4-square boss (area 0.16, perimeter 1.6) of height 0.5.
        // Volume = 1 + 0.16·0.5 = 1.08; area = 6 + 1.6·0.5 = 6.8. This is the
        // *absolute* test of the inner-loop subtraction (the hole is not filled).
        let m = cube_then_pad(1.0, 0.2, 0.5);
        assert!(close(m.volume, 1.08), "vol {}", m.volume);
        assert!(close(m.area, 6.8), "area {}", m.area);
    }

    /// Carve a `2·hw` square pocket of depth `dist` into a `size` cube's top face.
    fn cube_then_pocket(size: f64, hw: f64, dist: f64) -> MassProps {
        let sq = |s: f64| {
            Profile2d::polygon(
                [[0.0, 0.0], [s, 0.0], [s, s], [0.0, s]]
                    .iter()
                    .map(|&p| Point2::from_array(p))
                    .collect(),
            )
            .unwrap()
        };
        let mut m = Model::new();
        let __w0 = SketchFrame::world(&m, Axis::Z);
        let OpOutput::Extrude { faces, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame: __w0,
                profile: sq(size),
                dist: size,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        let pocket = centred_on_the_lid(size, hw);
        let OpOutput::PocketOnFace { solid, .. } = apply(
            &mut m,
            &Operation::PocketOnFace {
                face: faces[1],
                profile: pocket,
                dist,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        mass_props(&m, solid).unwrap()
    }

    #[test]
    fn pocket_mass() {
        // Unit cube − a 0.4-square pocket of depth 0.5. Volume = 1 − 0.16·0.5 =
        // 0.92 (the inward walls contribute negatively); area = 6 + 1.6·0.5 = 6.8
        // (same as the boss). Exercises the hole subtraction with inward walls.
        let m = cube_then_pocket(1.0, 0.2, 0.5);
        assert!(close(m.volume, 0.92), "vol {}", m.volume);
        assert!(close(m.area, 6.8), "area {}", m.area);
    }

    /// Build a hollow solid: an `outer`-cube with a concentric `inner`-cube void,
    /// the inner cube's shell reversed inward (M5 containment). Returns its mass.
    fn cube_in_cube(min: Point3, outer: f64, inner: f64) -> MassProps {
        let mut m = Model::new();
        let ext = |s: f64| Vector3::from_array([s, s, s]);
        let a = m.add_cuboid(min, min + ext(outer));
        let gap = 0.5 * (outer - inner); // centered ⇒ strictly interior on all sides
        let inner_min = min + ext(gap);
        let b = m.add_cuboid(inner_min, inner_min + ext(inner));

        let b_outer = m.solids.get(b).outer;
        let void = m.reversed_shell(b_outer);
        let a_outer = m.solids.get(a).outer;
        let hollow = m.push_solid(Solid {
            outer: a_outer,
            cavities: vec![void],
        });
        m.live_solids.retain(|&s| s == hollow); // supersede the two source cubes
        mass_props(&m, hollow).unwrap()
    }

    #[test]
    fn cube_in_cube_subtracts_the_void() {
        // 4-cube with a concentric 2-cube void: V = 4³ − 2³ = 56; total surface
        // = outer 6·4² + void 6·2² = 96 + 24 = 120 (both surfaces bound material).
        let props = cube_in_cube(Point3::origin(), 4.0, 2.0);
        assert!(close(props.volume, 56.0), "vol {}", props.volume);
        assert!(close(props.area, 120.0), "area {}", props.area);
    }

    /// Unsigned area of a 2D polygon (independent check for the concave test).
    fn shoelace(pts: &[[f64; 2]]) -> f64 {
        let mut two_area = 0.0;
        for i in 0..pts.len() {
            let a = pts[i];
            let b = pts[(i + 1) % pts.len()];
            two_area += a[0] * b[1] - b[0] * a[1];
        }
        0.5 * two_area.abs()
    }

    // --- proptest ---

    use proptest::prelude::*;

    fn wide_point() -> impl Strategy<Value = Point3> {
        prop::array::uniform3(-1e6f64..1e6).prop_map(Point3::from_array)
    }

    proptest! {
        /// Random box, possibly far from the origin — exercises the R
        /// cancellation guard as well as the polygon path.
        #[test]
        fn prop_cuboid_mass(
            min in wide_point(),
            a in 0.1f64..1e3, b in 0.1f64..1e3, c in 0.1f64..1e3,
        ) {
            let max = min + Vector3::from_array([a, b, c]);
            let mut m = Model::new();
            let s = m.add_cuboid(min, max);
            let props = mass_props(&m, s).unwrap();
            prop_assert!(close(props.volume, a * b * c), "vol {} vs {}", props.volume, a*b*c);
            prop_assert!(close(props.area, 2.0 * (a*b + b*c + c*a)), "area {}", props.area);
        }

        /// Random cylinder with an arbitrary axis direction — exercises the
        /// non-axis-aligned frame, the lateral sign, and the 2π band height.
        #[test]
        fn prop_cylinder_mass(
            base in wide_point(),
            dir in prop::array::uniform3(-1e3f64..1e3)
                .prop_map(Vector3::from_array)
                .prop_filter("axis must be clearly nonzero", |v| v.norm() > 1e-3),
            r in 0.5f64..10.0, h in 0.5f64..1e3,
        ) {
            let mut m = Model::new();
            let s = m.add_cylinder(base, dir, r, h);
            let props = mass_props(&m, s).unwrap();
            prop_assert!(close(props.volume, PI * r * r * h), "vol {} vs {}", props.volume, PI*r*r*h);
            prop_assert!(close(props.area, 2.0 * PI * r * (r + h)), "area {}", props.area);
        }

        /// A boss adds `A_p·dist` of volume and `P·dist` of surface (the top hole
        /// area cancels the boss cap). For a `2·hw` square: `A_p = 4hw²`, `P = 8hw`.
        /// This absolutely exercises the inner-loop subtraction (uncancelled hole).
        #[test]
        fn prop_pad_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..5.0) {
            let m = cube_then_pad(size, hw, dist);
            let vol = size * size * size + 4.0 * hw * hw * dist;
            let area = 6.0 * size * size + 8.0 * hw * dist;
            prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
            prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
        }

        /// A pocket *removes* `A_p·dist` of volume (inward walls) while adding the
        /// same `P·dist` of surface as a boss. `dist ≤ 0.9 < size` keeps it blind.
        #[test]
        fn prop_pocket_volume(size in 1.0f64..5.0, hw in 0.05f64..0.3, dist in 0.1f64..0.9) {
            let m = cube_then_pocket(size, hw, dist);
            let vol = size * size * size - 4.0 * hw * hw * dist;
            let area = 6.0 * size * size + 8.0 * hw * dist;
            prop_assert!(close(m.volume, vol), "vol {} vs {}", m.volume, vol);
            prop_assert!(close(m.area, area), "area {} vs {}", m.area, area);
        }

        /// A concentric cube void of any interior size: V = outer³ − inner³, and
        /// the total surface is the sum of both. `min` ranges far from the origin
        /// to exercise the R-cancellation guard over the *full* (outer + void)
        /// boundary, not just the outer shell.
        #[test]
        fn prop_cube_in_cube(
            min in wide_point(),
            outer in 2.0f64..1e3,
            ratio in 0.1f64..0.85,
        ) {
            let inner = outer * ratio;
            let props = cube_in_cube(min, outer, inner);
            let vol = outer * outer * outer - inner * inner * inner;
            let area = 6.0 * (outer * outer + inner * inner);
            prop_assert!(close(props.volume, vol), "vol {} vs {}", props.volume, vol);
            prop_assert!(close(props.area, area), "area {} vs {}", props.area, area);
        }
    }
}
