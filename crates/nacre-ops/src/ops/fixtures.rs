//! **Solids built through the product's own operations** — test-only (`test-util`).
//!
//! ★ A fixture here states nothing the product does not: it is a sequence of calls an application
//! makes (a datum plane, a circle sketched on it, an extrude). That is the difference from a
//! constructor that writes cells directly — such a door skips the realization funnel, and a
//! population built through it can hold what no operation produces (a tilted cap lifted from
//! realized `f64`, whose rim is an ellipse while its edge cache says circle).
//!
//! ★ **A box is a floor and a height**, the way an application states one (kit's `cuboid`
//! extrudes a rectangle by its size): a datum plane, a rectangle on it, an extrude. Two corners
//! are not a statement the product takes — the extrude reads its distance as the shortest decimal
//! that round-trips (design 「구성 산술은 사용자가 쓴 십진수로 한다」), so a box whose two heights
//! differ by no such decimal cannot be built with its far cap where the corners say.

use super::*;

/// A cylinder built by [`cylinder`] or [`cylinder_with_seam`], with its faces **named** — the
/// extrude's `[base, top, sides…]` order is a fact of `OpOutput::Extrude`, and a caller reading
/// a face by index would be guessing it.
#[derive(Clone, Copy, Debug)]
pub struct CylinderSolid {
    pub solid: Handle<Solid>,
    /// The cap on the sketch plane, facing against the axis.
    pub base: Handle<Face>,
    /// The far cap, facing along the axis.
    pub top: Handle<Face>,
    pub lateral: Handle<Face>,
}

/// A closed cylinder from `base` along `axis` for `height`, of `radius`: the plane through `base`
/// with normal `axis` stated as a datum ([`SketchPlane::from_origin_normal`]), a whole circle
/// about its origin, extruded. The seam sits where the product puts it — the frame's `+x̂`.
///
/// Panics on a statement the product refuses (a zero axis, a value outside the decimal window):
/// a fixture's input is the test's own.
pub fn cylinder(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    let plane = SketchPlane::from_origin_normal(base, axis).expect("a nonzero cylinder axis");
    extruded_circle(m, plane, radius, height)
}

/// [`cylinder`] with the seam written out — for a test whose proposition involves **where the
/// seam is**. The sketch frame is `+x̂ = seam`, `+ŷ = axis × seam`
/// ([`SketchPlane::from_axes`]), so the circle's seam lands on the `seam` side of the axis and
/// the frame's normal is `axis`. `seam` must be perpendicular to `axis`; `axis.any_perpendicular()`
/// is the usual choice. Exact when the two and their cross product are the decimals written (an
/// axis-aligned pair always is).
pub fn cylinder_with_seam(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    seam: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    extruded_circle(
        m,
        SketchPlane::from_axes(base, seam, axis.cross(seam)),
        radius,
        height,
    )
}

fn extruded_circle(m: &mut Model, plane: SketchPlane, radius: f64, height: f64) -> CylinderSolid {
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating the cylinder's base plane: {other:?}"),
    };
    let ring = Ring2d::circle(Point2::from_array([0.0, 0.0]), radius).expect("a stated radius");
    let profile = Profile2d::from_normalized_rings(ring, vec![]);
    match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: height,
        },
    ) {
        Ok(OpOutput::Extrude { solid, faces }) => {
            let [base, top, lateral] = faces[..] else {
                panic!("a circle extrudes to three faces, got {}", faces.len())
            };
            CylinderSolid {
                solid,
                base,
                top,
                lateral,
            }
        }
        other => panic!("extruding the cylinder's circle: {other:?}"),
    }
}

/// An axis-aligned box from `min` to `max` — [`cuboid_on`] with the floor at `min.z` and the
/// height `max.z − min.z`, **exactly**: the height is the difference of the two corners'
/// decimals, and it must be a decimal an `f64` carries (its shortest decimal), or the far cap
/// would not stand at `max.z`. Panics otherwise, saying so — state the floor and the height
/// instead ([`cuboid_on`]); a fixture's input is the test's own.
pub fn cuboid(m: &mut Model, min: Point3, max: Point3) -> Handle<Solid> {
    let [x0, y0, z0] = min.as_array();
    let [x1, y1, z1] = max.as_array();
    // ★ Reversed corners are refused rather than read: a lower `max.z` would become a negative
    // height — a box below the floor — and a reversed `x`/`y` a clockwise rectangle the extrude
    // quietly winds the other way.
    assert!(
        x1 > x0 && y1 > y0 && z1 > z0,
        "cuboid needs max > min on every axis: {min:?} .. {max:?}"
    );
    let dec = |v: f64| Rat::from_decimal(v).expect("a cuboid corner inside the decimal window");
    let diff = dec(z1)
        .checked_sub(dec(z0))
        .expect("a cuboid height inside `Rat`");
    let height = diff.to_f64();
    assert_eq!(
        Rat::from_decimal(height),
        Some(diff),
        "the height {z0} .. {z1} is no decimal an f64 carries — state the floor and the height \
         (`cuboid_on`)"
    );
    cuboid_on(m, [x0, y0], [x1, y1], z0, height)
}

/// An axis-aligned box over the rectangle `xy_min .. xy_max`, standing on the plane `z` and
/// reaching `height` along `+z` — or, for a negative `height`, hanging below the plane. The plane
/// is a datum through `(0, 0, z)`; below it the frame's normal is `−z` (its `ŷ` turned, and the
/// rectangle's `y` with it), so two boxes on one plane in opposite directions **share** it — one
/// statement of `z`, which a stack needs and two heights summed in `f64` do not give.
///
/// The extrude's `[floor cap, far cap, walls…]` order is the fixture's: `faces[0]` is the cap on
/// the plane `z`, `faces[1]` the far one — for [`cuboid`] the bottom and the top. Asserted here,
/// with the box's six faces, twelve edges and eight corners, because tests read the caps by
/// position. The walls' order is the extrude's and nothing reads it.
pub fn cuboid_on(
    m: &mut Model,
    xy_min: [f64; 2],
    xy_max: [f64; 2],
    z: f64,
    height: f64,
) -> Handle<Solid> {
    let ([x0, y0], [x1, y1]) = (xy_min, xy_max);
    assert!(
        x1 > x0 && y1 > y0 && height != 0.0,
        "cuboid_on needs a rectangle and a nonzero height: {xy_min:?} .. {xy_max:?}, {height}"
    );
    let s = if height < 0.0 { -1.0 } else { 1.0 };
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.0, 0.0, z]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, s, 0.0]),
    );
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating the box's floor plane: {other:?}"),
    };
    let p = |x: f64, y: f64| Point2::from_array([x, s * y]);
    let profile = Profile2d::polygon(vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)])
        .expect("a rectangle inside the decimal window");
    let (solid, faces) = match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: height.abs(),
        },
    ) {
        Ok(OpOutput::Extrude { solid, faces }) => (solid, faces),
        other => panic!("extruding the box's rectangle: {other:?}"),
    };
    // What the box is, and where its caps sit — the choices a positional reader leans on.
    let dec = |v: f64| Rat::from_decimal(v).expect("decimal");
    let far = dec(z)
        .checked_add(dec(height))
        .expect("a box height inside `Rat`")
        .to_f64();
    assert_eq!(faces.len(), 6, "a rectangle extrudes to six faces");
    let corners_of = |f: Handle<Face>| -> Vec<[u64; 3]> {
        m.face(f)
            .outer
            .half_edges
            .iter()
            .flat_map(|he| m.edge(he.edge).vertices)
            .map(|v| m.vertex_point(v).as_array().map(f64::to_bits))
            .collect()
    };
    for (f, at) in [(faces[0], z), (faces[1], far)] {
        assert!(
            corners_of(f).iter().all(|c| c[2] == at.to_bits()),
            "faces[0] is the cap on the plane, faces[1] the far one"
        );
    }
    let mut got: Vec<[u64; 3]> = faces.iter().flat_map(|&f| corners_of(f)).collect();
    got.sort_unstable();
    got.dedup();
    let mut want: Vec<[u64; 3]> = [x0, x1]
        .into_iter()
        .flat_map(|x| {
            [y0, y1]
                .into_iter()
                .flat_map(move |y| [z, far].map(|zz| [x, y, zz]))
        })
        .map(|c| c.map(f64::to_bits))
        .collect();
    want.sort_unstable();
    assert_eq!(
        got, want,
        "the box's eight corners are the stated ones, realized"
    );
    let edges: std::collections::HashSet<_> = faces
        .iter()
        .flat_map(|&f| m.face(f).outer.half_edges.iter().map(|he| he.edge))
        .collect();
    assert_eq!(edges.len(), 12, "a box has twelve edges");
    solid
}
