//! Shared fixtures for the boolean **capability coverage** suite.
//!
//! Every test here drives the kernel through its **public API** only
//! (`apply`/`boolean`/`Operation`, `nacre_topo::Model`, `nacre_props`,
//! `nacre_validate`) — the same front door a user calls. That is the point of
//! this suite: it proves the public path produces the right topology, not that
//! an internal shortcut does. Finer-grained checks of private helpers live in
//! each module's own `#[cfg(test)] mod tests`.
#![allow(dead_code)] // shared helpers: any one coverage module uses only some.

use nacre_exact::{Axis, Isometry};
use nacre_geom::intersect::planes_coplanar;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::DatumDef;
use nacre_ops::SketchFrame;
use nacre_ops::{
    BoolError, BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean, replay,
};
use nacre_store::Handle;
use nacre_topo::{Face, Model, Solid};

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
pub fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

/// The signed volume of a solid, via the public mass-properties crate — the
/// independent oracle every capability test asserts against.
pub fn volume(m: &Model, s: Handle<Solid>) -> f64 {
    nacre_props::mass_props(m, s).unwrap().volume
}

/// A boolean whose result is exactly one solid (the common single-body case),
/// over the public `boolean`. Asserts the single-solid shape so call sites read
/// as before.
pub fn boolean_one(
    model: &mut Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<Handle<Solid>, BoolError> {
    let solids = boolean(model, kind, a, b)?;
    assert_eq!(
        solids.len(),
        1,
        "boolean_one: expected one solid, got {}",
        solids.len()
    );
    Ok(solids[0])
}

pub fn p2(x: f64, y: f64) -> Point2 {
    Point2::from_array([x, y])
}

pub fn square() -> Profile2d {
    Profile2d::polygon(vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)]).unwrap()
}

pub fn regular_ngon(n: usize, r: f64) -> Profile2d {
    let points = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            p2(r * a.cos(), r * a.sin())
        })
        .collect();
    Profile2d::polygon(points).unwrap()
}

/// A world-XY extrude **for `m`** — the frame names that model's seeded plane.
///
/// ★ It takes the model because an extrude now *names* its plane rather than carrying it, and a
/// handle is only valid in its own arena. For a log that will be `replay`ed rather than applied,
/// see [`extrude_log_op`].
pub fn extrude_op(m: &Model, profile: Profile2d, dist: f64) -> Operation {
    Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile,
        dist,
    }
}

/// The same, for a log that is **replayed**: the frame names a throwaway model's seed, which
/// `replay` re-anchors onto the model it builds. Applying one of these to a
/// live model is the cross-model misuse the debug guard catches.
pub fn extrude_log_op(profile: Profile2d, dist: f64) -> Operation {
    extrude_op(&Model::new(), profile, dist)
}

/// Is there an outer-shell face on the plane through `pt` with normal `n`,
/// oriented that way? A capability test that wants to say "a face sits on z = 1.5
/// facing +z" has no face handle in hand — asserting geometry from coordinates is
/// exactly what an acceptance test may do (public `Model`/`Plane` only).
pub fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
    let Some(target) = Plane::from_point_normal(pt, n) else {
        return false;
    };
    let shell = m.solid(solid).outer;
    m.shell(shell).faces.iter().any(|&fh| {
        let f = m.face(fh);
        let Surface::Plane(plane) = m.surface_cache(f.surface) else {
            return false;
        };
        let sign = f64::from(f.orientation.sign());
        planes_coplanar(plane, &target) && (plane.normal() * sign).dot(n) > 0.0
    })
}

/// Two unit boxes overlapping in a corner: A = [0,1]³, B = [0.5,1.5]³.
pub fn two_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.5]),
        Point3::from_array([1.5, 1.5, 1.5]),
    );
    (m, a, b)
}

/// Two unit boxes stacked sharing their z=1 interface plane (a face-contact pair,
/// not overlapping volume): A = [0,1]³ below, B = z∈[1,2] above.
pub fn stacked_cubes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([1.0, 1.0, 2.0]),
    );
    (m, a, b)
}

// ---- rigid transforms (public: apply(Operation::Transform)) ----

/// Apply a rigid transform to `s` via the public `apply(Transform)` path and
/// return the image, rebuilding adjacency (as every call site does).
pub fn xf(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    )
    .unwrap() else {
        panic!("expected Transform output");
    };
    m.rebuild_adjacency();
    solid
}

/// 30° about Z through (1,1,0) — the standard oblique tilt for rotation-invariance.
pub fn rot30() -> Isometry {
    use nacre_exact::{Angle, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// `deg`° about `axis` through the origin.
pub fn rot_iso(axis: Axis, deg: i128) -> Isometry {
    use nacre_exact::{Angle, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
    })
}

// ---- L-prism shape family (public: replay(Extrude) + add_cuboid) ----

/// The canonical L-prism (footprint area 3, height 1 ⇒ volume 3).
pub fn l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    let s = m.live_solids()[0];
    (m, s)
}

/// The L-prism with its profile wound the other way (vertex order rotated) — the
/// same solid, used to prove a profile's winding cannot change the boolean.
pub fn rotated_l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d::polygon(vec![
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
        p2(0.0, 0.0),
        p2(2.0, 0.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(l, 1.0)]).unwrap();
    let s = m.live_solids()[0];
    (m, s)
}

pub fn l_and_corner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([1.3, -0.3, 0.2]),
        Point3::from_array([2.4, 0.4, 1.4]),
    );
    (m, l, bx)
}

pub fn l_and_rod() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let rod = m.add_cuboid(
        Point3::from_array([0.3, 0.3, -0.5]),
        Point3::from_array([0.5, 0.6, 1.5]),
    );
    (m, l, rod)
}

pub fn l_and_inner_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(Point3::from_array([0.1; 3]), Point3::from_array([0.9; 3]));
    (m, l, bx)
}

/// The distinct boundary-vertex coordinates of a solid's outer shell.
pub fn outer_points(m: &Model, s: Handle<Solid>) -> Vec<[f64; 3]> {
    let mut seen = std::collections::HashSet::new();
    let mut pts = Vec::new();
    let sh = m.solid(s).outer;
    for &fh in &m.shell(sh).faces {
        for he in &m.face(fh).outer.half_edges {
            {
                for vh in m.edge(he.edge).vertices {
                    if seen.insert(vh) {
                        pts.push(m.vertex_point(vh).as_array());
                    }
                }
            }
        }
    }
    pts
}

/// The L with a box straddling its reflex corner (1,1): a single chord with one
/// reflex bend — the box vertex (1.6,1.6,·) sits in the notch (inside the convex
/// hull, outside the L), where a convex half-space test would misclassify it.
pub fn l_and_reflex_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([0.6, 0.6, 0.2]),
        Point3::from_array([1.6, 1.6, 1.4]),
    );
    (m, l, bx)
}

// ---- U-prism and multi-chord shapes (public: apply(Extrude) + add_cuboid) ----

/// A U-prism (two prongs + a bridge), footprint spanning a reflex-rich outline.
pub fn u_prism() -> (Model, Handle<Solid>) {
    let u = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(3.0, 0.0),
        p2(3.0, 2.3),
        p2(2.0, 2.3),
        p2(2.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 2.0),
        p2(0.0, 2.0),
    ])
    .unwrap();
    let m = replay(&[extrude_log_op(u, 1.0)]).unwrap();
    let s = m.live_solids()[0];
    (m, s)
}

pub fn u_and_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, u) = u_prism();
    let slab = m.add_cuboid(
        Point3::from_array([-0.5, 1.5, -0.5]),
        Point3::from_array([3.5, 2.5, 1.5]),
    );
    (m, u, slab)
}

/// Extrude a profile on a plane offset up the z-axis (a bar sitting above z=0.5).
fn extrude_at_z(m: &mut Model, profile: Profile2d, z: f64, dist: f64) -> Handle<Solid> {
    let __g200 = datum_frame(
        m,
        SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
    );
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g200,
            profile,
            dist,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    solid
}

/// The unit cube, and a bar across it spun `deg`° about Y through `(0.5, ·, pivot_z)`.
///
/// The bar is `0.4` wide in x and `2·half_z` tall, centred on the cube's top plane `z = 1`, so a
/// **square** cross-section (`half_z == 0.2`) spun **45°** sends two corners exactly `half_z·√2`
/// sideways and back to `z = 1`: a tool *edge* lying in the target's face plane, whose endpoints
/// are named by three bar planes while a fourth — the cube's top — passes through them.
///
/// Every argument is a knob on that coincidence, which is what makes the family worth sharing:
/// change the angle, break the squareness, or lift the pivot off the plane and the same model
/// builds. Rotation about Y ignores the pivot's y, so only `pivot_z` matters here.
/// The same family with the bar's far x edge nudged by `ulps` — the non-degenerate neighbours of
/// the four-plane case, which the concurrency needs to be bit-exact to survive.
pub fn cube_and_spun_bar_ulp(
    half_z: f64,
    deg: i128,
    pivot_z: nacre_exact::Rat,
    ulps: i64,
) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut x_hi = 0.7f64;
    for _ in 0..ulps.unsigned_abs() {
        x_hi = if ulps > 0 {
            f64::from_bits(x_hi.to_bits() + 1)
        } else {
            f64::from_bits(x_hi.to_bits() - 1)
        };
    }
    cube_and_spun_bar_x(half_z, deg, pivot_z, x_hi)
}

pub fn cube_and_spun_bar(
    half_z: f64,
    deg: i128,
    pivot_z: nacre_exact::Rat,
) -> (Model, Handle<Solid>, Handle<Solid>) {
    cube_and_spun_bar_x(half_z, deg, pivot_z, 0.7)
}

pub fn cube_and_spun_bar_x(
    half_z: f64,
    deg: i128,
    pivot_z: nacre_exact::Rat,
    x_hi: f64,
) -> (Model, Handle<Solid>, Handle<Solid>) {
    use nacre_exact::{Angle, Rat, Rotation};
    let mut m = Model::new();
    let __op = extrude_op(&m, square(), 1.0);
    let OpOutput::Extrude { solid: cube, .. } = apply(&mut m, &__op).unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    let bar = Profile2d::polygon(vec![
        p2(0.3, -0.5),
        p2(x_hi, -0.5),
        p2(x_hi, 1.5),
        p2(0.3, 1.5),
    ])
    .unwrap();
    let bar = extrude_at_z(&mut m, bar, 1.0 - half_z, 2.0 * half_z);
    m.rebuild_adjacency();
    let bar = xf(
        &mut m,
        bar,
        Isometry::rotation(Rotation {
            axis: Axis::Y,
            pivot: [Rat::new(1, 2).unwrap(), Rat::from_int(0), pivot_z],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        }),
    );
    (m, cube, bar)
}

pub fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bar = Profile2d::polygon(vec![
        p2(1.8, 0.8),
        p2(2.1, 0.8),
        p2(2.1, 2.1),
        p2(0.8, 2.1),
        p2(0.8, 1.8),
        p2(1.8, 1.8),
    ])
    .unwrap();
    let b = extrude_at_z(&mut m, bar, 0.5, 1.0);
    (m, l, b)
}

pub fn l_and_ell_stub() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let ell = Profile2d::polygon(vec![
        p2(0.2, 0.25),
        p2(0.85, 0.25),
        p2(0.85, 0.4),
        p2(0.35, 0.4), // reflex
        p2(0.35, 0.9),
        p2(0.2, 0.9),
    ])
    .unwrap();
    let stub = extrude_at_z(&mut m, ell, 0.5, 1.0);
    (m, l, stub)
}

pub fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let staple = Profile2d::polygon(vec![
        p2(0.1, 0.5),
        p2(0.6, 0.5),
        p2(0.6, 1.3),
        p2(0.8, 1.3),
        p2(0.8, 0.45),
        p2(1.4, 0.45),
        p2(1.4, 1.5),
        p2(0.1, 1.5),
    ])
    .unwrap();
    let __frame0 = datum_frame(
        &mut m,
        SketchPlane::from_origin_normal(
            Point3::from_array([0.0, 1.3, 0.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
        )
        .expect("a unit normal"),
    );
    let OpOutput::Extrude { solid: st, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __frame0,
            profile: staple,
            dist: 0.65,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, st)
}

// ---- pad / pocket features (public: apply(Operation::Pad/PocketOnFace)) ----

/// A unit cube with its top face handle — the standard target for pad/pocket.
pub fn cube_with_top() -> (Model, Handle<Face>) {
    let mut m = Model::new();
    let __op = extrude_op(&m, square(), 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &__op).unwrap() else {
        unreachable!()
    };
    let top = faces[1]; // base, top, sides…
    (m, top)
}

/// The `0.4` boss/pocket footprint on `[0.3, 0.7]²` of the unit cube's lid.
///
/// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
/// onto the face's plane, and the axes are the arbitrary-axis convention's, which for `n = ẑ` are
/// `u = +x̂`, `v = +ŷ` — so a frame point `(a, b)` is world `(a, b, 1)`, the identity.
pub fn small_square() -> Profile2d {
    Profile2d::polygon(vec![p2(0.3, 0.7), p2(0.3, 0.3), p2(0.7, 0.3), p2(0.7, 0.7)]).unwrap()
}

pub fn pad_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
    Operation::PadOnFace {
        face,
        profile,
        dist,
    }
}

pub fn pocket_op(face: Handle<Face>, profile: Profile2d, dist: f64) -> Operation {
    Operation::PocketOnFace {
        face,
        profile,
        dist,
    }
}

// ---- more shape fixtures (all public: add_cuboid / apply(Extrude|Pocket)) ----

pub fn cube_and_notch() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let y = m.add_cuboid(
        Point3::from_array([3.0, -1.0, -1.0]),
        Point3::from_array([7.0, 1.4, 1.2]),
    );
    (m, a, y)
}

pub fn nested_boxes() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([3.0; 3]));
    let b = m.add_cuboid(Point3::from_array([1.0; 3]), Point3::from_array([2.0; 3]));
    (m, a, b)
}

pub fn l_and_popup_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let bx = m.add_cuboid(
        Point3::from_array([0.5, 0.5, 0.2]),
        Point3::from_array([2.5, 1.5, 1.2]),
    );
    (m, l, bx)
}

pub fn l_and_dimple() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, l) = l_prism();
    let stub = m.add_cuboid(
        Point3::from_array([0.3, 0.3, 0.5]),
        Point3::from_array([0.7, 0.7, 1.5]),
    );
    (m, l, stub)
}

/// A unit cube with a small blind pocket carved in its top.
pub fn pocketed_cube() -> (Model, Handle<Solid>) {
    let (mut m, top) = cube_with_top();
    let OpOutput::PocketOnFace { solid, .. } =
        apply(&mut m, &pocket_op(top, small_square(), 0.5)).unwrap()
    else {
        unreachable!()
    };
    (m, solid)
}

pub fn top_pocketed_cube() -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let __op = extrude_op(&m, square(), 1.0);
    let OpOutput::Extrude { faces, .. } = apply(&mut m, &__op).unwrap() else {
        unreachable!()
    };
    let OpOutput::PocketOnFace { solid, .. } =
        apply(&mut m, &pocket_op(faces[1], small_square(), 0.5)).unwrap()
    else {
        unreachable!()
    };
    (m, solid)
}

pub fn pocket_and_slab(z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, pc) = pocketed_cube();
    let slab = m.add_cuboid(
        Point3::from_array([-0.2, -0.25, z0]),
        Point3::from_array([1.3, 1.2, 1.5]),
    );
    (m, slab, pc)
}

/// A slot that runs off **one** edge of the unit cube's lid — world `x ∈ [0.25, 0.75]`,
/// `y ∈ [-0.25, 0.75]`. On a lid, frame coordinates are world coordinates (see [`small_square`]).
pub fn edge_overhang_profile() -> Profile2d {
    Profile2d::polygon(vec![
        p2(0.25, 0.75),
        p2(0.25, -0.25),
        p2(0.75, -0.25),
        p2(0.75, 0.75),
    ])
    .unwrap()
}

/// A channel that runs off **both** opposite edges of the unit cube's lid — world
/// `x ∈ [0.25, 0.75]`, `y ∈ [-0.25, 1.25]`.
pub fn spanning_slab_profile() -> Profile2d {
    Profile2d::polygon(vec![
        p2(0.25, 1.25),
        p2(0.25, -0.25),
        p2(0.75, -0.25),
        p2(0.75, 1.25),
    ])
    .unwrap()
}

/// f64 approximate equality for coordinate/volume comparisons.
pub fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

pub fn test_iso() -> (Isometry, [f64; 3]) {
    use nacre_exact::Rat;
    (
        Isometry::translation([
            Rat::new(7, 2).unwrap(),
            Rat::from_int(-4),
            Rat::from_int(11),
        ]),
        [3.5, -4.0, 11.0],
    )
}

pub fn translate_iso(off: [i128; 3]) -> Isometry {
    use nacre_exact::Rat;
    Isometry::translation([
        Rat::from_int(off[0]),
        Rat::from_int(off[1]),
        Rat::from_int(off[2]),
    ])
}

pub fn rigid_iso(axis: Axis, deg: i128, off: [i128; 3]) -> Isometry {
    use nacre_exact::{Angle, Rat, Rotation};
    Isometry::rigid(
        Rotation {
            axis,
            pivot: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
        },
        [
            Rat::from_int(off[0]),
            Rat::from_int(off[1]),
            Rat::from_int(off[2]),
        ],
    )
}
