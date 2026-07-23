//! Shared fixtures for the boolean **capability coverage** suite.
//!
//! Every test here drives the kernel through its **public API** only
//! (`apply`/`boolean`/`Operation`, `nacre_topo::Model`, `nacre_props`,
//! `nacre_validate`) — the same front door a user calls. That is the point of
//! this suite: it proves the public path produces the right topology, not that
//! an internal shortcut does. Finer-grained checks of private helpers live in
//! each module's own `#[cfg(test)] mod tests`.
#![allow(dead_code)] // shared helpers: any one coverage module uses only some.

use nacre_geom::intersect::planes_coplanar;
use nacre_geom::{Plane, Surface};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{
    BoolError, BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean, replay,
};
use nacre_scalar::{Axis, Isometry};
use nacre_store::Handle;
use nacre_topo::{Model, Orientation, Solid};

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
    Profile2d {
        points: vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, 1.0), p2(0.0, 1.0)],
    }
}

pub fn regular_ngon(n: usize, r: f64) -> Profile2d {
    let points = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            p2(r * a.cos(), r * a.sin())
        })
        .collect();
    Profile2d { points }
}

pub fn extrude_op(profile: Profile2d, dist: f64) -> Operation {
    Operation::Extrude {
        plane: SketchPlane::world_xy(),
        profile,
        dist,
    }
}

/// Is there an outer-shell face on the plane through `pt` with normal `n`,
/// oriented that way? A capability test that wants to say "a face sits on z = 1.5
/// facing +z" has no face handle in hand — asserting geometry from coordinates is
/// exactly what an acceptance test may do (public `Model`/`Plane` only).
pub fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
    let Some(target) = Plane::from_point_normal(pt, n) else {
        return false;
    };
    let shell = m.solids.get(solid).outer;
    m.shells.get(shell).faces.iter().any(|&fh| {
        let f = m.faces.get(fh);
        let Surface::Plane(plane) = m.surfaces.get(f.surface) else {
            return false;
        };
        let sign = match f.orientation {
            Orientation::Forward => 1.0,
            Orientation::Reversed => -1.0,
        };
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
    let OpOutput::Transform { solid } = apply(m, &Operation::Transform { solid: s, isometry: iso })
        .unwrap()
    else {
        panic!("expected Transform output");
    };
    m.rebuild_adjacency();
    solid
}

/// 30° about Z through (1,1,0) — the standard oblique tilt for rotation-invariance.
pub fn rot30() -> Isometry {
    use nacre_scalar::{Angle, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis: Axis::Z,
        point: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// `deg`° about `axis` through the origin.
pub fn rot_iso(axis: Axis, deg: i128) -> Isometry {
    use nacre_scalar::{Angle, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis,
        point: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
    })
}

// ---- L-prism shape family (public: replay(Extrude) + add_cuboid) ----

/// The canonical L-prism (footprint area 3, height 1 ⇒ volume 3).
pub fn l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d {
        points: vec![
            p2(0.0, 0.0),
            p2(2.0, 0.0),
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
        ],
    };
    let m = replay(&[extrude_op(l, 1.0)]).unwrap();
    let s = m.live_solids[0];
    (m, s)
}

/// The L-prism with its profile wound the other way (vertex order rotated) — the
/// same solid, used to prove a profile's winding cannot change the boolean.
pub fn rotated_l_prism() -> (Model, Handle<Solid>) {
    let l = Profile2d {
        points: vec![
            p2(2.0, 1.0),
            p2(1.0, 1.0),
            p2(1.0, 2.0),
            p2(0.0, 2.0),
            p2(0.0, 0.0),
            p2(2.0, 0.0),
        ],
    };
    let m = replay(&[extrude_op(l, 1.0)]).unwrap();
    let s = m.live_solids[0];
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
    let sh = m.solids.get(s).outer;
    for &fh in &m.shells.get(sh).faces {
        for he in &m.faces.get(fh).outer.half_edges {
            if let Some(bd) = m.edges.get(he.edge).bounds {
                for vh in bd {
                    if seen.insert(vh) {
                        pts.push(m.vertices.get(vh).point.as_array());
                    }
                }
            }
        }
    }
    pts
}
