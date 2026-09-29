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
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::DatumDef;
use nacre_ops::SketchFrame;
use nacre_ops::{
    BoolError, BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply, boolean, replay,
};
use nacre_store::Handle;
use nacre_topo::{Face, Model, Solid, Vertex};

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

/// **A tilted sketch frame that is exact in decimals** — the Pythagorean axes `u = (0.6, 0.8, 0)`,
/// `v = (−0.48, 0.36, 0.8)` (normal `(0.64, −0.48, 0.6)`) at `origin`. Every axis is a decimal
/// the kernel lifts to an exact rational, so a solid sketched here is stated in the world — a
/// cylinder on it has a world axis — while its caches round.
pub fn pythagorean_frame(m: &mut Model, origin: Point3) -> SketchFrame {
    datum_frame(
        m,
        SketchPlane::from_axes(
            origin,
            Vector3::from_array([0.6, 0.8, 0.0]),
            Vector3::from_array([-0.48, 0.36, 0.8]),
        ),
    )
}

/// **A line whose true direction has a zero component the caches do not.** On
/// [`pythagorean_frame`] at the origin: a cylinder `A` over the unit circle about the origin,
/// height 2, and a prism `B` over `(−1.2,−1.5)·(3,−1.5)·(3,1.5)·(1.2,1.5)`, height 1. `B`'s
/// slanted wall runs along the sketch direction `(2.4, 3)`, which is `2.4·u + 3·v = (0, 3, 2.4)`
/// in the world — the `x` component is exactly `0` (`2.4·0.6 = 3·0.48`) — and passes through the
/// origin, so it contains `A`'s axis. Where it crosses `A`'s bottom cap, the line pierces the rim
/// twice at points with one `x`, and a cross product of the rounded plane caches leaves that
/// component a few ulps off `0`. The wall halves the cylinder, so `A ∩ B` is a half cylinder of
/// volume `π/2`. Returns `(model, a, b)`.
pub fn a_tilted_bore_and_an_axial_wall() -> (Model, Handle<Solid>, Handle<Solid>) {
    a_bore_and_a_prism(
        |m| pythagorean_frame(m, Point3::from_array([0.0; 3])),
        &[[-1.2, -1.5], [3.0, -1.5], [3.0, 1.5], [1.2, 1.5]],
    )
}

/// **A wall along a tilted cylinder's axis that cuts it off the axis** — the one
/// [`a_tilted_bore_and_an_axial_wall`] is not: its world normal has no exact zero component, so
/// the two rulings it cuts the lateral in are found where no cache coordinate agrees with the
/// truth. On [`pythagorean_frame`] at the origin, `B` is the prism over
/// `(−1,1.1)·(1,−0.4)·(3,−0.4)·(3,3.1)`, height 1: its slanted wall runs `0.28` from `A`'s axis
/// (`0.6x + 0.8y = 0.28`), and its other walls miss the unit circle (the top one runs `1.43` from
/// the axis, the bottom one starts at `x = 1`). `A ∩ B` is the circular segment beyond the wall,
/// area `acos(0.28) − 0.28·√(1 − 0.28²)`, height 1; `B`'s area is `8.5`. Returns `(model, a, b)`.
pub fn a_tilted_bore_and_a_secant_wall() -> (Model, Handle<Solid>, Handle<Solid>) {
    a_bore_and_a_prism(
        |m| pythagorean_frame(m, Point3::from_array([0.0; 3])),
        &[[-1.0, 1.1], [1.0, -0.4], [3.0, -0.4], [3.0, 3.1]],
    )
}

/// **A unit cylinder `A` (height 2) and a prism `B` over `quad` (height 1), both sketched on the
/// frame `frame` states** — the pair the bore-and-wall fixtures share, `quad` and the frame being
/// all that differs. Both stand on the frame's plane, so `B` overlaps `A` over at most height 1.
/// Returns `(model, a, b)`.
pub fn a_bore_and_a_prism(
    frame: impl FnOnce(&mut Model) -> SketchFrame,
    quad: &[[f64; 2]],
) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let frame = frame(&mut m);
    let extrude = |m: &mut Model, profile: Profile2d, dist: f64| {
        let Ok(OpOutput::Extrude { solid, .. }) = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        ) else {
            panic!("extrude on the fixture's frame")
        };
        m.rebuild_adjacency();
        solid
    };
    let disk = nacre_ops::from_paths(vec![
        nacre_ops::Ring2d::circle(p2(0.0, 0.0), 1.0).expect("a unit circle"),
    ])
    .expect("a disk")
    .remove(0);
    let a = extrude(&mut m, disk, 2.0);
    let quad =
        Profile2d::polygon(quad.iter().map(|q| p2(q[0], q[1])).collect()).expect("a quadrilateral");
    let b = extrude(&mut m, quad, 1.0);
    (m, a, b)
}

/// **A corner that lies on a wall in rationals and off it in `f64`.** The box `x ∈ [0, xs]`,
/// `y ∈ [0, ys]`, `z ∈ [−1, 2]` (the common of two slabs, so `x = xs` and `y = ys` are caps with
/// unit normals) and the prism over the triangle `(0,0),(1,0),(1,k)`, `z ∈ [0, 1]`, whose wall
/// `y = k·x` runs through integer points. With `ys = k·xs` in decimals the box's corner edge is on
/// that wall exactly (`3·0.1 = 0.3`), while `3·fl(0.1) ≠ fl(0.3)`. Returns `(model, box, prism)`.
pub fn decimal_coincidence(xs: f64, ys: f64, k: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    fn extrude(m: &mut Model, axis: Axis, profile: Profile2d, dist: f64) -> Handle<Solid> {
        let op = Operation::Extrude {
            frame: SketchFrame::world(m, axis),
            profile,
            dist,
        };
        let Ok(OpOutput::Extrude { solid, .. }) = apply(m, &op) else {
            panic!("extrude along {axis:?}")
        };
        m.rebuild_adjacency();
        solid
    }
    let rect = |lo: [f64; 2], hi: [f64; 2]| {
        Profile2d::polygon(vec![
            p2(lo[0], lo[1]),
            p2(hi[0], lo[1]),
            p2(hi[0], hi[1]),
            p2(lo[0], hi[1]),
        ])
        .unwrap()
    };
    let mut m = Model::new();
    // `x ∈ [0, xs]`: the YZ frame (`u = ŷ`, `v = ẑ`) swept along `+x̂`.
    let a = extrude(&mut m, Axis::X, rect([-1.0, -1.0], [4.0, 2.0]), xs);
    // `y ∈ [0, ys]`: the ZX frame (`u = ẑ`, `v = x̂`) swept along `+ŷ`.
    let b = extrude(&mut m, Axis::Y, rect([-1.0, -1.0], [2.0, 2.0]), ys);
    let bx = boolean_one(&mut m, BoolKind::Common, a, b).expect("two slabs cross in a box");
    m.rebuild_adjacency();
    let tri = Profile2d::polygon(vec![p2(0.0, 0.0), p2(1.0, 0.0), p2(1.0, k)]).unwrap();
    let prism = extrude(&mut m, Axis::Z, tri, 1.0);
    (m, bx, prism)
}

/// A prism over the world-`XY` polygon `pts`, `z ∈ [0, 1]`.
fn prism_z(m: &mut Model, pts: &[[f64; 2]]) -> Handle<Solid> {
    let profile = Profile2d::polygon(pts.iter().map(|&[x, y]| p2(x, y)).collect()).unwrap();
    let op = Operation::Extrude {
        frame: SketchFrame::world(m, Axis::Z),
        profile,
        dist: 1.0,
    };
    let Ok(OpOutput::Extrude { solid, .. }) = apply(m, &op) else {
        panic!("extrude along z")
    };
    m.rebuild_adjacency();
    solid
}

/// **Two different walls whose `f64` plane caches are the same bits.** Prism `B` over
/// `(0.1,0.3)·(1.1,3.3)·(0.1,3.3)` has its wall on `y = 3x`; prism `C` over `(0,−c)·(1,3)·(0,3)`,
/// `c = 5.551115123125783e-17`, has its wall on `y = (3+c)x − c` — a different plane, crossing
/// `B`'s at `x = 1`. A wall's cache is its producer's row with the truth's first point as anchor:
/// both rows round to `(3, −1, 0)`, and the `d` of both is `fl(3·fl(0.1)) − fl(0.3) = 2⁻⁵⁴ =
/// fl(c)`. So `[3, −1, 0, −c]` describes both, bit for bit, and a class merge that read the cache
/// made them one plane. Returns `(model, b, c)`.
pub fn rounded_twin_walls() -> (Model, Handle<Solid>, Handle<Solid>) {
    let c = 5.551115123125783e-17;
    let mut m = Model::new();
    let b = prism_z(&mut m, &[[0.1, 0.3], [1.1, 3.3], [0.1, 3.3]]);
    let cc = prism_z(&mut m, &[[0.0, -c], [1.0, 3.0], [0.0, 3.0]]);
    (m, b, cc)
}

/// **Two different walls whose face corners' `f64` caches lie on one `f64` plane**, on two
/// solids that do not touch. Prism `B` over `(0,0)·(0.1,0.3)·(−1,1)` has its wall on `y = 3x`
/// with corner caches `(0,0)` and `(fl(0.1), fl(0.3))`; `C` is the prism of
/// [`rounded_twin_walls`] cut to the slab `x ∈ [0.2, 0.4]`, whose wall `y = (3+c)x − c` has its
/// corners at `x = 0.2, 0.4`, `c·0.8` and `c·0.6` below `0.6` and `1.2` — which round to `fl(0.6)`
/// and `fl(1.2)`, and `(fl(0.2), fl(0.6)) = 2·(fl(0.1), fl(0.3))`, `(fl(0.4), fl(1.2)) = 4·(…)` lie
/// exactly on `B`'s corner line. An exact `orient3d` over the corner caches calls the two walls one plane.
/// Returns `(model, b, c)`.
pub fn rounded_corner_walls() -> (Model, Handle<Solid>, Handle<Solid>) {
    let c = 5.551115123125783e-17;
    let mut m = Model::new();
    let b = prism_z(&mut m, &[[0.0, 0.0], [0.1, 0.3], [-1.0, 1.0]]);
    let full = prism_z(&mut m, &[[0.0, -c], [1.0, 3.0], [0.0, 3.0]]);
    let slab = prism_z(
        &mut m,
        &[[0.2, -10.0], [0.4, -10.0], [0.4, 10.0], [0.2, 10.0]],
    );
    let cc = boolean_one(&mut m, BoolKind::Common, full, slab).expect("a slab of the prism");
    m.rebuild_adjacency();
    (m, b, cc)
}

/// Is there an outer-shell face on the plane through `pt` with normal `n`, oriented that way?
/// A capability test that wants to say "a face sits on z = 1.5 facing +z" has no face handle in
/// hand, and asserting geometry from what the test wrote is exactly what an acceptance test may do.
///
/// ★ **Asked of the truth.** The written point and normal are lifted exactly and the plane is
/// named from three points on it — the point and two basis crosses of the normal, the
/// construction `SketchPlane::from_origin_normal` states a plane with (no products, so a long
/// decimal still names it) — and compared with each face's world name; the facing is the exact
/// sign of the name's normal against `n`, times the name's sense and the face's orientation.
pub fn has_face_on_plane(m: &Model, solid: Handle<Solid>, pt: Point3, n: Vector3) -> bool {
    use nacre_exact::Rat;
    let lift = |v: [f64; 3]| v.map(|x| Rat::from_decimal(x).expect("a written decimal"));
    let (o, nn) = (lift(pt.as_array()), lift(n.as_array()));
    let zero = Rat::from_int(0);
    let neg = |x: Rat| zero.checked_sub(x).expect("a lifted decimal negates");
    // Two independent directions square to `n`, each a component shuffle of the written decimals.
    let (u, w) = if nn[0] == zero && nn[1] == zero {
        ([nn[2], zero, zero], [zero, nn[2], zero])
    } else {
        ([neg(nn[1]), nn[0], zero], [neg(nn[2]), zero, nn[0]])
    };
    let add = |a: [Rat; 3], b: [Rat; 3]| -> [Rat; 3] {
        core::array::from_fn(|k| a[k].checked_add(b[k]).expect("a small sum"))
    };
    let Some(target) = nacre_exact::plane_name_exact(o, add(o, u), add(o, w)) else {
        return false;
    };
    let shell = m.solid(solid).outer;
    m.shell(shell).faces.iter().any(|&fh| {
        let f = m.face(fh);
        if m.world_plane_name(f.surface).as_ref() != Some(&target) {
            return false;
        }
        let [a, b, c, _] = target.coeff_ints();
        let along =
            i8::from(nacre_exact::normal_sense(&[a, b, c], nn) == nacre_exact::Orient::Positive)
                * 2
                - 1;
        let sense = m
            .world_plane_name_sense(f.surface)
            .expect("a face with a world name has its sense");
        along * sense.sign() * f.orientation.sign() > 0
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

/// A 2 × 2 block `height` tall on the world XY plane, moved up by `lift` — carried into its
/// statements (every one moves exactly), so the top is the world-stated plane `z = height + lift`
/// even where that has no short decimal (`1/3`). [`through_lifted_block`] is the recorded twin.
pub fn lifted_block(m: &mut Model, height: f64, lift: nacre_exact::Rat) -> Handle<Solid> {
    let square = Profile2d::polygon(vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(2.0, 2.0), p2(0.0, 2.0)])
        .expect("a square");
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: SketchFrame::world(m, Axis::Z),
            profile: square,
            dist: height,
        },
    )
    .expect("a block") else {
        unreachable!()
    };
    let zero = nacre_exact::Rat::from_int(0);
    xf(m, solid, Isometry::translation([zero, zero, lift]))
}

/// [`lifted_block`], then reflected in `x = 0` — carried as well, so its faces are world-stated
/// and take the world's right-handed frames. [`through_mirrored_block`] is the recorded twin.
pub fn mirrored_lifted_block(m: &mut Model, height: f64, lift: nacre_exact::Rat) -> Handle<Solid> {
    let solid = lifted_block(m, height, lift);
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid,
            axis: Axis::X,
            offset: nacre_exact::Rat::from_int(0),
        },
    )
    .expect("a mirror") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// **A block whose every motion is recorded** — the product road to a *folding* motion chain.
///
/// A motion is carried into a solid's statements only when every face's statement moves with
/// it, and a `Through` plane is stated by vertex handles, which no transport moves — so a solid
/// with one `Through` face records whatever moves it, and its other faces gain chains of
/// translations, quarter turns and reflections that fold to a world statement. A fresh block
/// never records such a chain: its moves are carried.
///
/// The `Through` face is the base of a prism raised on a datum through three corners of a helper
/// box `[0,4]×[0,3]×[0,1]` — `(4,0,0)`, `(0,3,0)`, `(4,0,1)`, on no face of the box, so the datum
/// is pushed as a statement of handles rather than interning onto a face's `Known` plane — which
/// is `3x + 4y = 12`. Its normal is a Pythagorean `(3,4,0)/5`, so the datum's frame is rational
/// (`û = ẑ × n̂` horizontal, `v̂ = ẑ`) and the prism is built in the world. The profile is
/// `[0,2] × [0,height]` in that frame, swept `2` along the normal: the top is the horizontal
/// plane `z = height`, `2 × 2` like [`lifted_block`]'s, and the volume is `4·height`.
///
/// ⚠ The helper box stays live — the datum's statement names its vertices — so an assertion over
/// every live vertex or solid sees it too.
pub fn through_block(m: &mut Model, height: f64) -> Handle<Solid> {
    let footprint =
        Profile2d::polygon(vec![p2(0.0, 0.0), p2(4.0, 0.0), p2(4.0, 3.0), p2(0.0, 3.0)])
            .expect("a rectangle");
    let OpOutput::Extrude { solid: helper, .. } =
        apply(m, &extrude_op(m, footprint, 1.0)).expect("the helper box")
    else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let corner = |m: &Model, at: [f64; 3]| {
        m.shell(m.solid(helper).outer)
            .faces
            .iter()
            .flat_map(|&f| m.face(f).outer.half_edges.clone())
            .map(|he| m.he_start(he))
            .find(|&v| m.vertex_point(v).as_array() == at)
            .expect("a corner of the helper box")
    };
    let corners = [[4.0, 0.0, 0.0], [0.0, 3.0, 0.0], [4.0, 0.0, 1.0]].map(|at| corner(m, at));
    let OpOutput::DatumPlane { plane, frame } = apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices(corners),
        },
    )
    .expect("a datum through the helper's corners") else {
        unreachable!()
    };
    assert!(
        matches!(
            m.surface(plane),
            nacre_topo::Surface::Plane {
                points: nacre_topo::PlanePoints::Through(_),
                ..
            }
        ),
        "the datum must be a statement of handles, or nothing here records"
    );
    let side = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(2.0, 0.0),
        p2(2.0, height),
        p2(0.0, height),
    ])
    .expect("a rectangle");
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile: side,
            dist: 2.0,
        },
    )
    .expect("a prism on the datum") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        m.shell(m.solid(solid).outer)
            .faces
            .iter()
            .any(|&f| m.face(f).surface == plane),
        "the prism's base must intern onto the datum's statement"
    );
    let top = top_surface(m, solid);
    assert!(
        m.plane_motion(top).is_none(),
        "the prism is built in the world: its top has no chain before it moves"
    );
    solid
}

/// [`through_block`] moved up by `lift` — recorded whatever `lift` is, so the top carries a
/// translation that folds to `z = height + lift`.
pub fn through_lifted_block(m: &mut Model, height: f64, lift: nacre_exact::Rat) -> Handle<Solid> {
    let solid = through_block(m, height);
    let zero = nacre_exact::Rat::from_int(0);
    let solid = xf(m, solid, Isometry::translation([zero, zero, lift]));
    let top = top_surface(m, solid);
    assert!(
        m.plane_motion(top).is_some() && m.world_plane_name(top).is_some(),
        "the lifted top carries a chain that folds"
    );
    solid
}

/// [`through_lifted_block`], then reflected in `x = 0` — its top carries a folding chain with a
/// reflection in it, the population whose frame is carried out left-handed.
pub fn through_mirrored_block(m: &mut Model, height: f64, lift: nacre_exact::Rat) -> Handle<Solid> {
    let solid = through_lifted_block(m, height, lift);
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid,
            axis: Axis::X,
            offset: nacre_exact::Rat::from_int(0),
        },
    )
    .expect("a mirror") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let top = top_surface(m, solid);
    assert!(
        m.plane_motion(top).is_some() && m.world_plane_name(top).is_some(),
        "the mirrored top carries a chain that folds"
    );
    solid
}

/// The top cap of `solid` — the face highest in `z`.
pub fn top_face(m: &Model, solid: Handle<Solid>) -> Handle<Face> {
    m.shell(m.solid(solid).outer)
        .faces
        .iter()
        .copied()
        .max_by(|&a, &b| {
            let z = |f| {
                nacre_props::face_props(m, f)
                    .expect("props")
                    .centroid
                    .as_array()[2]
            };
            z(a).total_cmp(&z(b))
        })
        .expect("a solid has faces")
}

/// The surface of [`top_face`].
pub fn top_surface(m: &Model, solid: Handle<Solid>) -> Handle<nacre_topo::Surface> {
    m.face(top_face(m, solid)).surface
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

/// Every live vertex of `m`, once each, in walk order (solids, then shells, faces, loops).
pub fn live_vertices(m: &Model) -> Vec<Handle<Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &s in m.live_solids() {
        let sol = m.solid(s);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let f = m.face(fh);
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for &he in &lp.half_edges {
                        let vh = m.he_start(he);
                        if seen.insert(vh) {
                            out.push(vh);
                        }
                    }
                }
            }
        }
    }
    out
}
