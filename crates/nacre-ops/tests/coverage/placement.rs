//! **Placing a part: an exact rational move and the f64 images it leaves.**
//!
//! A motion's translation is an exact rational (`Isometry::translate: [Rat; 3]`), but the cache
//! spends it in f64 — `to_f64(t)` and then an f64 add, two roundings. These tests pin that the
//! model does not take that result for the truth:
//!
//! - two parts placed so they share a wall exactly are **one body**, even where the wall's two
//!   f64 images are 1 ULP apart (a `Fuse` returning two disjoint bodies there would let a caller
//!   taking `[0]` silently keep half the part);
//! - a part that is **rotated and then placed** is described exactly, so a boolean on it builds;
//! - decimal dimensions whose f64 images break a rational coincidence (`3·0.1 = 0.3`, but
//!   `3·fl(0.1) ≠ fl(0.3)`) are **judged on the truth**, so a corner that lies on a wall is on it;
//! - two different walls whose rounded images coincide stay **two planes**, each result face on
//!   its own;
//! - a line whose true direction has a `0` component — which its plane caches leave a few ulps
//!   off — is **ordered on the truth**, so two points on it are not read as one;
//! - a wall whose world normal has no exact zero component is **arranged on its rulings** like one
//!   that has.

use crate::common::*;
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, boolean};
use nacre_topo::Model;

/// Two unit cubes placed at `n/d` and `n/d + 1`, so they share the wall `x = (n+d)/d` exactly.
/// Returns the volumes of the solids the `Fuse` produced, and whether the first cube's `+x` wall
/// and the second's `−x` wall came out as **one surface handle** before the fuse.
fn place_and_fuse(n: i128, d: i128) -> Result<(Vec<f64>, bool), BoolError> {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    let b = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([1.0; 3]));
    m.rebuild_adjacency();
    let shift = |m: &mut Model, s, x: Rat| {
        let s = xf(
            m,
            s,
            Isometry::translation([x, Rat::from_int(0), Rat::from_int(0)]),
        );
        m.rebuild_adjacency();
        s
    };
    let a = shift(&mut m, a, Rat::new(n, d).expect("offset"));
    let b = shift(&mut m, b, Rat::new(n + d, d).expect("offset"));
    // The wall on each side: the face whose centroid is furthest along `±x`.
    let wall = |m: &Model, s, sign: f64| {
        m.shell(m.solid(s).outer)
            .faces
            .iter()
            .copied()
            .max_by(|&f, &g| {
                let x = |f| {
                    sign * nacre_props::face_props(m, f)
                        .expect("props")
                        .centroid
                        .as_array()[0]
                };
                x(f).total_cmp(&x(g))
            })
            .map(|f| m.face(f).surface)
            .expect("a cube has faces")
    };
    let one_handle = wall(&m, a, 1.0) == wall(&m, b, -1.0);
    let out = boolean(&mut m, BoolKind::Fuse, a, b)?;
    m.rebuild_adjacency();
    Ok((out.into_iter().map(|s| volume(&m, s)).collect(), one_handle))
}

/// **One body, and one plane, though the two walls' f64 coordinates differ.**
///
/// `x = 1` moved by `7/11` and `x = 0` moved by `18/11` are the same real number, and their f64
/// images differ in the last place. Both moves are carried into the statements, so each wall is
/// stated in the world and the two statements name the same plane — they intern onto **one**
/// handle before the boolean ever asks. The cache is not the thing to fix.
#[test]
fn a_shared_wall_interns_one_plane_though_its_f64_images_differ() {
    for (n, d) in [(7i128, 11i128), (13, 23)] {
        let (vols, one_handle) = place_and_fuse(n, d).expect("fuse");
        assert!(one_handle, "{n}/{d}: the shared wall is one surface");
        assert_eq!(vols.len(), 1, "{n}/{d}: one part, got {vols:?}");
        assert!((vols[0] - 2.0).abs() < 1e-9, "{n}/{d}: {}", vols[0]);
    }
}

/// **Every placement's wall is one plane, and the part one body.**
///
/// The control for the test above, over offsets whose two roundings agree and ones where they
/// disagree: the statements do not round, so the answer does not depend on which.
#[test]
fn most_placements_still_build_one_body() {
    for (n, d) in [
        (1i128, 2i128),
        (1, 3),
        (2, 5),
        (3, 10),
        (1, 7),
        (5, 13),
        (9, 17),
        (11, 19),
        (1, 6),
        (5, 9),
    ] {
        let (vols, one_handle) = place_and_fuse(n, d).expect("fuse");
        assert!(one_handle, "{n}/{d}: the shared wall is one surface");
        assert_eq!(vols.len(), 1, "{n}/{d}: {vols:?}");
        assert!((vols[0] - 2.0).abs() < 1e-9, "{n}/{d}: {}", vols[0]);
    }
}

/// **Rotate, then place — builds, every time.**
///
/// The chain names both motions, in order. A history recording rotations only has no node to
/// name for `R` then `T`, and "make a feature, turn it, put it where it goes" — ordinary
/// modelling — is then refused 15/15 (measured).
#[test]
fn rotate_then_place_builds() {
    let mut built = 0;
    for deg in [7i128, 17, 30, 45, 63] {
        for off in [[3i128, 0, 0], [5, -3, 2], [0, 4, 0]] {
            let mut m = Model::new();
            let base = m.add_cuboid(
                Point3::from_array([-2.0, -2.0, 0.0]),
                Point3::from_array([2.0, 2.0, 1.0]),
            );
            let tool = m.add_cuboid(
                Point3::from_array([-0.5, -0.5, -1.0]),
                Point3::from_array([0.5, 0.5, 2.0]),
            );
            m.rebuild_adjacency();
            let tool = xf(
                &mut m,
                tool,
                Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                }),
            );
            let place = Isometry::translation(off.map(Rat::from_int));
            let tool = xf(&mut m, tool, place);
            let base = xf(&mut m, base, place);
            let out = boolean(&mut m, BoolKind::Cut, base, tool)
                .unwrap_or_else(|e| panic!("{deg}° then {off:?}: {e:?}"));
            m.rebuild_adjacency();
            // The tool is a through-cut of a 4×4×1 slab by a 1×1 post: 16 − 1, whatever the angle.
            let v = volume(&m, out[0]);
            assert!((v - 15.0).abs() < 1e-9, "{deg}° then {off:?}: volume {v}");
            built += 1;
        }
    }
    assert_eq!(built, 15, "the whole sweep builds");
}

/// **A wall reached by reflection and a wall reached by translation are the same wall.**
///
/// `AxisMirror` reflects as `2·offset − x` with `offset` already dropped to `f64`. The doubling is
/// exact (a power of two), so the split needs a *second, independent* route to the same plane:
/// here one wall arrives by reflection and the other by translation. Two roundings each, taken in
/// a different order, and the `f64` images disagree in the last place — the same shape as
/// [`a_shared_wall_interns_one_plane_though_its_f64_images_differ`], through the mirror.
///
/// Both motions are carried into the statements — `2·(1/3) − 1` and `0 − 1/3` are the same
/// rational — so the two walls state one plane and the part is one body. The `f64` coordinates
/// differ in the last place, and they should.
#[test]
fn a_mirrored_wall_and_a_placed_wall_merge_the_part() {
    // Reflect `x = 1` in `x = 1/3`: the image is `−1/3`, computed as `2·fl(1/3) − 1`.
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let nacre_ops::OpOutput::Mirror { solid: a } = nacre_ops::apply(
        &mut m,
        &nacre_ops::Operation::Mirror {
            solid: a,
            axis: Axis::X,
            offset: Rat::new(1, 3).expect("plane"),
        },
    )
    .expect("mirror") else {
        unreachable!("mirror yields Mirror output")
    };
    m.rebuild_adjacency();
    // Reach the same wall by translation instead: `x = 0` moved by `−1/3`.
    let b = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let b = xf(
        &mut m,
        b,
        Isometry::translation([
            Rat::new(-1, 3).expect("x"),
            Rat::from_int(0),
            Rat::from_int(0),
        ]),
    );

    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the fuse itself succeeds");
    m.rebuild_adjacency();
    let vols: Vec<f64> = out.iter().map(|&s| volume(&m, s)).collect();
    assert_eq!(vols.len(), 1, "one body, not a split: got {vols:?}");
    // Two unit cubes meeting on the shared wall.
    assert!((vols[0] - 2.0).abs() < 1e-9, "volume {vols:?}");
}

/// The box's corner edge lies on the prism's wall `y = k·x` in rationals and off it in `f64`, so
/// `common` is the triangular prism over `(0,0),(xs,0),(xs,ys)` — five faces, volume `xs·ys/2`.
/// A judge that answered for the rounded model refused the first case as `ZeroLengthEdge`; the
/// last two are controls whose images keep the coincidence (`2·fl(0.1) = fl(0.2)` is a binary
/// shift, `0.5` and `1.5` are dyadic).
#[test]
fn a_coincidence_the_f64_image_breaks_is_judged_on_the_truth() {
    for (xs, ys, k) in [
        (0.1, 0.3, 3.0),
        (0.7, 2.1, 3.0),
        (0.1, 0.2, 2.0),
        (0.5, 1.5, 3.0),
    ] {
        let case = format!("x = {xs}, y = {ys}, wall y = {k}x");
        let (mut m, bx, prism) = decimal_coincidence(xs, ys, k);
        let v = volume(&m, bx);
        assert!((v - 3.0 * xs * ys).abs() < 1e-12, "{case}: box volume {v}");
        let out = boolean_one(&mut m, BoolKind::Common, bx, prism)
            .unwrap_or_else(|e| panic!("{case}: {e:?}"));
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{case}: {vs:?}");
        let faces = m.shell(m.solid(out).outer).faces.len();
        assert_eq!(faces, 5, "{case}: faces");
        let v = volume(&m, out);
        assert!((v - xs * ys / 2.0).abs() < 1e-12, "{case}: volume {v}");
    }
}

/// The one wall of `s` that is neither axis-aligned nor shared — found by its name, the truth.
fn slanted_wall(
    m: &Model,
    s: nacre_store::Handle<nacre_topo::Solid>,
) -> nacre_store::Handle<nacre_topo::Surface> {
    let walls: Vec<_> = surfaces(m, &[s])
        .into_iter()
        .filter(|&h| {
            let name = m.world_plane_name(h).expect("an unmoved plane is named");
            let ints = name.coeff_ints();
            ints[..3].iter().filter(|c| **c != 0.into()).count() == 2
        })
        .collect();
    assert_eq!(walls.len(), 1, "one slanted wall");
    walls[0]
}

/// Every surface the outer faces of `solids` lie on.
fn surfaces(
    m: &Model,
    solids: &[nacre_store::Handle<nacre_topo::Solid>],
) -> std::collections::BTreeSet<nacre_store::Handle<nacre_topo::Surface>> {
    solids
        .iter()
        .flat_map(|&s| {
            m.shell(m.solid(s).outer)
                .faces
                .iter()
                .map(|&f| m.face(f).surface)
        })
        .collect()
}

/// **Two different walls with one rounded image stay two planes** — each result face lies on the
/// plane of the face it came from, whichever operand comes first.
///
/// * Walls whose face corners' caches lie on one `f64` plane, on two solids that do not touch
///   (`rounded_corner_walls`): the fuse is the two solids, face for face.
/// * Walls whose plane caches are the same bits (`rounded_twin_walls`): the common's slanted face
///   is `B`'s wall `y = 3x` — the binding one over the overlap — and never `C`'s. With the cache
///   admitted as evidence, `Common(C, B)` put it on `C`'s plane (`validate` saw nothing).
#[test]
fn a_wall_keeps_its_plane_when_another_walls_rounded_image_matches() {
    for swap in [false, true] {
        let (mut m, b, c) = rounded_corner_walls();
        let want = surfaces(&m, &[b, c]);
        let (x, y) = if swap { (c, b) } else { (b, c) };
        let out = boolean(&mut m, BoolKind::Fuse, x, y).expect("two apart solids fuse");
        assert_eq!(out.len(), 2, "swap = {swap}: two bodies");
        assert_eq!(
            surfaces(&m, &out),
            want,
            "swap = {swap}: faces on their own planes"
        );
    }
    let (mut m, b, c) = rounded_twin_walls();
    let (b_wall, c_wall) = (slanted_wall(&m, b), slanted_wall(&m, c));
    assert_ne!(b_wall, c_wall);
    let out = boolean_one(&mut m, BoolKind::Common, c, b).expect("common");
    let got = surfaces(&m, &[out]);
    assert!(got.contains(&b_wall), "the binding wall is B's");
    assert!(!got.contains(&c_wall), "C's wall is not on the common");
}

/// **The seam point of two walls whose plane caches are the same bits is where the walls meet.**
/// `Common(B, C)` over [`rounded_twin_walls`] asks for the corner where the two walls cross
/// (`x = 1`). Solved from the three classes' plane caches, two of them read as parallel and the
/// corner was refused (`ThreePlanes`); realized from its definition it is the corner, and the
/// result is the one the other operand order returns — same faces, same volume.
#[test]
fn the_seam_point_of_two_twin_cached_walls_is_realized() {
    let (mut m, b, c) = rounded_twin_walls();
    let bc = boolean_one(&mut m, BoolKind::Common, b, c).expect("the walls meet");
    let (mut m2, b2, c2) = rounded_twin_walls();
    let cb = boolean_one(&mut m2, BoolKind::Common, c2, b2).expect("common");
    assert_eq!(surfaces(&m, &[bc]), surfaces(&m2, &[cb]), "the same faces");
    let (vbc, vcb) = (volume(&m, bc), volume(&m2, cb));
    assert!((vbc - vcb).abs() < 1e-12, "volume {vbc} against {vcb}");
}

/// **Two points on a line are ordered along its true direction**, not along an axis the rounding
/// invents. In [`a_tilted_bore_and_an_axial_wall`] the wall's line across the cylinder's cap runs
/// along `(0, 3, 2.4)`: an axis picked off the plane caches, where that `0` is a few ulps, compares
/// two rim points on `x` — which they share — and every boolean was refused
/// (`TraceDeclined { CoincidentFeatures }`). On the truth each is built, and the volumes are the
/// analytic ones: the wall halves the unit cylinder of height 2, and `B` has area 9 and height 1.
#[test]
fn a_line_is_ordered_along_its_true_direction() {
    use std::f64::consts::PI;
    let (common, a_only, b_only) = (PI / 2.0, 2.0 * PI - PI / 2.0, 9.0 - PI / 2.0);
    for (kind, swapped, want) in [
        (BoolKind::Fuse, false, 2.0 * PI + 9.0 - PI / 2.0),
        (BoolKind::Fuse, true, 2.0 * PI + 9.0 - PI / 2.0),
        (BoolKind::Common, false, common),
        (BoolKind::Common, true, common),
        (BoolKind::Cut, false, a_only),
        (BoolKind::Cut, true, b_only),
    ] {
        let (mut m, a, b) = a_tilted_bore_and_an_axial_wall();
        let (x, y) = if swapped { (b, a) } else { (a, b) };
        let out = boolean_one(&mut m, kind, x, y)
            .unwrap_or_else(|e| panic!("{kind:?} swapped = {swapped}: {e:?}"));
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} swapped = {swapped}: {vs:?}");
        let v = volume(&m, out);
        assert!(
            (v - want).abs() < 1e-9,
            "{kind:?} swapped = {swapped}: volume {v}, want {want}"
        );
    }
}

/// **A wall along a tilted cylinder's axis that cuts it off the axis is arranged on its rulings.**
/// [`a_tilted_bore_and_a_secant_wall`]'s wall has no exact zero in its world normal, so where it
/// meets the lateral no coordinate of any cache agrees with the truth; each boolean builds, and the
/// volumes are the analytic ones — the circular segment beyond the wall, `B`'s area `8.5`, `A`'s
/// `2π`, all at the overlap height 1.
#[test]
fn a_secant_wall_along_a_tilted_axis_is_arranged_on_its_rulings() {
    use std::f64::consts::PI;
    let d: f64 = 0.28;
    let s = d.acos() - d * (1.0 - d * d).sqrt();
    for (kind, swapped, want) in [
        (BoolKind::Fuse, false, 2.0 * PI + 8.5 - s),
        (BoolKind::Fuse, true, 2.0 * PI + 8.5 - s),
        (BoolKind::Common, false, s),
        (BoolKind::Common, true, s),
        (BoolKind::Cut, false, 2.0 * PI - s),
        (BoolKind::Cut, true, 8.5 - s),
    ] {
        let (mut m, a, b) = a_tilted_bore_and_a_secant_wall();
        let (x, y) = if swapped { (b, a) } else { (a, b) };
        let out = boolean_one(&mut m, kind, x, y)
            .unwrap_or_else(|e| panic!("{kind:?} swapped = {swapped}: {e:?}"));
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(vs.is_empty(), "{kind:?} swapped = {swapped}: {vs:?}");
        let v = volume(&m, out);
        assert!(
            (v - want).abs() < 1e-9,
            "{kind:?} swapped = {swapped}: volume {v}, want {want}"
        );
    }
}

/// **Every such wall, over a grid** — the slanted wall of `B` from `(−1, y₀)` to `(1, y₁)`, `y₀` in
/// `0.5…1.6` and `y₁` in `−1.5…0.5` without `0`, on [`pythagorean_frame`]: 240 placements whose
/// wall crosses the lateral off the axis (at small `y₀` the top wall crosses too). Each of the six
/// booleans builds and is valid, and the volumes add up — `Fuse + Common = |A| + |B|`,
/// `(A − B) + Common = |A|`, `(B − A) + Common = |B|` — which a wrong answer `validate` cannot see
/// would break. `y₁ = 0` is left out: that corner lies on the cylinder, the family
/// `edge_on_a_ruling` holds.
#[test]
#[ignore = "1,440 booleans (run with --ignored)"]
fn every_secant_wall_on_a_tilted_bore_builds() {
    let tau = 2.0 * std::f64::consts::PI;
    let mut placements = 0;
    for i0 in 5..=16 {
        for i1 in (-15..=5).filter(|&i| i != 0) {
            let (y0, y1) = (f64::from(i0) / 10.0, f64::from(i1) / 10.0);
            let quad = [[-1.0, y0], [1.0, y1], [3.0, y1], [3.0, y0 + 2.0]];
            // The premise, per placement: the wall's line crosses the unit circle.
            let reach = (y0 + y1).abs() / (4.0 + (y1 - y0) * (y1 - y0)).sqrt();
            assert!(
                reach < 1.0,
                "y₀ = {y0}, y₁ = {y1}: the wall misses the cylinder"
            );
            let area = (0..4)
                .map(|k| {
                    let (p, q) = (quad[k], quad[(k + 1) % 4]);
                    p[0] * q[1] - q[0] * p[1]
                })
                .sum::<f64>()
                .abs()
                / 2.0;
            let run = |kind: BoolKind, swapped: bool| -> f64 {
                let (mut m, a, b) = a_bore_and_a_prism(
                    |m| pythagorean_frame(m, Point3::from_array([0.0; 3])),
                    &quad,
                );
                let (x, y) = if swapped { (b, a) } else { (a, b) };
                let out = boolean(&mut m, kind, x, y).unwrap_or_else(|e| {
                    panic!("y₀ = {y0}, y₁ = {y1}, {kind:?} swapped = {swapped}: {e:?}")
                });
                m.rebuild_adjacency();
                let vs = nacre_validate::validate(&m);
                assert!(
                    vs.is_empty(),
                    "y₀ = {y0}, y₁ = {y1}, {kind:?} swapped = {swapped}: {vs:?}"
                );
                out.iter().map(|&s| volume(&m, s)).sum()
            };
            let common = run(BoolKind::Common, false);
            assert!((run(BoolKind::Common, true) - common).abs() < 1e-9);
            let fuse = run(BoolKind::Fuse, false);
            assert!((run(BoolKind::Fuse, true) - fuse).abs() < 1e-9);
            let (a_only, b_only) = (run(BoolKind::Cut, false), run(BoolKind::Cut, true));
            for (sum, want, what) in [
                (fuse + common, tau + area, "Fuse + Common"),
                (a_only + common, tau, "(A − B) + Common"),
                (b_only + common, area, "(B − A) + Common"),
            ] {
                assert!(
                    (sum - want).abs() < 1e-9,
                    "y₀ = {y0}, y₁ = {y1}: {what} = {sum}, want {want}"
                );
            }
            placements += 1;
        }
    }
    assert_eq!(placements, 240, "the grid");
}
