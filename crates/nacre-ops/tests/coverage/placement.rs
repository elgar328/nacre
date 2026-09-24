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
//!   `3·fl(0.1) ≠ fl(0.3)`) are **judged on the truth**, so a corner that lies on a wall is on it.

use crate::common::*;
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, boolean};
use nacre_topo::Model;

/// Two unit cubes placed at `n/d` and `n/d + 1`, so they share the wall `x = (n+d)/d` exactly.
/// Returns the solids the `Fuse` produced, with their volumes.
fn place_and_fuse(n: i128, d: i128) -> Result<Vec<f64>, BoolError> {
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
    let out = boolean(&mut m, BoolKind::Fuse, a, b)?;
    m.rebuild_adjacency();
    Ok(out.into_iter().map(|s| volume(&m, s)).collect())
}

/// **One body, though the two walls' f64 coordinates differ.**
///
/// `x = 1` moved by `7/11` and `x = 0` moved by `18/11` are the same real number, and their f64
/// images differ in the last place. Each wall carries its motion in its *definition*, so the
/// judgment realizes both at high precision, finds them the same plane, and merges. The cache is
/// not the thing to fix.
#[test]
fn a_shared_wall_merges_though_its_f64_images_differ() {
    for (n, d) in [(7i128, 11i128), (13, 23)] {
        let vols = place_and_fuse(n, d).expect("fuse");
        assert_eq!(vols.len(), 1, "{n}/{d}: one part, got {vols:?}");
        assert!((vols[0] - 2.0).abs() < 1e-9, "{n}/{d}: {}", vols[0]);
    }
}

/// **Offsets whose wall survives the two roundings still build one body.**
///
/// The control for the test above: most non-dyadic offsets round alike on both walls; the ones
/// where the two roundings disagree are the ones that test.
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
        let vols = place_and_fuse(n, d).expect("fuse");
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
/// exact (a power of two), so the split needed a *second, independent* route to the same plane:
/// here one wall arrives by reflection and the other by translation. Two roundings each, taken in
/// a different order, and the `f64` images disagree in the last place — the same shape as
/// [`a_shared_wall_merges_though_its_f64_images_differ`], through the mirror.
///
/// The chain holds reflections — `1/3` is not dyadic, so the node is recorded — and the two
/// definitions realize to the same plane. (Carried instead by *conjugating* an existing chain, a
/// surface with no chain has nothing to conjugate, its mirror image is declared exact, and the
/// part fuses into **two** bodies.) The `f64` coordinates differ in the last place, and they
/// should.
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
