//! **Placing a part: what an exact rational move loses on the way in.**
//!
//! A motion's translation is an exact rational (`Isometry::translate: [Rat; 3]`), but the kernel
//! spends it in f64 — `to_f64(t)` and then an f64 add, two roundings — and then treats the result
//! as the truth. These tests pin what that costs, as it stands today:
//!
//! - two parts placed so they share a wall exactly can end up with the wall **1 ULP apart**, so a
//!   `Fuse` returns two disjoint bodies where one is right (and a caller taking `[0]` silently
//!   keeps half the part);
//! - a part that is **rotated and then placed** cannot be described exactly at all, so every
//!   boolean on it is refused.
//!
//! Both assertions describe today's behaviour and are expected to flip.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, RejectReason, boolean};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
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

/// **Two bodies where there is one — the shared wall lands 1 ULP apart.**
///
/// `x = 1` moved by `7/11` and `x = 0` moved by `18/11` are the same real number. But
/// `1.0 + to_f64(7/11)` and `to_f64(18/11)` differ in the last place, so the wall becomes two
/// planes and nothing merges across it. OCCT builds one body of volume 2 from the same operands
/// (`nacre-oracle`); this asserts what we build instead.
#[test]
fn a_shared_wall_one_ulp_apart_splits_the_part() {
    for (n, d) in [(7i128, 11i128), (13, 23)] {
        let vols = place_and_fuse(n, d).expect("the fuse itself succeeds");
        assert_eq!(
            vols.len(),
            2,
            "{n}/{d}: expected today's split, got {vols:?}"
        );
        for v in &vols {
            assert!((v - 1.0).abs() < 1e-9, "{n}/{d}: half a part each, got {v}");
        }
    }
}

/// **Offsets whose wall survives the two roundings still build one body.**
///
/// Pinned beside the failure so the next reader can see the split is not "non-dyadic offsets are
/// unsupported" — most of them are fine. It is the ones where the two roundings disagree.
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

/// **Rotate, then place — refused, every time.**
///
/// The rotation history records rotations only, so `R` then `T` has no node to name and the plane
/// cannot be stated exactly. "Make a feature, turn it, put it where it goes" is ordinary
/// modelling, and none of it builds.
#[test]
fn rotate_then_place_is_refused() {
    let mut refused = 0;
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
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                }),
            );
            let place = Isometry::translation(off.map(Rat::from_int));
            let tool = xf(&mut m, tool, place);
            let base = xf(&mut m, base, place);
            assert_eq!(
                boolean(&mut m, BoolKind::Cut, base, tool).unwrap_err(),
                BoolError::Unsupported {
                    reason: RejectReason::InexactSurface
                },
                "{deg}° then {off:?}"
            );
            refused += 1;
        }
    }
    assert_eq!(refused, 15, "the whole sweep is refused today");
}

/// **A mirror has the same disease — reached by a second route.**
///
/// `AxisMirror` reflects as `2·offset − x` with `offset` already dropped to f64. The doubling is
/// exact (a power of two), so the split needs a *second, independent* route to the same plane:
/// here one wall arrives by reflection and the other by translation. Two roundings each, taken in
/// a different order, and they disagree in the last place — the same failure as
/// [`a_shared_wall_one_ulp_apart_splits_the_part`], through the mirror.
///
/// This is what puts the mirror path in this cell's scope rather than a backlog.
#[test]
fn a_mirrored_wall_and_a_placed_wall_split_the_part() {
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
    assert_eq!(vols.len(), 2, "expected today's split, got {vols:?}");
}
