//! **A fin array — three walls landing on what the previous two built.**
//!
//! A hub with bars fused around it is the shape a code-CAD user writes in a loop, and it is where
//! the rotated arrangement runs out of cases. Three fins reproduce it: fuse a diametric pair, then
//! a third at an angle, and the trace declines. Sweeping the pair's angle shows it is not one
//! unlucky number — a fraction of a smooth sweep fails, for **three different reasons**.
//!
//! The rejects here are asserted as they stand today. Each is a coverage gap with its own cause,
//! and fixing one flips exactly its own assertion — which is the point of naming them apart.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, DeclineKind, RejectReason, boolean};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

/// The hub, and one bar per angle: `[0.5,4]` long, `0.4` wide, spun about Z through the origin.
///
/// The bar reaches *into* the hub (`x` from `0.5`, the hub spans `±1`), so each fin genuinely
/// intersects it and its neighbours rather than merely touching.
fn hub_and_fins(angles: &[f64]) -> (Model, Handle<Solid>, Vec<Handle<Solid>>) {
    let mut m = Model::new();
    let hub = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    m.rebuild_adjacency();
    let fins = angles
        .iter()
        .map(|&a| {
            let f = m.add_cuboid(
                Point3::from_array([0.5, -0.2, 1.0]),
                Point3::from_array([4.0, 0.2, 3.0]),
            );
            m.rebuild_adjacency();
            // Rational degrees: the angle is exact, only its cos/sin are not.
            let deg = Rat::new((a * 100.0).round() as i128, 100).expect("angle");
            xf(
                &mut m,
                f,
                Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    point: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(deg).expect("angle"),
                }),
            )
        })
        .collect();
    (m, hub, fins)
}

/// Fuse the fins onto the hub in order; `Ok(volume)` or the first reject.
fn fuse_fins(angles: &[f64]) -> Result<f64, BoolError> {
    let (mut m, hub, fins) = hub_and_fins(angles);
    let mut part = hub;
    for f in fins {
        part = boolean(&mut m, BoolKind::Fuse, part, f)?[0];
        m.rebuild_adjacency();
    }
    Ok(nacre_props::mass_props(&m, part).expect("props").volume)
}

/// **The three roots, each with its own name.**
///
/// All three are the same shape of input — a diametric pair, then a third fin — and all three fail
/// on that third fuse. They used to be one reject name (`LoopOrientMismatch`, raised from 22
/// places), which said the engine was unhappy without saying about what.
#[test]
fn three_fins_reject_by_cause() {
    // A run — an edge of the face lying *on* the class plane — with the third fin's wall crossing
    // it strictly inside. The trace requires a run's nodes to be adjacent after sorting. (The
    // face handle rides along in the reason and is not asserted: it is an index into an arena
    // whose numbering is not this test's business.)
    let e = fuse_fins(&[54.0, 234.0, 252.0]).unwrap_err();
    assert!(
        matches!(
            e,
            BoolError::Unsupported {
                reason: RejectReason::TraceDeclined {
                    kind: DeclineKind::RunSplit,
                    ..
                }
            }
        ),
        "{e:?}"
    );
    // Two plane triples naming one point — a fourth plane through it, so an edge has no direction.
    assert_eq!(
        fuse_fins(&[16.0, 196.0, 214.0]).unwrap_err(),
        BoolError::Unsupported {
            reason: RejectReason::CoincidentNodes
        }
    );
    // A component of the class's subdivision that the inside/outside labels never reach.
    assert_eq!(
        fuse_fins(&[20.0, 200.0, 218.0]).unwrap_err(),
        BoolError::Unsupported {
            reason: RejectReason::UnreachedCell
        }
    );
    // …and the same arrangement with exact (90°-family) rotations builds, so none of this is
    // about the fin *shape*.
    assert!(fuse_fins(&[0.0, 180.0, 198.0]).is_ok());
}

/// **The census: how much of a smooth sweep fails, and of what.**
///
/// `#[ignore]`: 90 booleans of three fuses each. The number this prints is what a fix is judged
/// by — "the script builds" is one point, this is the population.
#[test]
#[ignore = "slow census sweep (run with --ignored)"]
fn theta_sweep_census() {
    let mut census: std::collections::BTreeMap<String, Vec<i64>> = Default::default();
    for t in 0..90 {
        let th = t as f64;
        if let Err(e) = fuse_fins(&[th, th + 180.0, th + 198.0]) {
            let name = match e {
                BoolError::Unsupported { reason } => reason.as_str().to_string(),
                other => format!("{other:?}"),
            };
            census.entry(name).or_default().push(t);
        }
    }
    let total: usize = census.values().map(|v| v.len()).sum();
    eprintln!("[fin census] {total}/90 reject: {census:?}");
    assert_eq!(total, 8, "the failing population changed: {census:?}");
}
