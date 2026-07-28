//! **A fin array — three walls landing on what the previous two built.**
//!
//! A hub with bars fused around it is the shape a code-CAD user writes in a loop, and it is where
//! the rotated arrangement used to run out of cases. Three fins reproduce it: fuse a diametric
//! pair, then a third at an angle. A sweep of the pair's angle used to fail 8 times in 90, under
//! five different reject names.
//!
//! All five were one root: a boolean's *result* carried no rotation provenance, so its faces were
//! described by their rounded coordinates and one wall became two plane classes on the next fuse.
//! Surfaces now state their own provenance (`SurfaceDef`), so a chained boolean keeps its faces
//! exact and the whole population builds. This file is the population, kept as the regression.

use crate::common::*;
use nacre_math::Point3;
use nacre_ops::{BoolError, BoolKind, boolean};
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

/// **The three arrangements that used to decline, each by a different name.**
///
/// All three are the same shape of input — a diametric pair, then a third fin — and all three
/// failed on that third fuse, under `TraceDeclined{RunSplit}`, `CoincidentNodes` and
/// `UnreachedCell`. Three symptoms, one cause; they build now, and their volumes agree with the
/// answer a (throwaway) forced plane merge produced independently, bit for bit.
#[test]
fn three_fins_build() {
    // Each fin adds `2.4` to the hub's `12` where it does not overlap its neighbours; the pair at
    // (θ, θ+180) is disjoint, the third at +18° shares a wedge near the hub — hence just under
    // the 19.2 upper bound, smoothly in θ.
    for (angles, want) in [
        ([54.0, 234.0, 252.0], 18.773_345_671_087_59),
        ([16.0, 196.0, 214.0], 18.960_779_720_063_556),
        ([20.0, 200.0, 218.0], 18.877_847_194_779_46),
    ] {
        let v = fuse_fins(&angles).unwrap_or_else(|e| panic!("{angles:?}: {e:?}"));
        assert!((v - want).abs() < 1e-12, "{angles:?}: {v} != {want}");
    }
    // …and the same arrangement with exact (90°-family) rotations still builds, unchanged.
    assert!(fuse_fins(&[0.0, 180.0, 198.0]).is_ok());
}

/// **The census: how much of a smooth sweep fails, and of what.**
///
/// `#[ignore]`: 90 booleans of three fuses each. This population — not "the script builds" — is
/// what the fix was judged by: 8 rejects under five names before, none after.
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
    assert_eq!(total, 0, "the failing population changed: {census:?}");
}

/// **A rotation followed by a translation is still refused — at the surface now.**
///
/// The rotation forest records rotations only, so `R` then `T` has no node to name and the plane
/// cannot be stated exactly. That used to be caught on the *vertices* (`TranslateInterleaved`,
/// surfacing as `RotatedUnderdetermined`); the exact plane now comes from the surface, so the
/// refusal has to live there too — `SurfaceDef::Inexact` → `InexactSurface`. Deleting the vertex
/// path must not have quietly turned an honest reject into a rounded answer, which is what this
/// pins.
#[test]
fn rotate_then_translate_is_still_refused() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    m.rebuild_adjacency();
    let a = xf(
        &mut m,
        a,
        Isometry::rotation(Rotation {
            axis: Axis::Z,
            point: [Rat::from_int(0); 3],
            angle: Angle::from_deg(Rat::from_int(30)).expect("angle"),
        }),
    );
    m.rebuild_adjacency();
    let a = xf(
        &mut m,
        a,
        Isometry::translation([Rat::from_int(5), Rat::from_int(-3), Rat::from_int(2)]),
    );
    m.rebuild_adjacency();
    let b = m.add_cuboid(
        Point3::from_array([5.0, -3.0, 3.0]),
        Point3::from_array([8.0, 0.0, 5.0]),
    );
    m.rebuild_adjacency();
    assert_eq!(
        boolean(&mut m, BoolKind::Fuse, a, b).unwrap_err(),
        BoolError::Unsupported {
            reason: nacre_ops::RejectReason::InexactSurface
        }
    );
}

/// **The shape that started this: a hub, twenty fins, two bores.**
///
/// This is the playground script a user actually wrote — a loop of `fuse`, then `cut` — and it
/// died on the fourteenth fin. Twenty fins is not twenty chances to hit one unlucky angle; it is
/// twenty *chained* booleans, each one reasoning about the faces the previous nineteen produced.
/// That is precisely what a result with no provenance could not survive.
///
/// `#[ignore]`: 22 booleans over a solid that grows to ~90 faces.
#[test]
#[ignore = "the full 22-boolean shape (run with --ignored)"]
fn the_whole_fin_array_with_bores_builds() {
    let mut m = Model::new();
    let mut part = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    m.rebuild_adjacency();
    for i in 0..20 {
        let fin = m.add_cuboid(
            Point3::from_array([0.5, -0.2, 1.0]),
            Point3::from_array([4.0, 0.2, 3.0]),
        );
        m.rebuild_adjacency();
        let fin = xf(
            &mut m,
            fin,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::new(i * 18, 1).unwrap()).expect("angle"),
            }),
        );
        m.rebuild_adjacency();
        match boolean(&mut m, BoolKind::Fuse, part, fin) {
            Ok(v) => part = v[0],
            Err(e) => panic!("fin {i} at {}°: {e:?}", i * 18),
        }
        m.rebuild_adjacency();
    }
    let after_fins = nacre_props::mass_props(&m, part).expect("props").volume;
    assert!(
        (after_fins - 57.918_146_700_212_01).abs() < 1e-9,
        "hub + 20 fins: {after_fins}"
    );
    for (n, (a, b)) in [
        ([-0.5, -0.5, -1.0], [0.5, 0.5, 4.0]),
        ([-3.5, -0.1, 1.5], [3.5, 0.1, 2.5]),
    ]
    .into_iter()
    .enumerate()
    {
        let tool = m.add_cuboid(Point3::from_array(a), Point3::from_array(b));
        m.rebuild_adjacency();
        match boolean(&mut m, BoolKind::Cut, part, tool) {
            Ok(v) => part = v[0],
            Err(e) => panic!("bore {n}: {e:?}"),
        }
        m.rebuild_adjacency();
    }
    let final_volume = nacre_props::mass_props(&m, part).expect("props").volume;
    assert!(
        (final_volume - 53.718_146_700_212).abs() < 1e-9,
        "the finished part: {final_volume}"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}
