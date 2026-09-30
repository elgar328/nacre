//! Models taken from user scripts: the rounded plate with bores, the rib, the turned cylinder.

use super::*;

/// The user's rib: a stepped profile standing on the plate's top face, extruded ±20 in y.
fn user_rib(m: &mut Model, x: f64) -> Handle<Solid> {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let line = |a: [f64; 2], b: [f64; 2]| line(p2(a[0], a[1]), p2(b[0], b[1]));
    let profile = stated(vec![
        line([12.0, 0.0], [12.0, 7.5]),
        line([12.0, 7.5], [27.0, 7.5]),
        line([27.0, 7.5], [27.0, 5.5]),
        line([27.0, 5.5], [62.0, 5.5]),
        line([62.0, 5.5], [62.0, -5.5]),
        line([62.0, -5.5], [27.0, -5.5]),
        line([27.0, -5.5], [27.0, -7.5]),
        line([27.0, -7.5], [12.0, -7.5]),
        line([12.0, -7.5], [12.0, 0.0]),
    ])
    .unwrap()
    .remove(0);
    let frame = SketchFrame::world(m, Axis::Y);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 40.0,
        },
    )
    .expect("the rib extrudes") else {
        unreachable!()
    };
    let r = |v: f64| nacre_exact::Rat::from_decimal(v).unwrap();
    transform(
        m,
        solid,
        &nacre_exact::Isometry::translation([r(x), r(-20.0), r(0.0)]),
    )
    .expect("the rib moves")
}

/// The plate's own volume, by hand: a 90 × 50 rectangle less four fillet corners, twelve thick,
/// less four bores.
fn rounded_plate_volume() -> f64 {
    let pi = std::f64::consts::PI;
    (90.0 * 50.0 - 4.0 * (25.0 - pi * 25.0 / 4.0)) * 12.0 - 4.0 * pi * 3.5 * 3.5 * 12.0
}

/// ★ **A fully rounded plate enters a boolean**, in every operation and either operand
/// order. Its cap ring has no three-plane corner at all, so an arm that asks for one against a
/// bore's disk finds no witness and refuses; the engine asks for the witnesses the ring actually
/// has. The partner stands on the plate's top face and overlaps it nowhere, so every
/// answer is arithmetic.
#[test]
fn a_fully_rounded_plate_with_bores_enters_a_boolean() {
    let plate_v = rounded_plate_volume();
    let box_v = 15.0 * 40.0 * 50.0;
    for swapped in [false, true] {
        for (kind, want) in [
            (BoolKind::Fuse, Some(plate_v + box_v)),
            (BoolKind::Cut, Some(if swapped { box_v } else { plate_v })),
            (BoolKind::Common, None),
        ] {
            let mut m = Model::new();
            let plate = rounded_plate(&mut m, 4, 4);
            let boss = m.add_cuboid(
                Point3::from_array([10.0, -20.0, 12.0]),
                Point3::from_array([25.0, 20.0, 62.0]),
            );
            m.rebuild_adjacency();
            let (a, b) = if swapped {
                (boss, plate)
            } else {
                (plate, boss)
            };
            let out = boolean(&mut m, kind, a, b)
                .unwrap_or_else(|e| panic!("{kind:?} swapped {swapped}: {e:?}"));
            m.rebuild_adjacency();
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{kind:?} swapped {swapped}: {:?}",
                nacre_validate::validate(&m)
            );
            match want {
                None => assert!(out.is_empty(), "{kind:?} swapped {swapped}: {out:?}"),
                Some(want) => {
                    assert_eq!(out.len(), 1, "{kind:?} swapped {swapped}: one body");
                    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                    assert!(
                        (v - want).abs() < 1e-9,
                        "{kind:?} swapped {swapped}: {v} vs {want}"
                    );
                }
            }
        }
    }
}

/// ★ **The user's script fuses**: the rounded, bored plate and two stepped ribs standing
/// on it. This is the model that reported the wall (a `witness_not_rational` on the very first
/// fuse, with the cylinder that comes later having nothing to do with it). The ribs stand on the
/// plate and overlap nothing, so the volume is the sum of the parts.
///
/// ★ The `cut` that follows in the script is **not** here: it meets the next wall, a class that
/// carries a circle and rulings at once — measured to be the plate's own side plane `x = 45`,
/// tangent to its corner fillets and crossed by the cylinder. The census row `roundplate user
/// script` holds that hand-off by name.
#[test]
fn the_users_rib_plate_fuses() {
    let mut m = Model::new();
    let plate = rounded_plate(&mut m, 4, 4);
    let r1 = user_rib(&mut m, 17.5);
    m.rebuild_adjacency();
    let ab = boolean(&mut m, BoolKind::Fuse, plate, r1).expect("plate + rib")[0];
    let r2 = user_rib(&mut m, -17.5);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, ab, r2).expect("+ the second rib");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    // Each rib is its profile's area (15 × 15 + 35 × 11) forty long.
    let want = rounded_plate_volume() + 2.0 * (15.0 * 15.0 + 35.0 * 11.0) * 40.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★ **The witnesses a ring offers are a fact about the ring, not about its bores.**
///
/// The supply is read from the corners' own names and solves, so under a rigid motion the rounded
/// plate must fare exactly as it did — and, more to the point, **exactly as the same plate without
/// bores does**: the two sets of refusals are equal.
///
/// ★★ **Both sets are empty.** Every ray from a *filleted* outline's three-plane corners grazes
/// under **six** of the group's rotations, with and without bores alike; giving a component the
/// points its **edges** name (not only its corners) decides those depths, and all six build.
/// **They are checked, not merely accepted** —
/// the `Ok` arm below asserts one body, a clean `validate`, and the exact volume, for every motion.
///
/// So the number is `0`, and this lock says two things at once: «the bores add nothing» and
/// «no motion of a filleted outline is refused».
///
/// (Why no rigid-motion oracle row: see `wall_and_gusset_operand` — a slanted plane's realization
/// differs in the last bit between the two paths, which a bit-exact digest cannot carry.)
#[test]
fn bores_change_nothing_about_a_rounded_plate_under_rigid_motion() {
    let outcome = |bores: usize| -> Vec<String> {
        let mut out = Vec::new();
        for (mn, iso, _) in motion_group() {
            let mut m = Model::new();
            let plate = rounded_plate(&mut m, 4, bores);
            let boss = m.add_cuboid(
                Point3::from_array([10.0, -20.0, 12.0]),
                Point3::from_array([25.0, 20.0, 62.0]),
            );
            let a = transform(&mut m, plate, &iso).expect("the plate moves");
            let b = transform(&mut m, boss, &iso).expect("the box moves");
            m.rebuild_adjacency();
            match boolean(&mut m, BoolKind::Fuse, a, b) {
                Ok(o) => {
                    assert_eq!(o.len(), 1, "{mn}: one body");
                    m.rebuild_adjacency();
                    assert!(
                        nacre_validate::validate(&m).is_empty(),
                        "{mn}: {:?}",
                        nacre_validate::validate(&m)
                    );
                    let v = nacre_props::mass_props(&m, o[0]).expect("props").volume;
                    let want = if bores == 0 {
                        (90.0 * 50.0 - 4.0 * (25.0 - std::f64::consts::PI * 25.0 / 4.0)) * 12.0
                    } else {
                        rounded_plate_volume()
                    } + 15.0 * 40.0 * 50.0;
                    assert!((v - want).abs() < 1e-9, "{mn}: {v} vs {want}");
                }
                Err(e) => out.push(format!("{mn}: {e:?}")),
            }
        }
        out
    };
    let bored = outcome(4);
    let plain = outcome(0);
    assert_eq!(
        bored.len(),
        plain.len(),
        "the bores must change nothing: bored {bored:?} vs plain {plain:?}"
    );
    for (a, b) in bored.iter().zip(plain.iter()) {
        let name = |s: &String| s.split(':').next().unwrap_or("").to_string();
        assert_eq!(name(a), name(b), "the same motions, either way: {a} vs {b}");
        assert!(
            a.contains("NoClearRay"),
            "and the wall that remains is the filleted outline's own, not this cell's: {a}"
        );
    }
    // ★ The population, held as a number so a change in it is loud: **none**.
    assert_eq!(bored.len(), 0, "{bored:?}");
}

/// **The user's script up to its last step** — the rounded, bored plate with both ribs fused on,
/// and the cylinder that bores across them. Two tests read it: the one that names the wall the
/// `cut` stops at, and the one that asks the same question of a moved copy.
fn users_model(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let plate = rounded_plate(m, 4, 4);
    let r1 = user_rib(m, 17.5);
    m.rebuild_adjacency();
    let ab = boolean(m, BoolKind::Fuse, plate, r1).expect("plate + rib")[0];
    let r2 = user_rib(m, -17.5);
    m.rebuild_adjacency();
    let abc = boolean(m, BoolKind::Fuse, ab, r2).expect("+ the second rib")[0];
    let tool = turned_cylinder(m, 10.0, 90.0, [-45.0, 0.0, 47.0]);
    m.rebuild_adjacency();
    (abc, tool)
}

/// ★★★★★ **The user's script builds, `cut` and all.**
///
/// The rounded, bored plate with both ribs fused on, and the cylinder that bores across them. It
/// rests on the footprint reader reading a **disk** — the tool's own cap,
/// which sits on the plate's side plane and clears the corner fillet's tangent line by twice its
/// radius. Unread, that clearance folds to "did not clear", a tangency row
/// is written for a contact that is not there, and the gate refuses the boolean on it.
///
/// The volume is the oracle: the fuse's `100695.2213` less the two rib bores `2·π·10²·11`.
///
/// ★ **The audit's denominator rides along.** «No circle meets a ruling» is worth nothing
/// unless the instrument looked: this asserts it measured pairs, failed to measure none, and found
/// no circle sitting inside a ruling strip — the shape the walk has never been handed.
#[test]
fn the_users_script_builds() {
    let mut m = Model::new();
    let (abc, tool) = users_model(&mut m);
    let out = boolean(&mut m, BoolKind::Cut, abc, tool).expect("the bore across the ribs");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = rounded_plate_volume() + 2.0 * (15.0 * 15.0 + 35.0 * 11.0) * 40.0
        - 2.0 * std::f64::consts::PI * 100.0 * 11.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    #[cfg(debug_assertions)]
    {
        use std::sync::atomic::Ordering;
        let audit = &crate::arrangement::MIXED_CLASS_AUDIT;
        let (pairs, unmeasured, inside) = (
            audit.pairs.load(Ordering::Relaxed),
            audit.unmeasured.load(Ordering::Relaxed),
            audit.inside.load(Ordering::Relaxed),
        );
        assert!(pairs > 0, "the mixed-class audit looked at nothing");
        assert_eq!(
            unmeasured, 0,
            "the audit could not measure {unmeasured} pairs"
        );
        assert_eq!(
            inside, 0,
            "a circle now sits inside a ruling strip — the walk's first time"
        );
    }
}

/// ★★★★ **And it builds the same body wherever it stands.**
///
/// Two newly-admitted paths meet in this model (the mixed class, arranged for real, and the
/// disk read from its own plane's exact statement), and a newly-admitted path is exactly
/// where a frame-dependent slip would hide. So the whole assembly is carried somewhere else by an
/// exact translation and asked again — now with a volume to answer with, not just a wall's name.
#[test]
fn the_users_script_builds_the_same_body_after_a_motion() {
    let mut m = Model::new();
    let (abc, tool) = users_model(&mut m);
    let d = |v: f64| nacre_exact::Rat::from_decimal(v).expect("a short decimal");
    let iso = nacre_exact::Isometry::translation([d(7.5), d(-13.0), d(4.25)]);
    let abc = transform(&mut m, abc, &iso).expect("the body moves");
    let tool = transform(&mut m, tool, &iso).expect("the tool moves");
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Cut, abc, tool).expect("the moved bore");
    assert_eq!(out.len(), 1, "one body");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let want = rounded_plate_volume() + 2.0 * (15.0 * 15.0 + 35.0 * 11.0) * 40.0
        - 2.0 * std::f64::consts::PI * 100.0 * 11.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}
