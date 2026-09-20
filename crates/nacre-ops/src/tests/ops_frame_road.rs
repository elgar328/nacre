use super::*;
use crate::exact::RatFrame;

/// The frame road's realized basis, wrapped so `exact()` can be asked of it — the same two
/// lines `extrude_on_frame` will run.
fn realized(model: &Model, frame: &SketchFrame) -> Option<SketchPlane> {
    let (o, u, v, _) = crate::rotated_vertex::frame_world_basis(
        model,
        frame.plane(),
        frame.placement(),
        frame.flip(),
    )?;
    Some(realized_plane(
        Point3::from_array(o),
        Vector3::from_array(u),
        Vector3::from_array(v),
    ))
}

/// State `sp` as a datum and hand back both roads' frames for it.
fn both(model: &mut Model, sp: SketchPlane) -> (Option<RatFrame>, Option<RatFrame>) {
    let OpOutput::DatumPlane { frame, .. } = apply(
        model,
        &Operation::DatumPlane {
            def: DatumDef::Stated(sp),
        },
    )
    .expect("every fixture here states a plane the kernel can hold") else {
        unreachable!()
    };
    let by_value = sp.exact();
    // ★ The repaired road: ask in `Rat`, never realize.
    let by_frame = exact_frame(model, &frame);
    (by_value, by_frame)
}

/// The road as it was before the repair — realize the axes, then lift them back. Kept so the
/// repair is visibly a different *question* rather than a refactor of the same one.
fn both_realized(model: &mut Model, sp: SketchPlane) -> (Option<RatFrame>, Option<RatFrame>) {
    let OpOutput::DatumPlane { frame, .. } = apply(
        model,
        &Operation::DatumPlane {
            def: DatumDef::Stated(sp),
        },
    )
    .expect("every fixture here states a plane the kernel can hold") else {
        unreachable!()
    };
    (sp.exact(), realized(model, &frame).and_then(|p| p.exact()))
}

fn population() -> Vec<(&'static str, SketchPlane)> {
    let p3 = Point3::from_array;
    let v3 = Vector3::from_array;
    vec![
        ("world_xy", SketchPlane::world_xy()),
        ("world_yz", SketchPlane::world_yz()),
        // ★ The axis whose derived frame is not the convention.
        ("world_zx", SketchPlane::world_zx()),
        (
            "offset_xy",
            SketchPlane::world_xy().with_origin(p3([0.0, 0.0, 0.5])),
        ),
        (
            "far_offset_xy",
            SketchPlane::world_xy().with_origin(p3([50.0, -37.25, 0.5])),
        ),
        // ★ Rational tilt: axes lift exactly, so this takes the world road today — the only
        // family where the two routes could disagree and it would matter.
        (
            "rational_tilt_wf",
            SketchPlane::from_axes(
                p3([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                v3([0.6, 0.8, 0.0]),
                v3([-0.48, 0.36, 0.8]),
            ),
        ),
        (
            "rational_tilt_345",
            SketchPlane::from_axes(
                p3([1.0, 2.0, 3.0]),
                v3([0.6, 0.8, 0.0]),
                v3([0.0, 0.0, 1.0]),
            ),
        ),
        // Irrational tilt: `exact()` declines on both roads — the agreement that matters here
        // is that they decline *together*.
        (
            "irrational_tilt",
            SketchPlane::from_origin_normal(p3([0.0; 3]), v3([1.0, 1.0, 1.0])).expect("a plane"),
        ),
        // ★ A real population: three call sites extrude on a plane facing −ŷ.
        (
            "negative_normal",
            SketchPlane::from_origin_normal(p3([0.0, 1.3, 0.0]), v3([0.0, -1.0, 0.0]))
                .expect("a plane"),
        ),
        (
            "through_points",
            SketchPlane::through_points(
                p3([1.0, 0.0, 0.0]),
                p3([1.0, 2.0, 0.0]),
                p3([1.0, 0.0, 3.0]),
            )
            .expect("a plane"),
        ),
    ]
}

/// ★★ **They agree — once the question is asked in rationals.**
///
/// Before the repair two families diverged (`rational_tilt_wf`, `rational_tilt_345`): an axis
/// survived the frame road only when it needed no normalizing, because `reduce_direction`
/// turns a `(0.6, 0.8, 0)` axis into the primitive `(3, 4, 0)` with `uu = 25` and a
/// *realization* multiplies by a numerically computed `1/5`, landing on `0.6000000000000001`.
/// `Rat::from_decimal` lifted that, orthonormality failed, and a perfectly rational plane
/// took the frame-node road — a different arena, not an ulp.
///
/// ★ Two hypotheses died on the way, and both were mine. **Origin cancellation** (`1.6 − 1.0`)
/// is not the cause: a unit axis at `(1, 2, 3)` survives, which
/// [`a_realized_axis_survives_only_when_it_needs_no_normalizing`] pins. And the planning
/// argument "`|u_raw|² = 1`, so `inv_sqrt_exact` returns exactly one" held only for axes that
/// are *already* unit; I generalized it to every rational axis and predicted agreement.
///
/// [`RatFrame::of_plane_frame`] asks in `Rat` and never realizes — the same rule
/// `plane_frame_named` already states for `v̂`, one level up.
#[test]
fn the_two_roads_take_the_same_road() {
    for (name, sp) in population() {
        let mut m = Model::new();
        let (by_value, by_frame) = both(&mut m, sp);
        println!(
            "stat frame_road {name} exact_by_value={} exact_by_frame={}",
            by_value.is_some(),
            by_frame.is_some()
        );
        assert_eq!(
            by_value.is_some(),
            by_frame.is_some(),
            "{name}: the two roads disagree about whether the axes are rational. That is not \
                 an ulp — it moves the plane between the world road and the frame-node road, and \
                 the arena differs by whole motion nodes."
        );
    }
}

/// ★★ **The road the repair replaced, kept as a positive control.**
///
/// Realizing the axes and lifting them back still loses exactly the two families whose axes
/// need normalizing. Without this the repair would read as a refactor of one question; with
/// it, the two routes are visibly different questions and the fix is visibly load-bearing.
#[test]
fn realizing_the_axes_first_still_loses_them() {
    let mut lost = Vec::new();
    for (name, sp) in population() {
        let mut m = Model::new();
        let (by_value, by_realized) = both_realized(&mut m, sp);
        if by_value.is_some() != by_realized.is_some() {
            lost.push(name);
        }
    }
    assert_eq!(
        lost,
        vec!["rational_tilt_wf", "rational_tilt_345"],
        "the realized route is what the repair stopped using; if this list changed, the \
             measurement the repair was built on moved with it"
    );
}

/// **The mechanism, isolated** — so the repair is checked against the cause, not the symptom.
///
/// On the realized route an axis survives exactly when its reduced form is already unit;
/// moving the frame's origin changes nothing, which is what ruled out the "differencing two
/// realized points cancels the origin" explanation. Asked in `Rat`, all four are rational.
#[test]
fn a_realized_axis_survives_only_when_it_needs_no_normalizing() {
    let p3 = Point3::from_array;
    let v3 = Vector3::from_array;
    for (name, origin, u) in [
        ("unit_axis_at_origin", p3([0.0; 3]), v3([1.0, 0.0, 0.0])),
        ("unit_axis_moved", p3([1.0, 2.0, 3.0]), v3([1.0, 0.0, 0.0])),
        ("scaled_axis_at_origin", p3([0.0; 3]), v3([0.6, 0.8, 0.0])),
        (
            "scaled_axis_moved",
            p3([1.0, 2.0, 3.0]),
            v3([0.6, 0.8, 0.0]),
        ),
    ] {
        let mut m = Model::new();
        let sp = SketchPlane::from_axes(origin, u, v3([0.0, 0.0, 1.0]));
        assert!(sp.exact().is_some(), "{name}: caller axes lift");
        let (_, by_realized) = both_realized(&mut m, sp);
        let (_, by_rat) = both(&mut m, sp);
        let unit_axis = u.as_array().iter().filter(|c| **c != 0.0).count() == 1;
        println!(
            "stat frame_axis_loss {name} realized={} rational={}",
            by_realized.is_some(),
            by_rat.is_some()
        );
        assert!(
            by_rat.is_some(),
            "{name}: asked in Rat, every one of these frames is rational"
        );
        assert_eq!(
            by_realized.is_some(),
            unit_axis,
            "{name}: the realized route keeps an axis it does not have to normalize and \
                 loses one it does — whatever the origin is"
        );
    }
}

/// Where both roads *do* reach the exact road, the rational frames are identical — so the
/// loss above is the only thing standing between the two vocabularies.
#[test]
fn the_two_roads_realize_the_same_axes_where_they_agree_on_the_road() {
    let mut agreed = 0;
    for (name, sp) in population() {
        let mut m = Model::new();
        let (by_value, by_frame) = both(&mut m, sp);
        let (Some(a), Some(b)) = (by_value, by_frame) else {
            println!("stat frame_axes {name} both_declined");
            continue;
        };
        assert_eq!(
            a, b,
            "{name}: same road, different axes — the prism would be built on a different \
                 rational frame and every vertex would move"
        );
        agreed += 1;
        println!("stat frame_axes {name} identical");
    }
    assert!(
        agreed >= 5,
        "only {agreed} fixtures reached the exact road; the population stopped measuring what \
             it was chosen to measure"
    );
}

/// ★★ **The invariant the agreement above quietly rests on, stated.**
///
/// `exact()` lifts the frame's *realized* origin with `Rat::from_decimal`, so a `Named`
/// placement's origin has to survive a round trip through `f64` to come back as the same
/// rational. Every `Named` origin today came from `Rat::from_decimal` in the first place
/// (`SketchFrame::named` lifts an `f64`; a datum takes `PlaneDef`'s lifted points), so it
/// does — but that is a property of the current producers, **not of the type**. A future
/// producer that *computes* an origin would break the round trip, and the plane would move to
/// the frame-node road without anything failing.
#[test]
fn every_named_origin_survives_the_round_trip_it_is_relied_on_for() {
    for (name, sp) in population() {
        let mut m = Model::new();
        let OpOutput::DatumPlane { frame, .. } = apply(
            &mut m,
            &Operation::DatumPlane {
                def: DatumDef::Stated(sp),
            },
        )
        .expect("stated") else {
            unreachable!()
        };
        let nacre_topo::FramePlacement::Named { origin, .. } = frame.placement() else {
            panic!("{name}: a stated datum always names its placement");
        };
        for (k, r) in origin.iter().enumerate() {
            let round = nacre_exact::Rat::from_decimal(r.to_f64());
            assert_eq!(
                round,
                Some(*r),
                "{name}: origin[{k}] does not survive f64 — `exact()` lifts the realized \
                     origin, so this is what keeps the frame road on the world road"
            );
        }
    }
}
