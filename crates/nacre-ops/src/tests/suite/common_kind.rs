//! The `Common` boolean, the concurrency audit and producer-side labels.

use super::*;

/// The L-prism and an L-shaped bar lying in its notch, biting two convex corners of
/// the L's top face. The bar spans `z ∈ [0.5, 1.5]`, so its body clears the cap.
///
/// Each bite crosses **two different** edges of the cap, which is exactly why no edge
/// is pierced twice — the bar takes corners, not edges. Two chords, no closed loop.
fn l_and_notch_bar() -> (Model, Handle<Solid>, Handle<Solid>) {
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
    let __w6 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: b, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w6,
            profile: bar,
            dist: 1.0,
        },
    )
    .unwrap() else {
        unreachable!("extrude yields Extrude output")
    };
    (m, l, b)
}

/// The L-prism and a П-shaped staple straddling the L's reflex corner. The profile
/// lives in the **XZ** sketch plane and extrudes along `−y`, so the L's cap (`z = 1`)
/// is *parallel* to the extrusion axis and the staple's section there falls into two
/// pieces: one wholly inside the cap, one wrapping the corner `(1,1)`.
///
/// That parallelism is the whole point. A prism cut by a plane **perpendicular** to
/// its axis meets a face in the profile, which is connected — so every component of
/// `profile ∩ f` reaches `∂f`, and a face can never carry both an arc and a loop. Every
/// earlier attempt at such a fixture died on that.
///
/// Leg bottoms sit at `z = 0.5` and `z = 0.45`: two coplanar faces of *one* operand
/// tripped the pre-cutover door guard, exactly as `u_prism`'s staggered prongs avoid.
/// And the legs span `y ∈ [0.65, 1.3]`, not `[0.7, 1.3]`, because `(1.4, 0.7)` lies on
/// the cap's fan diagonal `y = x/2` and `segment_crosses_face` would graze it.
fn l_and_staple() -> (Model, Handle<Solid>, Handle<Solid>) {
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

/// **What the corpus actually contains in the way of concurrent vertices — and whether the
/// engine's rule for noticing them is complete.**
///
/// The four-plane work rests on one premise: a point's identity can be made a function of the
/// *set* of planes through it, because **every producer that meets the point derives the same
/// set**. Until now that was checked on a single model. This measures it against ground truth
/// (every plane asked, not the rule asking itself) over the whole fixture corpus.
///
/// Two things are asserted, and the second is the load-bearing one:
///
/// 1. **Exactly four.** No corpus point has five or more planes through it. That matters
///    because it is what makes both discovery rules complete: the trace learns `{wc} ∪ t`, and
///    with `|S| = 4` that *is* `S`. A five-plane point would leave it one short — so if this
///    ever fires, the identity rule needs the full set from somewhere else, and the message
///    says which model found it.
/// 2. **The trace's rule reproduces ground truth.** Wherever the trace would notice a
///    concurrency (a name on class `wc` that does not mention `wc`), the set it would record
///    equals the set every plane agrees on.
///
/// The count is printed rather than pinned: this is a *measurement*, and a number here would
/// only pin today's fixture list.
#[test]
fn concurrent_vertices_are_four_planes_and_the_trace_sees_all_of_them() {
    let mut boxed: Vec<(String, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
    macro_rules! fixture {
        ($name:ident) => {{
            let (m, a, b) = $name();
            boxed.push((stringify!($name).to_string(), m, a, b));
        }};
    }
    fixture!(two_boxes);
    fixture!(nested_boxes);
    fixture!(cube_and_notch);
    fixture!(stacked_cubes);
    fixture!(l_and_corner_box);
    fixture!(l_and_reflex_box);
    fixture!(l_and_inner_box);
    fixture!(l_and_popup_box);
    fixture!(l_and_notch_bar);
    fixture!(l_and_ell_stub);
    fixture!(l_and_staple);
    fixture!(l_and_dimple);
    fixture!(l_and_rod);
    fixture!(u_and_slab);
    // ...and the same shapes tilted, since a rotated operand is where concurrencies actually
    // turn up: an axis-aligned corpus would measure the easy half and call it the whole.
    macro_rules! tilted {
        ($name:ident) => {{
            let (mut m, a, b) = $name();
            let iso = rot_iso(nacre_exact::Axis::Z, 30);
            let a = transform(&mut m, a, &iso).unwrap();
            m.rebuild_adjacency();
            let b = transform(&mut m, b, &iso).unwrap();
            m.rebuild_adjacency();
            boxed.push((format!("{} (tilted)", stringify!($name)), m, a, b));
        }};
    }
    tilted!(two_boxes);
    tilted!(nested_boxes);
    tilted!(cube_and_notch);
    tilted!(stacked_cubes);
    tilted!(l_and_corner_box);
    tilted!(l_and_reflex_box);
    tilted!(l_and_inner_box);
    tilted!(l_and_popup_box);
    tilted!(l_and_notch_bar);
    tilted!(l_and_ell_stub);
    tilted!(l_and_staple);
    tilted!(l_and_dimple);
    tilted!(l_and_rod);
    tilted!(u_and_slab);

    // ★ And the models that actually have one. Without these the sweep asserts nothing: the
    // corpus above turns out to carry **no** concurrent vertex at all, so it can measure the
    // blast radius of the coming stages but not the discovery rule. These are the reported
    // model (a bar spun 45° whose bottom corner edge lands in the block's x = 0.5 plane) and
    // its variants; `rejects.rs` pins the same shape from outside.
    for (label, half_z, deg) in [
        ("four_plane (the reported model)", 0.2, 45),
        ("four_plane at 44 deg (a near miss)", 0.2, 44),
        ("four_plane, taller bar (a near miss)", 0.3, 45),
    ] {
        let (m, t, bar) = four_plane_model(half_z, deg);
        boxed.push((label.to_string(), m, t, bar));
    }

    let (mut with_any, mut total, mut trace_rule_checked) = (0usize, 0usize, 0usize);
    let mut lines_found = 0usize;
    for (name, m, a, b) in &boxed {
        let found = arrangement::concurrency_audit(m, *a, *b).unwrap();
        if !found.is_empty() {
            with_any += 1;
        }
        total += found.len();
        for c in &found {
            lines_found += c.lines.len();
            assert_eq!(
                c.planes.len(),
                4,
                "{name}: a {}-plane point at {:?} — the trace learns only `{{wc}} ∪ t`, \
                     which is four names at most, so this one would be discovered incomplete",
                c.planes.len(),
                c.planes
            );
            let triple = combinatorics::three_plane_name(c.triple).expect("a three-plane node");
            if !triple.contains(&c.wc) {
                trace_rule_checked += 1;
                let mut derived = triple.to_vec();
                derived.push(c.wc);
                derived.sort_unstable();
                assert_eq!(
                    derived, c.planes,
                    "{name}: on class {} the trace would record {derived:?} for the point \
                         named {:?}, but every plane says {:?}",
                    c.wc, c.triple, c.planes
                );
            }
        }
    }
    println!(
        "concurrency audit: {with_any}/{} fixtures carry a concurrent vertex, {total} in total",
        boxed.len()
    );
    // A sweep that quietly stops finding anything reads as agreement, so say what was actually
    // exercised — and fail if the load-bearing assertion never ran.
    assert!(
        trace_rule_checked > 0,
        "no observation reached the trace's discovery condition, so the rule was not measured"
    );
    // The same discovery must also yield the *line* aliases — three planes sharing a line show
    // up as a sub-triple of `S` that names no point. In this corpus every concurrency comes
    // from exactly that (a tool edge lying in a target plane), so each one carries one.
    assert!(
        lines_found > 0,
        "no line-sharing triple was derived, yet every concurrency here comes from one"
    );
    println!("  and {lines_found} carried a line-sharing triple");
    println!("  of which {trace_rule_checked} exercised the trace's discovery rule");
}

/// **Every producer states its side in the label frame** — the invariant family #2 restored,
/// swept over the whole two-solid corpus (prints the interesting classes with `--nocapture`).
///
/// The arrangement states its cell labels as `[*_above, *_below]` about one direction per plane
/// class: the class root's **stored surface normal** (`Seated{body_above}` and `emit_faces`'
/// `flip` are written against it). `combinatorics::side_of` answers in the root's **outward** frame
/// instead, and the two are opposite exactly when the root face is `Reversed`
/// (`orient_sign == -1`) — which no `add_cuboid` face ever is, but a face an earlier boolean
/// re-emitted flipped is. `graze_above` read `side_of` raw, so on a pocket wall it flipped the
/// wrong label bit. It no longer reads a point's side at all — [`arrangement::run_body_above`] derives
/// the occupied side from the ring's travel, and the frame term cancels there because
/// `order_along`'s direction and the label frame are defined by the same stored normal — but
/// this class stays the corpus's only crossed-frame witness, so it is what would catch a
/// producer that regresses to a raw `side_of`.
///
/// What this pins, measured before the fix:
/// - the pocket fixture has 5 `orient_sign == -1` classes carrying seated *and* graze segments
///   (4 walls + the floor); every other fixture has **no** `Reversed` root at all, which is why
///   the whole corpus passed with the frames crossed and why converting cannot regress it;
/// - the four wall classes stopped at `loop_orient_mismatch`; the floor class did **not** — its
///   four rim edges all carry a graze, so seated and graze were wrong *together*, consistently,
///   and the label survived verification while being inverted (a silent wrong, not a reject);
/// - the hand-derived sides on those classes, so a re-crossed frame fails here first.
#[test]
fn every_producer_states_its_side_in_the_label_frame() {
    let mut boxed: Vec<(&str, Model, Handle<Solid>, Handle<Solid>)> = Vec::new();
    macro_rules! fixture {
        ($name:ident) => {{
            let (m, a, b) = $name();
            boxed.push((stringify!($name), m, a, b));
        }};
    }
    fixture!(two_boxes);
    fixture!(nested_boxes);
    fixture!(cube_and_notch);
    fixture!(stacked_cubes);
    fixture!(l_and_corner_box);
    fixture!(l_and_reflex_box);
    fixture!(l_and_inner_box);
    fixture!(l_and_popup_box);
    fixture!(l_and_notch_bar);
    fixture!(l_and_ell_stub);
    fixture!(l_and_staple);
    fixture!(l_and_dimple);
    fixture!(l_and_rod);
    fixture!(u_and_slab);
    {
        // The pocket family: `pocketed_cube` is itself a boolean result, so its pocket walls
        // are `Reversed` faces. Box coordinates are `pocket_corner_cut`'s.
        let (mut m, pc) = pocketed_cube();
        let bx = m.add_cuboid(
            Point3::from_array([0.85, 0.85, 0.85]),
            Point3::from_array([1.15, 1.15, 1.15]),
        );
        boxed.push(("pocket_corner_cut", m, pc, bx));
    }

    let mut reversed_with_graze: Vec<String> = Vec::new();
    let mut reversed_seated_only: Vec<String> = Vec::new();
    for (name, m, a, b) in &boxed {
        let audits = arrangement::frame_audit(m, BoolKind::Cut, *a, *b).unwrap();
        for au in &audits {
            let interesting = au.orient_sign < 0 || au.failed_at.is_some();
            if !interesting {
                continue;
            }
            let where_ = format!(
                "{name}: wc={} pt={:?} n={:?} orient_sign={} seated={:?} graze={:?} trans={} \
                     declined={:?} failed_at={:?}",
                au.wc,
                au.root_point,
                au.root_normal,
                au.orient_sign,
                au.seated,
                au.grazes,
                au.transversals,
                au.declined,
                au.failed_at,
            );
            println!("{where_}");
            if au.orient_sign < 0 {
                if au.grazes.is_empty() {
                    if !au.seated.is_empty() {
                        reversed_seated_only.push(where_);
                    }
                } else {
                    reversed_with_graze.push(where_);
                }
            }
        }
    }
    println!("--- reversed-root classes carrying a graze (the set the fix moves) ---");
    for r in &reversed_with_graze {
        println!("  {r}");
    }
    println!("--- reversed-root classes with seated but no graze (the alternative's risk) ---");
    for r in &reversed_seated_only {
        println!("  {r}");
    }
    // The pocket fixture must keep supplying such classes, or this test has stopped exercising
    // the crossed-frame configuration and would pass vacuously.
    assert_eq!(
        reversed_with_graze.len(),
        5,
        "the pocket's 4 walls + floor are the corpus's only reversed-root classes with a graze"
    );
    assert!(
        reversed_with_graze.iter().all(|r| r.starts_with("pocket")),
        "no other fixture may have one: {reversed_with_graze:?}"
    );

    // Hand-derived sides on the pocket wall class `x = 0.7` (root = the pocket's +x wall, whose
    // outward normal points into the void, so the stored normal `+x` makes "above" the material
    // side `x > 0.7`). The wall is seated with its body above; the two side walls and the floor
    // graze it from `x < 0.7`, i.e. below. Crossed frames invert the grazes.
    //
    // The **box's top face** grazes it too, from `x > 0.7`: the pocket's opening makes that face
    // a notched region whose edge rides this plane with the material outside the pocket. That is
    // a run whose flanks *differ*, which the engine used to read as a straddling transversal —
    // this class is the corpus's only `Reversed` root, so it is also the only place the frame
    // handling of `arrangement::run_body_above` is exercised against a crossed frame: the three
    // `false` entries below are the pre-existing answers, unchanged by the new rule.
    let (m, pc, bx) = boxed
        .iter()
        .find_map(|(n, m, a, b)| (*n == "pocket_corner_cut").then_some((m, *a, *b)))
        .unwrap();
    let wall = arrangement::frame_audit(m, BoolKind::Cut, pc, bx)
        .unwrap()
        .into_iter()
        .find(|au| au.root_point == [0.7, 0.7, 1.0] && au.root_normal == [1.0, 0.0, 0.0])
        .expect("the x=0.7 pocket wall class");
    assert_eq!(wall.orient_sign, -1, "the pocket wall is a Reversed face");
    assert_eq!(
        wall.seated,
        vec![true; 4],
        "body above = material at x > 0.7"
    );
    assert_eq!(
        wall.grazes,
        vec![true, false, false, false],
        "the top face grazes from x > 0.7 (its pocket-opening edge, material outside); \
             the two side walls and the floor graze from x < 0.7"
    );
    assert_eq!(wall.failed_at, None, "the class labels consistently");
}

/// The oblique twin of the seated `Common` in `bands`: the same box, but the cylinder stands on a
/// tilted rational frame (the Pythagorean axes `u = (0.6, 0.8, 0)`, `v = (−0.48, 0.36, 0.8)`,
/// normal `(0.64, −0.48, 0.6)`), so none of the box's planes is either ⊥ or ∥ to its axis and its
/// lateral face **meets** them. Every crossing is an ellipse — the population
/// this door still names. ★ The oblique arm asks the faces first, so the fixture
/// has to be one whose faces the planes actually cross — this one's base rim straddles `z = 0` —
/// and stated on the production road: the test door's `add_cylinder` lifts irrational caps, and
/// a class with no world description is `CylinderGateUndecided`, an honest but different fact.
#[test]
fn common_rejects_an_oblique_cylinder() {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let frame = datum_frame(
        &mut m,
        crate::SketchPlane::from_axes(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.6, 0.8, 0.0]),
            Vector3::from_array([-0.48, 0.36, 0.8]),
        ),
    );
    let profile = stated(vec![circle(
        nacre_math::Point2::from_array([0.0, 0.0]),
        0.5,
    )])
    .unwrap()
    .remove(0);
    let OpOutput::Extrude { solid: cyl, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 2.0,
        },
    )
    .expect("a tilted cylinder on a rational frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert_rejects(
        || boolean_one(&mut m, BoolKind::Common, a, cyl),
        RejectReason::ObliqueCylinderCut,
    );
}

proptest! {
    #![proptest_config(proptest::test_runner::Config::with_failure_persistence(
        proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions")
    ))]
    /// Overlapping axis-aligned boxes: the intersection volume equals the
    /// independent AABB-overlap product (mixed A/B axis-aligned vertices).
    #[test]
    fn common_axis_boxes_volume_matches_aabb_overlap(
        amin in prop::array::uniform3(-5.0f64..5.0),
        aext in prop::array::uniform3(1.0f64..4.0),
        t in prop::array::uniform3(0.05f64..0.7),
        bext in prop::array::uniform3(1.0f64..4.0),
    ) {
        let amax: [f64; 3] = std::array::from_fn(|i| amin[i] + aext[i]);
        let bmin: [f64; 3] = std::array::from_fn(|i| amin[i] + t[i] * aext[i]);
        let bmax: [f64; 3] = std::array::from_fn(|i| bmin[i] + bext[i]);
        let expected: f64 = (0..3)
            .map(|i| (amax[i].min(bmax[i]) - bmin[i]).max(0.0))
            .product();
        prop_assume!(expected > 1e-3);

        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array(amin), Point3::from_array(amax));
        let b = m.add_cuboid(Point3::from_array(bmin), Point3::from_array(bmax));
        let res = boolean_one(&mut m, BoolKind::Common, a, b);
        prop_assume!(res.is_ok()); // skip rare coplanar/degenerate configs
        let r = res.unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        let vol = nacre_props::mass_props(&m, r).unwrap().volume;
        prop_assert!((vol - expected).abs() <= 1e-9 * expected.max(1.0), "{vol} vs {expected}");
    }

    /// A tilted square prism (oblique planes) intersected with a big enclosing
    /// box is the prism — exercises non-axis-aligned face normals (the in/out
    /// sign and CCW ordering) with an independent oracle (the prism's own mass).
    #[test]
    fn common_tilted_prism_with_enclosing_box_is_the_prism(
        nx in -0.5f64..0.5,
        ny in -0.5f64..0.5,
    ) {
        let plane = SketchPlane::from_origin_normal(
            Point3::origin(),
            Vector3::from_array([nx, ny, 1.0]),
        )
        .unwrap();
        let __plane = plane;
        let mut __scratch203 = Model::new();
        let __g202 = datum_frame(&mut __scratch203, __plane);
        let mut m = replay(&[
            Operation::DatumPlane { def: DatumDef::Stated(__plane) },
            Operation::Extrude {
                frame: __g202,
            profile: square(),
            dist: 1.0,
        }])
        .unwrap();
        let prism = *m.live_solids().first().unwrap();
        let vol_prism = nacre_props::mass_props(&m, prism).unwrap().volume;
        let c = m.add_cuboid(Point3::from_array([-10.0; 3]), Point3::from_array([10.0; 3]));
        let res = boolean_one(&mut m, BoolKind::Common, prism, c);
        prop_assume!(res.is_ok());
        let r = res.unwrap();
        m.rebuild_adjacency();
        prop_assert!(nacre_validate::validate(&m).is_empty());
        let vol_r = nacre_props::mass_props(&m, r).unwrap().volume;
        prop_assert!(
            (vol_r - vol_prism).abs() <= 1e-9 * vol_prism.max(1.0),
            "{vol_r} vs {vol_prism}"
        );
    }
}
