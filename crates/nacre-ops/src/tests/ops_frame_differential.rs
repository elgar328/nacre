use super::*;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
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
use nacre_topo::Surface;

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([x0, y0]),
        Point2::from_array([x1, y0]),
        Point2::from_array([x1, y1]),
        Point2::from_array([x0, y1]),
    ])
    .expect("a rectangle")
}

/// A ring with a hole, so the differential covers the inner-loop plumbing too.
fn washer() -> Profile2d {
    let ring = |a: f64, b: f64| {
        vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ]
    };
    Profile2d::with_holes(ring(0.0, 4.0), vec![ring(1.0, 3.0)]).expect("a washer")
}

/// Everything about a model that a road could change, named and indexed.
fn arena(m: &Model) -> Vec<(String, String)> {
    let mut out = vec![(
        "len".into(),
        format!(
            "{} {} {} {} {} {}",
            m.vertex_count(),
            m.edge_count(),
            m.face_count(),
            m.shell_count(),
            m.solid_count(),
            m.surface_count()
        ),
    )];
    let mut i = 0u32;
    while let Some(h) = m.vertex_handle_at(i) {
        i += 1;
        let p = m.vertex_point(h).as_array();
        out.push((
            format!("v{}", h.index()),
            format!(
                "{:x},{:x},{:x}",
                p[0].to_bits(),
                p[1].to_bits(),
                p[2].to_bits()
            ),
        ));
    }
    let mut i = 0u32;
    while let Some(h) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(h);
        out.push((
            format!("e{}", h.index()),
            format!(
                "{},{} {},{}",
                e.surfaces[0].index(),
                e.surfaces[1].index(),
                e.vertices[0].index(),
                e.vertices[1].index()
            ),
        ));
    }
    let mut i = 0u32;
    while let Some(h) = m.face_handle_at(i) {
        i += 1;
        let f = m.face(h);
        let loops: Vec<String> = std::iter::once(&f.outer)
            .chain(f.inner.iter())
            .map(|lp| {
                lp.half_edges
                    .iter()
                    .map(|he| format!("{}{}", he.edge.index(), if he.forward { "+" } else { "-" }))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect();
        out.push((
            format!("f{}", h.index()),
            format!(
                "s{} {:?} {}",
                f.surface.index(),
                f.orientation,
                loops.join(" | ")
            ),
        ));
    }
    out.push((
        "live".into(),
        m.live_solids()
            .iter()
            .map(|h| h.index().to_string())
            .collect::<Vec<_>>()
            .join(","),
    ));
    out
}

/// How many motion nodes a model holds — the "which road" signal. A world-road prism makes none.
fn nodes(m: &Model) -> usize {
    // Motion handles are only reachable through surfaces' truth; count the distinct leaves.
    let mut seen = std::collections::BTreeSet::new();
    for i in 0..m.surface_count() as u32 {
        let h = m.surface_handle_at(i).expect("in range");
        let motion = match m.surface(h) {
            Surface::Plane { motion, .. } | Surface::Cylinder { motion, .. } => *motion,
        };
        if let Some(node) = motion {
            seen.insert(node.index());
        }
    }
    seen.len()
}

fn population() -> Vec<(&'static str, SketchPlane)> {
    let p3 = Point3::from_array;
    let v3 = Vector3::from_array;
    vec![
        ("world_xy", SketchPlane::world_xy()),
        ("world_yz", SketchPlane::world_yz()),
        ("world_zx", SketchPlane::world_zx()),
        (
            "offset_xy",
            SketchPlane::world_xy().with_origin(p3([0.0, 0.0, 0.5])),
        ),
        (
            "far_offset_xy",
            SketchPlane::world_xy().with_origin(p3([50.0, -37.25, 0.5])),
        ),
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
        (
            "irrational_tilt",
            SketchPlane::from_origin_normal(p3([0.0; 3]), v3([1.0, 1.0, 1.0])).expect("a plane"),
        ),
        (
            "negative_normal",
            SketchPlane::from_origin_normal(p3([0.0, 1.3, 0.0]), v3([0.0, -1.0, 0.0]))
                .expect("a plane"),
        ),
    ]
}

/// ★★ **The gate.** Same plane, said two ways, one arena.
#[test]
fn the_two_roads_build_the_same_prism() {
    for profile_name in ["rect", "washer"] {
        for (name, sp) in population() {
            let profile = || {
                if profile_name == "rect" {
                    rect(0.1, 0.2, 2.3, 1.7)
                } else {
                    washer()
                }
            };
            let what = format!("{name}/{profile_name}");

            // The value road: today's operation.
            let mut by_value = Model::new();
            let __frame0 = datum_frame(&mut by_value, sp);
            apply(
                &mut by_value,
                &Operation::Extrude {
                    frame: __frame0,
                    profile: profile(),
                    dist: 1.5,
                },
            )
            .unwrap_or_else(|e| panic!("{what}: value road failed: {e:?}"));

            // The handle road: state the plane, then extrude on the frame it hands back.
            let mut by_frame = Model::new();
            let OpOutput::DatumPlane { frame, .. } = apply(
                &mut by_frame,
                &Operation::DatumPlane {
                    def: DatumDef::Stated(sp),
                },
            )
            .unwrap_or_else(|e| panic!("{what}: datum failed: {e:?}")) else {
                unreachable!()
            };
            extrude_on_frame(&mut by_frame, &frame, &profile(), 1.5)
                .unwrap_or_else(|e| panic!("{what}: frame road failed: {e:?}"));

            // ★ Which road, first. A frame-node road is a valid model too, so an equality check
            // alone would report the difference without naming its cause.
            assert_eq!(
                nodes(&by_frame),
                nodes(&by_value),
                "{what}: the two roads disagree about whether a motion node is needed — the \
                     arena differs by whole nodes, not by an ulp"
            );

            let (a, b) = (arena(&by_value), arena(&by_frame));
            for (x, y) in a.iter().zip(&b) {
                assert_eq!(x, y, "{what}: arenas first differ at {}", x.0);
            }
            assert_eq!(a.len(), b.len(), "{what}: different cell counts");
            println!(
                "stat frame_differential {what} identical nodes={}",
                nodes(&by_value)
            );
        }
    }
}

/// ★★ **The negative control: the ZX trap is real, and the sugar is what avoids it.**
///
/// The arbitrary-axis rule gives the ZX plane `+u = −x̂` while the convention — and
/// `SketchPlane::world_zx` — says `+u = +ẑ`. So reaching for `SketchFrame::canonical` on the
/// ZX seed puts a caller's profile a quarter turn from where they asked, and the test above,
/// which uses the datum's own `Named` frame, would never notice.
///
/// Without this, "the sugar is just `canonical`" is a simplification that passes everything.
#[test]
fn the_zx_seed_without_the_sugar_turns_the_sketch() {
    let profile = rect(0.0, 0.0, 2.0, 1.0);

    let mut stated = Model::new();
    let __w0 = SketchFrame::world(&stated, Axis::Y);
    apply(
        &mut stated,
        &Operation::Extrude {
            frame: __w0,
            profile: profile.clone(),
            dist: 1.0,
        },
    )
    .expect("the convention");

    let mut derived = Model::new();
    let zx = derived.world_plane(nacre_scalar::Axis::Y);
    extrude_on_frame(&mut derived, &SketchFrame::canonical(zx), &profile, 1.0)
        .expect("the derivation is a perfectly good frame — it is just a different one");

    let (a, b) = (arena(&stated), arena(&derived));
    assert_ne!(
        a, b,
        "the derived ZX frame must differ from the convention — if these ever agree, either \
             the seeding or the arbitrary-axis rule moved, and SketchFrame::world is dead weight"
    );
    // And the sugar is what closes it: same convention, same arena.
    let mut sugared = Model::new();
    let f = SketchFrame::world(&sugared, nacre_scalar::Axis::Y);
    extrude_on_frame(&mut sugared, &f, &profile, 1.0).expect("the sugar");
    let c = arena(&sugared);
    for (x, y) in a.iter().zip(&c) {
        assert_eq!(x, y, "the sugar must reproduce the convention: {}", x.0);
    }
    println!("stat frame_differential world_zx_without_sugar differs=true sugar_matches=true");
}
