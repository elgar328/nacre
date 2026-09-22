//! Sketched models: circles and annuli, slots, rounded rectangles, fillets, four-plane vertices.

use super::*;

// ─── A whole circle extrudes into the cylinder primitive's own solid ────────────────────────

/// **A solid's b-rep, position-canonically.** What two builders must agree on when they claim to
/// state one solid, with handle numbering left out — the two push their arenas in different
/// orders, so handles cannot be the lock; everything the handles *name* can.
#[derive(Debug, PartialEq)]
struct BrepDigest {
    vertex_bits: Vec<[u64; 3]>,
    vertex_defs: Vec<String>,
    plane_bits: Vec<[u64; 4]>,
    cylinder_defs: Vec<String>,
    /// `(surface kind, forward?, outer loop length, inner loop count)` per face.
    faces: Vec<(&'static str, bool, usize, usize)>,
    /// `(carrier kinds, curve kind, same vertex at both ends?)` per edge.
    edges: Vec<(String, &'static str, bool)>,
    volume_bits: u64,
    /// Every mesh triangle as its three positions' bits, each triangle and the list sorted.
    triangle_bits: Vec<[[u64; 3]; 3]>,
}

fn brep_digest(m: &Model, s: Handle<Solid>) -> BrepDigest {
    let bits3 = |p: Point3| p.as_array().map(f64::to_bits);
    let kind = |h: Handle<Surface>| match m.surface_cache(h) {
        nacre_geom::Surface::Plane(_) => "plane",
        nacre_geom::Surface::Cylinder(_) => "cylinder",
    };
    let sol = m.solid(s).clone();
    let (mut vertex_bits, mut vertex_defs, mut plane_bits, mut cylinder_defs) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut faces, mut edges) = (Vec::new(), Vec::new());
    let mut seen_edges: Vec<Handle<nacre_topo::Edge>> = Vec::new();
    let mut seen_vertices: Vec<Handle<Vertex>> = Vec::new();
    for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
        for fh in m.shell(sh).faces.clone() {
            let face = m.face(fh);
            match m.surface_cache(face.surface) {
                nacre_geom::Surface::Plane(pl) => {
                    plane_bits.push(pl.coefficients().map(f64::to_bits))
                }
                nacre_geom::Surface::Cylinder(_) => {
                    if let nacre_topo::Surface::Cylinder { def, .. } = m.surface(face.surface) {
                        cylinder_defs.push(format!("{def:?}"));
                    }
                }
            }
            faces.push((
                kind(face.surface),
                face.orientation == Orientation::Forward,
                face.outer.half_edges.len(),
                face.inner.len(),
            ));
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    if !seen_edges.contains(&he.edge) {
                        seen_edges.push(he.edge);
                        let e = m.edge(he.edge);
                        let mut carriers = [kind(e.surfaces[0]), kind(e.surfaces[1])];
                        carriers.sort_unstable();
                        let curve = match m.edge_curve(he.edge) {
                            nacre_geom::Curve::Line(_) => "line",
                            nacre_geom::Curve::Circle(_) => "circle",
                        };
                        edges.push((carriers.join("+"), curve, e.vertices[0] == e.vertices[1]));
                    }
                    for &vh in m.edge(he.edge).vertices.iter() {
                        if !seen_vertices.contains(&vh) {
                            seen_vertices.push(vh);
                            vertex_bits.push(bits3(m.vertex_point(vh)));
                            vertex_defs.push(match *m.vertex(vh) {
                                Vertex::ThreePlane(_) => "three-plane".to_string(),
                                Vertex::OnSeam(_) => "on-seam".to_string(),
                                Vertex::Pierce { root, .. } => format!("pierce {root:?}"),
                            });
                        }
                    }
                }
            }
        }
    }
    let mesh = nacre_tess::tessellate(m, &nacre_tess::TessConfig::default()).expect("tess");
    let mut triangle_bits: Vec<[[u64; 3]; 3]> = mesh
        .triangles
        .iter()
        .map(|(_, t)| {
            let mut tri = t.vertices.map(|h| bits3(mesh.vertices.get(h).pos));
            tri.sort_unstable();
            tri
        })
        .collect();
    vertex_bits.sort_unstable();
    vertex_bits.dedup();
    vertex_defs.sort_unstable();
    plane_bits.sort_unstable();
    plane_bits.dedup();
    cylinder_defs.sort_unstable();
    cylinder_defs.dedup();
    faces.sort_unstable();
    edges.sort_unstable();
    triangle_bits.sort_unstable();
    BrepDigest {
        vertex_bits,
        vertex_defs,
        plane_bits,
        cylinder_defs,
        faces,
        edges,
        volume_bits: nacre_props::mass_props(m, s).unwrap().volume.to_bits(),
        triangle_bits,
    }
}

/// ★★★★★ **The sugar claim, measured.** `cylinder()` is «a circle sketched, then extruded», so
/// the extrude of a whole-circle profile must state the very solid the cylinder primitive states
/// (`Model::add_cylinder_exact`, the test door): the
/// same vertex coordinates (bits), the same two seam definitions, the same
/// plane coefficients and cylinder definition, the same faces (kind, sense, loop shape), the same
/// edges (carriers, curve, the circle closing on one vertex), the same volume bits and the same
/// mesh. Handle numbering is the one thing left out — the two builders push their arenas in
/// different orders, and nothing downstream reads a handle's number.
#[test]
fn a_circle_profile_extrude_states_the_cylinder_primitive_solid() {
    let (center, radius, dist) = ([0.5, 0.25], 0.1, 2.0);
    let mut a = Model::new();
    let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
    let (sa, _) = a
        .add_cylinder_exact(
            [r(center[0]), r(center[1]), r(0.0)],
            [r(0.0), r(0.0), r(1.0)],
            [r(1.0), r(0.0), r(0.0)],
            r(radius),
            r(dist),
            None,
        )
        .expect("the primitive builds");
    let mut b = Model::new();
    let frame = SketchFrame::world(&b, Axis::Z);
    let OpOutput::Extrude { solid: sb, faces } = apply(
        &mut b,
        &Operation::Extrude {
            frame,
            profile: circle_profile(center, radius),
            dist,
        },
    )
    .expect("the circle extrudes") else {
        unreachable!()
    };
    a.rebuild_adjacency();
    b.rebuild_adjacency();
    assert!(nacre_validate::validate(&a).is_empty());
    assert!(
        nacre_validate::validate(&b).is_empty(),
        "{:?}",
        nacre_validate::validate(&b)
    );
    assert_eq!(faces.len(), 3, "base cap, top cap, lateral");
    let (da, db) = (brep_digest(&a, sa), brep_digest(&b, sb));
    assert_eq!(da.vertex_defs, vec!["on-seam", "on-seam"]);
    assert_eq!(da, db);
}

/// **A round hole and a ring**: a plate with a circular bore, and an annulus — the circle as a
/// hole ring (its lateral faces the material from outside) and as both rings of one profile.
#[test]
fn a_round_hole_and_an_annulus_extrude_exactly() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let pi = std::f64::consts::PI;
    // A 10 × 10 × 1 plate with a bore of radius 2 at (5, 5).
    let mut edges: Vec<Stated> = vec![
        line(p2(0.0, 0.0), p2(10.0, 0.0)),
        line(p2(10.0, 0.0), p2(10.0, 10.0)),
        line(p2(10.0, 10.0), p2(0.0, 10.0)),
        line(p2(0.0, 10.0), p2(0.0, 0.0)),
    ];
    edges.push(stated::circle(p2(5.0, 5.0), 2.0));
    let profiles = stated(edges).unwrap();
    assert_eq!(profiles.len(), 1);
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, faces } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: profiles[0].clone(),
            dist: 1.0,
        },
    )
    .expect("the bored plate extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((v - (100.0 - 4.0 * pi)).abs() < 1e-9, "{v}");
    assert_eq!(faces.len(), 7, "two caps, four walls, one bore");
    let bore = faces
        .iter()
        .copied()
        .find(|&f| {
            matches!(
                m.surface_cache(m.face(f).surface),
                nacre_geom::Surface::Cylinder(_)
            )
        })
        .expect("a bore face");
    assert_eq!(
        m.face(bore).orientation,
        Orientation::Forward.flipped(),
        "the material is outside the bore"
    );
    let caps_with_a_hole = faces
        .iter()
        .filter(|&&f| m.face(f).inner.len() == 1)
        .count();
    assert_eq!(caps_with_a_hole, 2);
    mesh_covers_faces("a plate with a round hole", &m, &[solid]);

    // An annulus: outer radius 5, inner 3, height 2.
    let profiles = stated(vec![
        stated::circle(p2(0.0, 0.0), 5.0),
        stated::circle(p2(0.0, 0.0), 3.0),
    ])
    .unwrap();
    assert_eq!((profiles.len(), profiles[0].holes().len()), (1, 1));
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, faces } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: profiles[0].clone(),
            dist: 2.0,
        },
    )
    .expect("the annulus extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((v - 32.0 * pi).abs() < 1e-9, "{v}");
    assert_eq!(faces.len(), 4, "two caps, the outer wall, the inner wall");
    let mut senses: Vec<bool> = faces
        .iter()
        .filter(|&&f| {
            matches!(
                m.surface_cache(m.face(f).surface),
                nacre_geom::Surface::Cylinder(_)
            )
        })
        .map(|&f| m.face(f).orientation == Orientation::Forward)
        .collect();
    senses.sort_unstable();
    assert_eq!(
        senses,
        vec![false, true],
        "the outer wall faces in, the inner wall faces out"
    );
    mesh_covers_faces("an annulus", &m, &[solid]);
}

/// **One cylinder surface under two solids is refused, not asserted.** Cylinders intern by
/// their exact statement, so two identical circle prisms — or two identical `cylinder()`
/// primitives — share one lateral surface handle, which the chart reads as one solid's. Measured
/// before the gate learned this: the fuse reached `cyl_chart`'s "one cylinder class carries rows
/// of both solids" assertion. Now it is the coaxial pair's refusal, by name.
#[test]
fn identical_cylinders_are_one_surface_and_their_boolean_is_refused_by_name() {
    let refused = |r: Result<Vec<Handle<Solid>>, BoolError>| {
        assert!(
            matches!(
                r,
                Err(BoolError::Rejected {
                    reason: RejectReason::CylinderPairContact,
                    ..
                })
            ),
            "{r:?}"
        );
    };
    // Two circle prisms.
    let mut m = Model::new();
    let extrude = |m: &mut Model| {
        let frame = SketchFrame::world(m, Axis::Z);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile: circle_profile([0.0, 0.0], 1.0),
                dist: 2.0,
            },
        )
        .unwrap() else {
            unreachable!()
        };
        solid
    };
    let (s1, s2) = (extrude(&mut m), extrude(&mut m));
    m.rebuild_adjacency();
    refused(crate::boolean(&mut m, BoolKind::Fuse, s1, s2));
    // Two primitives — the same hazard predates the sketch road.
    let mut m = Model::new();
    let prim = |m: &mut Model| {
        let r = |x: i128| nacre_exact::Rat::from_int(x);
        m.add_cylinder_exact(
            [r(0), r(0), r(0)],
            [r(0), r(0), r(1)],
            [r(1), r(0), r(0)],
            r(1),
            r(2),
            None,
        )
        .unwrap()
        .0
    };
    let (s1, s2) = (prim(&mut m), prim(&mut m));
    m.rebuild_adjacency();
    refused(crate::boolean(&mut m, BoolKind::Cut, s1, s2));
}

/// **A circle on a slanted wall's frame** — the measurement the plan asked for, locked at what it
/// found. The prism alone is a valid solid: the arc wall's f64 cache is realized through the
/// motion frame the way its vertices are, and the disk cap's winding check realizes the plane's
/// truth points the same way (a frame triangle against a world normal would assert). Padding it
/// onto the
/// wall then declines by the name a tilted `cylinder()` primitive gets today: the cylinder gate
/// cannot state a rotated cylinder against another body yet.
#[test]
fn a_circle_prism_on_a_slanted_wall_builds_and_its_pad_declines_by_name() {
    let (mut m, wall) = prism_with_a_slanted_wall();
    let sp = crate::ops::face_plane(&m, wall).expect("planar");
    let d = nacre_props::face_props(&m, wall).unwrap().centroid - sp.origin;
    let (cu, cv) = (d.dot(sp.x_axis), d.dot(sp.y_axis));
    let profile = circle_profile([cu, cv], 0.4);
    {
        let (mut m2, wall2) = prism_with_a_slanted_wall();
        let frame = crate::ops::face_sketch_frame(&m2, wall2).expect("a face frame");
        let OpOutput::Extrude { solid, .. } = apply(
            &mut m2,
            &Operation::Extrude {
                frame,
                profile: profile.clone(),
                dist: 1.0,
            },
        )
        .expect("the circle prism builds on the wall's frame") else {
            unreachable!()
        };
        m2.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m2).is_empty(),
            "{:?}",
            nacre_validate::validate(&m2)
        );
        let v = nacre_props::mass_props(&m2, solid).unwrap().volume;
        assert!((v - 0.16 * std::f64::consts::PI).abs() < 1e-9, "{v}");
    }
    let r = apply(
        &mut m,
        &Operation::PadOnFace {
            face: wall,
            profile,
            dist: 1.0,
        },
    );
    assert!(
        matches!(
            r,
            Err(OpError::Boolean(BoolError::Rejected {
                reason: RejectReason::CylinderGateUndecided,
                ..
            }))
        ),
        "{r:?}"
    );
}

// ─── Arcs between vertices — slot, rounded rectangle, D ─────────────────────────────────────

fn edges_profile(edges: Vec<Stated>) -> Profile2d {
    stated(edges).unwrap().remove(0)
}

fn slot_profile(cx0: f64, cx1: f64, cy: f64, r: f64) -> Profile2d {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    edges_profile(vec![
        line(p2(cx0, cy - r), p2(cx1, cy - r)),
        arc_turns(p2(cx1, cy), p2(cx1, cy - r), 2),
        line(p2(cx1, cy + r), p2(cx0, cy + r)),
        arc_turns(p2(cx0, cy), p2(cx0, cy + r), 2),
    ])
}

fn extrude_world_z(
    m: &mut Model,
    profile: Profile2d,
    dist: f64,
) -> (Handle<Solid>, Vec<Handle<Face>>) {
    let frame = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, faces } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("the profile extrudes") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    (solid, faces)
}

/// The definitions of a solid's vertices, as the digest spells them, sorted.
fn vertex_def_names(m: &Model, s: Handle<Solid>) -> Vec<String> {
    brep_digest(m, s).vertex_defs
}

/// ★★★★★ **A slot stands.** Two straight walls tangent to two half cylinders: every corner is a
/// `Pierce` whose root is the **double** one — the wall's plane touches the cylinder along the
/// ruling through that corner, which is what "tangent" says. Volume `(2rL + πr²)·h`, six faces,
/// valid, and the mesh covers it.
#[test]
fn a_slot_extrudes_with_tangent_pierce_corners() {
    let pi = std::f64::consts::PI;
    let mut m = Model::new();
    let (solid, faces) = extrude_world_z(&mut m, slot_profile(0.0, 30.0, 0.0, 5.0), 2.0);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((v - (300.0 + 25.0 * pi) * 2.0).abs() < 1e-9, "{v}");
    assert_eq!(
        faces.len(),
        6,
        "two caps, two straight walls, two half cylinders"
    );
    assert_eq!(vertex_def_names(&m, solid), vec!["pierce Double"; 8]);
    let cylinders = faces
        .iter()
        .filter(|&&f| {
            matches!(
                m.surface_cache(m.face(f).surface),
                nacre_geom::Surface::Cylinder(_)
            )
        })
        .count();
    assert_eq!(cylinders, 2);
    mesh_covers_faces("a slot", &m, &[solid]);
}

/// **A rounded rectangle**: four straight walls, four quarter cylinders, sixteen tangent corners.
/// Volume `(wh − (4 − π)r²)·h`.
#[test]
fn a_rounded_rectangle_extrudes() {
    let pi = std::f64::consts::PI;
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let l = |a: [f64; 2], b: [f64; 2]| line(p2(a[0], a[1]), p2(b[0], b[1]));
    let q = |c: [f64; 2], s: [f64; 2]| arc_turns(p2(c[0], c[1]), p2(s[0], s[1]), 1);
    let profile = edges_profile(vec![
        l([5.0, 0.0], [35.0, 0.0]),
        q([35.0, 5.0], [35.0, 0.0]),
        l([40.0, 5.0], [40.0, 15.0]),
        q([35.0, 15.0], [40.0, 15.0]),
        l([35.0, 20.0], [5.0, 20.0]),
        q([5.0, 15.0], [5.0, 20.0]),
        l([0.0, 15.0], [0.0, 5.0]),
        q([5.0, 5.0], [0.0, 5.0]),
    ]);
    let mut m = Model::new();
    let (solid, faces) = extrude_world_z(&mut m, profile, 1.0);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((v - (800.0 - (4.0 - pi) * 25.0)).abs() < 1e-9, "{v}");
    assert_eq!(faces.len(), 10);
    assert_eq!(vertex_def_names(&m, solid), vec!["pierce Double"; 16]);
    mesh_covers_faces("a rounded rectangle", &m, &[solid]);
}

/// **A D**: a half disk closed by its diameter. The chord's wall plane runs *through* the axis, so
/// it crosses the cylinder in two rulings — the corners are the pair's `Lo` and `Hi`, not a
/// tangency. Volume `πr²h/2`, four faces.
#[test]
fn a_half_disk_has_the_pair_roots_at_its_corners() {
    let pi = std::f64::consts::PI;
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let profile = edges_profile(vec![
        line(p2(0.0, 5.0), p2(0.0, -5.0)),
        arc_turns(p2(0.0, 0.0), p2(0.0, -5.0), 2), // through (5, 0)
    ]);
    let mut m = Model::new();
    let (solid, faces) = extrude_world_z(&mut m, profile, 2.0);
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).unwrap().volume;
    assert!((v - 25.0 * pi).abs() < 1e-9, "{v}");
    assert_eq!(
        faces.len(),
        4,
        "two caps, the chord wall, the half cylinder"
    );
    let defs = vertex_def_names(&m, solid);
    assert_eq!(defs.len(), 4);
    assert!(
        defs.contains(&"pierce Lo".to_string()) && defs.contains(&"pierce Hi".to_string()),
        "{defs:?}"
    );
    mesh_covers_faces("a half disk", &m, &[solid]);
}

/// **A slot padded onto and pocketed into a plate** — the tangent ruling is a carrier with
/// its own name (`side = 0`).
/// The slot sits on the plate's corner: `cx ∈ [−10, 10]`, `r = 2` about the
/// face's origin, so the part over the plate is a `10 × 2` strip plus a quarter disk. The pad
/// adds the whole prism (it touches the plate along that part and hangs past its edges), the
/// pocket removes the part inside: `8000 + 3·(80 + 4π)` and `8000 − 2·(20 + π)`.
#[test]
fn a_slot_pads_and_pockets_at_the_tangent_ruling() {
    let plate = || {
        let mut m = Model::new();
        let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let square = Profile2d::polygon(vec![
            p2(0.0, 0.0),
            p2(40.0, 0.0),
            p2(40.0, 40.0),
            p2(0.0, 40.0),
        ])
        .unwrap();
        let (solid, faces) = extrude_world_z(&mut m, square, 5.0);
        (m, solid, faces[1])
    };
    let slot = || slot_profile(-10.0, 10.0, 0.0, 2.0);
    let pi = std::f64::consts::PI;
    let built = |m: &mut Model, r: Result<OpOutput, OpError>, want: f64| {
        let solid = match r.expect("the slot builds at its tangent rulings") {
            OpOutput::PadOnFace { solid, .. } | OpOutput::PocketOnFace { solid, .. } => solid,
            other => panic!("{other:?}"),
        };
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(m).is_empty(),
            "{:?}",
            nacre_validate::validate(m)
        );
        let v = nacre_props::mass_props(m, solid).expect("props").volume;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    };
    let (mut m, _, top) = plate();
    let r = apply(
        &mut m,
        &Operation::PadOnFace {
            face: top,
            profile: slot(),
            dist: 3.0,
        },
    );
    built(&mut m, r, 8000.0 + 3.0 * (80.0 + 4.0 * pi));
    let (mut m, _, top) = plate();
    let r = apply(
        &mut m,
        &Operation::PocketOnFace {
            face: top,
            profile: slot(),
            dist: 2.0,
        },
    );
    built(&mut m, r, 8000.0 - 2.0 * (20.0 + pi));
}

/// **A half-disk boss pads onto a plate.** The D's chord wall crosses its own cylinder through
/// the axis (`Lo`/`Hi` corners) and the cylinder carries only half its circle — the
/// split (no phantom piece for the uncovered half) and the
/// chart (an uncovered piece is the face's end) build it. Half the D hangs past the plate's
/// `y = 0` edge; the pad adds the whole half disk: `8000 + 3·(π·25/2)`.
#[test]
fn a_half_disk_boss_pads_onto_a_plate() {
    let mut m = Model::new();
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let square = Profile2d::polygon(vec![
        p2(0.0, 0.0),
        p2(40.0, 0.0),
        p2(40.0, 40.0),
        p2(0.0, 40.0),
    ])
    .unwrap();
    let (_, faces) = extrude_world_z(&mut m, square, 5.0);
    let d_shape = edges_profile(vec![
        line(p2(0.0, 5.0), p2(0.0, -5.0)),
        arc_turns(p2(0.0, 0.0), p2(0.0, -5.0), 2),
    ]);
    let r = apply(
        &mut m,
        &Operation::PadOnFace {
            face: faces[1],
            profile: d_shape,
            dist: 3.0,
        },
    );
    let OpOutput::PadOnFace { solid, .. } = r.expect("the half disk pads") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, solid).expect("props").volume;
    let want = 8000.0 + 3.0 * (std::f64::consts::PI * 25.0 / 2.0);
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// **A sketched bored plate meets a box** — whole circles are the arc vocabulary the boolean
/// already reads (a bore's rim is a circle), so the plate with a round hole drawn as a sketch
/// fuses and cuts like a drilled one. Plate `10 × 10 × 1` minus `π·2²`, a `4 × 2 × 3` box
/// straddling its left edge (overlap `2 × 2 × 1`).
#[test]
fn a_sketched_bored_plate_meets_a_box() {
    let pi = std::f64::consts::PI;
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let bored = || {
        edges_profile(vec![
            line(p2(0.0, 0.0), p2(10.0, 0.0)),
            line(p2(10.0, 0.0), p2(10.0, 10.0)),
            line(p2(10.0, 10.0), p2(0.0, 10.0)),
            line(p2(0.0, 10.0), p2(0.0, 0.0)),
            stated::circle(p2(5.0, 5.0), 2.0),
        ])
    };
    for (kind, want) in [
        (BoolKind::Fuse, 100.0 - 4.0 * pi + 24.0 - 4.0),
        (BoolKind::Cut, 100.0 - 4.0 * pi - 4.0),
        (BoolKind::Common, 4.0),
    ] {
        let mut m = Model::new();
        let (plate, _) = extrude_world_z(&mut m, bored(), 1.0);
        let box_ = m.add_cuboid(
            Point3::from_array([-2.0, 4.0, -1.0]),
            Point3::from_array([2.0, 6.0, 2.0]),
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, kind, plate, box_).unwrap_or_else(|e| panic!("{kind:?}: {e:?}"));
        m.rebuild_adjacency();
        assert_eq!(out.len(), 1, "{kind:?}");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{kind:?}: {:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - want).abs() < 1e-9, "{kind:?}: {v} vs {want}");
        mesh_covers_faces("a sketched bored plate and a box", &m, &out);
    }
}

/// ★ **A filleted plate builds beside a far box, and the decline probe reads
/// nothing.** The wall faces that ride a tangent ruling get their outer loop: the tangent
/// ruling is a carrier with `side = 0`, the seated wall states its own piece, the split keeps no
/// phantom arc, the angular order and the winding walk read the smooth corner as a half turn, and
/// the chart reads an uncovered rim piece as the face's end. Two bodies, the plate's volume
/// `7·8 − 2·(1.5² − π·1.5²/4)` plus the box's 1.
#[test]
fn a_filleted_plate_builds_and_the_decline_probe_reads_nothing() {
    use crate::arrangement::decline_probe::ROWS;
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let edges = vec![
        line(p2(-2.0, -4.0), p2(2.0, -4.0)),
        arc_turns(p2(2.0, -2.5), p2(2.0, -4.0), 1),
        line(p2(3.5, -2.5), p2(3.5, 4.0)),
        line(p2(3.5, 4.0), p2(-3.5, 4.0)),
        line(p2(-3.5, 4.0), p2(-3.5, -2.5)),
        arc_turns(p2(-2.0, -2.5), p2(-3.5, -2.5), 1),
    ];
    let profile = stated(edges).unwrap().remove(0);
    let mut m = Model::new();
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: plate, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("the filleted plate extrudes") else {
        unreachable!()
    };
    let far = m.add_cuboid(
        Point3::from_array([20.0, 20.0, 20.0]),
        Point3::from_array([21.0, 21.0, 21.0]),
    );
    m.rebuild_adjacency();
    // Owned, so the declines this fixture would produce are this test's and not a neighbour's.
    let out = crate::ledger::owned(|| boolean(&mut m, BoolKind::Fuse, plate, far));
    let out = out.expect("the filleted plate and the far box fuse");
    assert_eq!(out.len(), 2, "two bodies, apart");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    let want = 56.0 - 2.0 * (2.25 - std::f64::consts::PI * 2.25 / 4.0) + 1.0;
    assert!((v - want).abs() < 1e-9, "plate + box: {v} vs {want}");
    assert!(
        ROWS.mine().is_empty(),
        "no decline was produced: {:?}",
        ROWS.mine()
    );
}

/// ★ **The user's plate, slot plate and gusset fold into one body.** The filleted,
/// bored plate (p1), the slot-windowed vertical plate standing on it (p2) and a triangular gusset
/// (p3) — every wall the three walls found: the fillets' tangent rulings, the slot's tangent
/// rulings, the plate's caps crossing the gusset's side plane inside the fillet's radius, the
/// nesting cell between a bore's and a fillet's rulings whose every corner is irrational. The
/// slot plate and the gusset stand **on** the plate (their bottoms lie in its top face), so the
/// volume is the sum of the parts: `48 + (34 − π/4) + 7.5` — the plate less its fillets' corners
/// and bores, the slot plate less its window, the gusset.
///
/// ★ The gusset stands at `x = 1.6` here; at the user's `1.5` its side plane runs **through the
/// fillet's axis** and so contains the fillet's tangent ruling — one line shared by two planes
/// and a cylinder. The second arm below is the plate and the
/// gusset at `1.5` fusing: `48 + 7.5`. The four-part fold at the script's own dimensions is
/// `the_users_fold_builds_at_the_scripts_own_dimensions`.
#[test]
fn the_users_plate_slot_and_gusset_fold() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let plate = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                line(p2(-1.5, -4.0), p2(1.5, -4.0)),
                arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                line(p2(3.5, -2.0), p2(3.5, 4.0)),
                line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                line(p2(-3.5, 4.0), p2(-3.5, -2.0)),
                arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1),
                stated::circle(p2(-1.5, -2.0), 1.0),
                stated::circle(p2(1.5, -2.0), 1.0),
            ],
            1.0,
        )
    };
    let slot_plate = |m: &mut Model| {
        let s = prism(
            m,
            Axis::Y,
            vec![
                line(p2(1.0, -3.5), p2(6.0, -3.5)),
                line(p2(6.0, -3.5), p2(6.0, 3.5)),
                line(p2(6.0, 3.5), p2(1.0, 3.5)),
                line(p2(1.0, 3.5), p2(1.0, -3.5)),
                line(p2(3.0, -0.5), p2(4.0, -0.5)),
                arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2),
                line(p2(4.0, 0.5), p2(3.0, 0.5)),
                arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2),
            ],
            1.0,
        );
        shift(m, s, [0.0, 3.0, 0.0])
    };
    let gusset = |m: &mut Model, x: f64| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.0, 1.0), p2(3.0, 1.0)),
                line(p2(3.0, 1.0), p2(3.0, 6.0)),
                line(p2(3.0, 6.0), p2(0.0, 1.0)),
            ],
            1.0,
        );
        shift(m, g, [x, 0.0, 0.0])
    };
    let mut m = Model::new();
    let a = plate(&mut m);
    let b = slot_plate(&mut m);
    let c = gusset(&mut m, 1.6);
    m.rebuild_adjacency();
    let ab = boolean(&mut m, BoolKind::Fuse, a, b).expect("plate + slot plate");
    assert_eq!(ab.len(), 1);
    m.rebuild_adjacency();
    let abc = boolean(&mut m, BoolKind::Fuse, ab[0], c).expect("+ gusset");
    assert_eq!(abc.len(), 1);
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, abc[0]).expect("props").volume;
    let pi = std::f64::consts::PI;
    let want = 48.0 + (34.0 - pi / 4.0) + 7.5;
    assert!((v - want).abs() < 1e-9, "the fold: {v} vs {want}");

    // The gusset on the fillet's axis plane: the shared line has one name, and the
    // gusset stands on the plate — the volume is the sum.
    let mut m = Model::new();
    let a = plate(&mut m);
    let c = gusset(&mut m, 1.5);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, c).expect("plate + gusset on the axis plane");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    assert!(
        (v - 55.5).abs() < 1e-9,
        "plate + gusset at 1.5: {v} vs 55.5"
    );
}

/// ★ **The gate reads the arc, not the circle.** A lateral face's footprint has an
/// angular extent, and the two clearance questions read it rather than a whole circle: the
/// oblique plane's reach (a quarter-cylinder sector beside a gusset whose slanted plane runs
/// within the radius but past the arc) and the parallel pair's cross-section (two half-cylinder
/// prisms whose infinite surfaces cross while their arcs face away). Each has its negative
/// control: the same plane through the arc, the same pair with the arcs facing each other.
#[test]
fn the_gate_reads_the_arc_not_the_circle() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let volume = |m: &Model, s: Handle<Solid>| nacre_props::mass_props(m, s).expect("props").volume;
    let clean = |m: &Model| {
        assert!(
            nacre_validate::validate(m).is_empty(),
            "{:?}",
            nacre_validate::validate(m)
        );
    };
    let rejects = |out: Result<Vec<Handle<Solid>>, BoolError>, want: RejectReason| {
        assert!(
            matches!(&out, Err(BoolError::Rejected { reason, .. }) if *reason == want),
            "{out:?} vs {want:?}"
        );
    };

    // (a) A quarter-cylinder sector — the third quadrant of a radius-1.5 disk about the z axis,
    // its two flat walls through the axis — and a gusset standing beside it at x ∈ [2, 3] whose
    // slanted face `2y − z = 1.5` sits at `y ∈ [0.75, 1.25]` over the sector's height: within the
    // radius of the axis, so the infinite cylinder is cut, but the sector's arc never reaches
    // `y > 0`. Two bodies that do not touch.
    let sector = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                line(p2(0.0, 0.0), p2(-1.5, 0.0)),
                arc_turns(p2(0.0, 0.0), p2(-1.5, 0.0), 1),
                line(p2(0.0, -1.5), p2(0.0, 0.0)),
            ],
            1.0,
        )
    };
    let gusset = |m: &mut Model, mirror_y: f64| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.75 * mirror_y, 0.0), p2(3.0 * mirror_y, 0.0)),
                line(p2(3.0 * mirror_y, 0.0), p2(3.0 * mirror_y, 4.5)),
                line(p2(3.0 * mirror_y, 4.5), p2(0.75 * mirror_y, 0.0)),
            ],
            1.0,
        );
        shift(m, g, [2.0, 0.0, 0.0])
    };
    let mut m = Model::new();
    let a = sector(&mut m);
    let b = gusset(&mut m, 1.0);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the plane runs past the arc");
    assert_eq!(out.len(), 2, "two bodies, apart");
    m.rebuild_adjacency();
    clean(&m);
    let v: f64 = out.iter().map(|&s| volume(&m, s)).sum();
    let want = 0.5625 * std::f64::consts::PI + 5.0625;
    assert!((v - want).abs() < 1e-9, "sector + gusset: {v} vs {want}");
    // Negative control: the gusset mirrored in y, its slanted face `−2y − z = 1.5` at
    // `y ∈ [−1.25, −0.75]` — through the sector's arc. An ellipse the kernel does not build.
    let mut m = Model::new();
    let a = sector(&mut m);
    let b = gusset(&mut m, -1.0);
    m.rebuild_adjacency();
    rejects(
        boolean(&mut m, BoolKind::Fuse, a, b),
        RejectReason::ObliqueCylinderCut,
    );

    // (b) Two half-cylinder prisms of radius 1.5 on parallel axes 2.5 apart — closer than the
    // radii's sum, so the infinite surfaces cross — with the flat sides facing each other and the
    // arcs facing away. Their rims' arcs share no point of the cross-section, and the fuse is two
    // bodies.
    let half = |m: &mut Model, cx: f64, bulge: f64| {
        // The chord is `x = cx`, `y ∈ [0, 3]`; the arc bulges toward `bulge · x`.
        let (top, bottom) = (p2(cx, 3.0), p2(cx, 0.0));
        let (start, chord) = if bulge < 0.0 {
            (top, line(bottom, top))
        } else {
            (bottom, line(top, bottom))
        };
        prism(
            m,
            Axis::Z,
            vec![arc_turns(p2(cx, 1.5), start, 2), chord],
            1.0,
        )
    };
    let mut m = Model::new();
    let a = half(&mut m, 0.0, -1.0);
    let b = half(&mut m, 2.5, 1.0);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the arcs face away");
    assert_eq!(out.len(), 2, "two bodies, apart");
    m.rebuild_adjacency();
    clean(&m);
    let v: f64 = out.iter().map(|&s| volume(&m, s)).sum();
    let want = 2.25 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "two half disks: {v} vs {want}");
    // Negative control: the second prism moved to the other side, its arc facing the first's —
    // the arcs cross, the bodies overlap, and the curve between two cylinders is not built yet.
    let mut m = Model::new();
    let a = half(&mut m, 0.0, -1.0);
    let b = half(&mut m, -2.5, 1.0);
    m.rebuild_adjacency();
    rejects(
        boolean(&mut m, BoolKind::Fuse, a, b),
        RejectReason::CylinderPairContact,
    );
}

/// ★ **The gate reads faces at every one of its four sites.** Three fixtures a gate reading
/// *surfaces* would refuse, and one the gate refuses for a fact about faces — then three shapes
/// such a refusal hides from the arrangement: a solid with two
/// coaxial cylinders sectioned by a wall (the annulus), a cap whose merged region is bounded by a
/// circle (the stacked pin's cut), and a component with no vertex anywhere (the tube).
#[test]
fn the_gate_reads_faces_at_every_site() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let volume = |m: &Model, s: Handle<Solid>| nacre_props::mass_props(m, s).expect("props").volume;
    let circle = |m: &mut Model, c: [f64; 2], r: f64, h: f64| {
        prism(m, Axis::Z, vec![stated::circle(p2(c[0], c[1]), r)], h)
    };
    let clean = |m: &Model| {
        assert!(
            nacre_validate::validate(m).is_empty(),
            "{:?}",
            nacre_validate::validate(m)
        );
    };

    // (a) Stacked coaxial cylinders of different radii: parallel axes, the smaller strictly
    // inside the larger as surfaces — the pair rule clears them outright (`cylinders_nested`),
    // and the fuse joins them across the shared cap plane.
    let mut m = Model::new();
    let boss = circle(&mut m, [0.0, 0.0], 2.0, 2.0);
    let pin = circle(&mut m, [0.0, 0.0], 1.0, 2.0);
    let pin = shift(&mut m, pin, [0.0, 0.0, 2.0]);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, boss, pin).expect("stacked cylinders fuse");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    clean(&m);
    let v = volume(&m, out[0]);
    let want = 10.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "boss + pin: {v} vs {want}");

    // (b) A plate with two holes and a gusset whose slanted plane runs past them: the oblique
    // arm asks the faces — each hole's reach along the gusset normal misses the plane's station —
    // and the fuse is the plate plus the gusset standing on it.
    let plate = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                line(p2(-3.5, -4.0), p2(3.5, -4.0)),
                line(p2(3.5, -4.0), p2(3.5, 4.0)),
                line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                line(p2(-3.5, 4.0), p2(-3.5, -4.0)),
                stated::circle(p2(-1.5, -2.0), 1.0),
                stated::circle(p2(1.5, -2.0), 1.0),
            ],
            1.0,
        )
    };
    let gusset = |m: &mut Model, at: [f64; 3]| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.0, 1.0), p2(3.0, 1.0)),
                line(p2(3.0, 1.0), p2(3.0, 6.0)),
                line(p2(3.0, 6.0), p2(0.0, 1.0)),
            ],
            1.0,
        );
        shift(m, g, at)
    };
    let mut m = Model::new();
    let a = plate(&mut m);
    let b = gusset(&mut m, [1.5, 0.0, 0.0]);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the gusset stands beside the holes");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    clean(&m);
    let v = volume(&m, out[0]);
    let want = 56.0 - 2.0 * std::f64::consts::PI + 7.5;
    assert!(
        (v - want).abs() < 1e-9,
        "plate − holes + gusset: {v} vs {want}"
    );

    // (c) Negative control: the same gusset moved so its slanted plane runs through a hole's
    // lateral (at z ∈ [0, 1] the plane sits at y ∈ [−2.6, −2], inside the hole's y ∈ [−3, −1]).
    // Where a plane meets a lateral obliquely the curve is an ellipse the kernel does not build,
    // and the name still says so.
    let mut m = Model::new();
    let a = plate(&mut m);
    let b = gusset(&mut m, [1.5, -2.0, 0.0]);
    m.rebuild_adjacency();
    assert!(matches!(
        boolean(&mut m, BoolKind::Fuse, a, b),
        Err(BoolError::Rejected {
            reason: RejectReason::ObliqueCylinderCut,
            ..
        })
    ));

    // (d) An annulus (one solid, two coaxial cylinders) and a box whose walls section both: the
    // annular cap's chord on a wall class is **two pieces**, carved by the bore — the disk outer
    // joined the parity sweep for this. No closed form for the volume (a circle clipped by a
    // rectangle), so the oracle is the algebra of the three results against the operands.
    let annulus = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                stated::circle(p2(2.0, 2.0), 3.0),
                stated::circle(p2(2.0, 2.0), 1.5),
            ],
            3.0,
        )
    };
    let box_b = |m: &mut Model| {
        m.add_cuboid(
            Point3::from_array([1.0, 0.0, 1.0]),
            Point3::from_array([6.0, 4.0, 5.0]),
        )
    };
    let mut vols = Vec::new();
    for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
        let mut m = Model::new();
        let a = annulus(&mut m);
        let b = box_b(&mut m);
        m.rebuild_adjacency();
        let (va, vb) = (volume(&m, a), volume(&m, b));
        let out = boolean(&mut m, kind, a, b).expect("the annulus and the box combine");
        assert_eq!(out.len(), 1, "{kind:?}");
        m.rebuild_adjacency();
        clean(&m);
        vols.push((va, vb, volume(&m, out[0])));
    }
    let (va, vb, fuse) = vols[0];
    let cut = vols[1].2;
    let common = vols[2].2;
    assert!(
        (va - 20.25 * std::f64::consts::PI).abs() < 1e-9 && (vb - 80.0).abs() < 1e-9,
        "operands: {va} {vb}"
    );
    assert!(
        (fuse + common - (va + vb)).abs() < 1e-9 && (cut + common - va).abs() < 1e-9,
        "fuse {fuse} cut {cut} common {common} vs {va} + {vb}"
    );
    assert!(common > 0.0 && cut < va, "the box does bite the annulus");

    // (e) The stacked pin **cut** from the boss: the pin only touches the boss's cap, so the
    // result is the boss — and on the shared plane the pin's disk fills the boss cap's hole, a
    // merged region whose outer bound is the boss's circle. Unmerged, the pin's seam vertex kept
    // naming a cylinder the result has no face on (`VertexNamesAbsentSurface`, measured).
    let mut m = Model::new();
    let boss = circle(&mut m, [0.0, 0.0], 2.0, 2.0);
    let pin = circle(&mut m, [0.0, 0.0], 1.0, 2.0);
    let pin = shift(&mut m, pin, [0.0, 0.0, 2.0]);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Cut, boss, pin).expect("the pin cuts nothing off the boss");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    clean(&m);
    let v = volume(&m, out[0]);
    let want = 8.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "the boss alone: {v} vs {want}");

    // (f) A bushing: a tube and the coaxial pin standing in its bore with a gap — nested as
    // surfaces, clear as faces. The fuse is two bodies, and telling which is inside which needs a
    // witness on a component that has **no vertex at all**: the tube's annular cap names a point
    // between its rims.
    let tube = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                stated::circle(p2(2.0, 2.0), 3.0),
                stated::circle(p2(2.0, 2.0), 1.5),
            ],
            3.0,
        )
    };
    let mut m = Model::new();
    let a = tube(&mut m);
    let pin = circle(&mut m, [2.0, 2.0], 1.0, 5.0);
    let pin = shift(&mut m, pin, [0.0, 0.0, -1.0]);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, a, pin).expect("the tube and the pin fuse");
    assert_eq!(out.len(), 2, "two bodies, neither inside the other");
    m.rebuild_adjacency();
    clean(&m);
    let v: f64 = out.iter().map(|&s| volume(&m, s)).sum();
    let want = 25.25 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-9, "tube + pin: {v} vs {want}");
    let mut m = Model::new();
    let a = tube(&mut m);
    let pin = circle(&mut m, [2.0, 2.0], 1.0, 5.0);
    let pin = shift(&mut m, pin, [0.0, 0.0, -1.0]);
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Common, a, pin).expect("a common of nothing is empty");
    assert!(out.is_empty(), "{out:?}");
}

/// ★★★★★ **The common perpendicular is the only direction that sees this pair apart.**
///
/// A 20 × 20 plate 8 thick, every corner filleted `5`, and a `d 7` drill laid along `x` through
/// it. Each fillet's axis stands `5` from the drill's — inside the radius sum `5 + 7/2`, so the
/// surface rule cannot clear them — and neither axis separates the faces: along the fillet's own
/// axis both cover the plate's thickness, along the drill's the fillet sits inside its span. What
/// does separate them is the rulings' cross product, `ŷ`: a fillet's quarter arc reaches
/// `y ∈ [−10, −5]` (or its mirror) and the drill only `[−7/2, 7/2]`.
///
/// ★ **The kernel's own corpus reaches this predicate nowhere else** — measured:
/// the whole suite made 67 face-pair questions and not one needed a third direction. So
/// this fixture is the kernel's only end-to-end hold on the rule, and the volume is the oracle:
/// the plate `(400 − 4(25 − 25π/4)) · 8` less the drill's `π(7/2)² · 20`.
#[test]
fn only_the_common_perpendicular_separates_a_fillet_from_a_crosswise_drill() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    // The plate sits at `x ∈ [2, 22]` so the drill, which starts on the YZ plane, runs clear
    // through it: fillet axes at `(7, ±5)` and `(17, ±5)`.
    let (lo, hi, r) = (2.0, 22.0, 5.0);
    let plate_edges = vec![
        line(p2(lo + r, -10.0), p2(hi - r, -10.0)),
        arc_turns(p2(hi - r, -10.0 + r), p2(hi - r, -10.0), 1),
        line(p2(hi, -10.0 + r), p2(hi, 10.0 - r)),
        arc_turns(p2(hi - r, 10.0 - r), p2(hi, 10.0 - r), 1),
        line(p2(hi - r, 10.0), p2(lo + r, 10.0)),
        arc_turns(p2(lo + r, 10.0 - r), p2(lo + r, 10.0), 1),
        line(p2(lo, 10.0 - r), p2(lo, -10.0 + r)),
        arc_turns(p2(lo + r, -10.0 + r), p2(lo, -10.0 + r), 1),
    ];
    let mut m = Model::new();
    let plate = prism(&mut m, Axis::Z, plate_edges, 8.0);
    // On the YZ plane `+u = ŷ`, `+v = ẑ`: the drill's axis is `(y, z) = (0, 4)`, the plate's
    // mid-thickness, and `7/2 < 4` keeps it clear of the caps — a wider drill would meet them
    // in a ruling at an irrational `y`, which is another wall entirely.
    let drill = prism(
        &mut m,
        Axis::X,
        vec![stated::circle(p2(0.0, 4.0), 3.5)],
        30.0,
    );
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Cut, plate, drill).expect("the drill crosses the plate");
    assert_eq!(out.len(), 1, "one plate with a hole through it");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
    let pi = std::f64::consts::PI;
    let want = (400.0 - 4.0 * (25.0 - 25.0 * pi / 4.0)) * 8.0 - pi * 3.5 * 3.5 * 20.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}

/// ★ **A four-plane operand vertex has one name.** A gusset whose apex lands on the
/// wall's top edge leaves the fused operand with a vertex where four faces meet. Named once per
/// face loop — its own plane and the two edges' walls — that vertex would reach the tracer under
/// **four names**, one of them (the top face's: front, top, slant) three planes sharing a line
/// that name no point, which the judge reads as lying on every class, folding everything onto
/// the wall's far corner. The vertex names itself from the classes its topology knows
/// (`canonical_triple`), so every loop hands the tracer the same name, no name is dependent, and
/// the alias fold of a full trace lands
/// on the vertex. Locked through the audit (`operand_vertex_audit`), the instrument that counts
/// the names.
#[test]
fn a_four_plane_operand_vertex_has_one_name() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let mut m = Model::new();
    let w = prism(
        &mut m,
        Axis::Y,
        vec![
            line(p2(1.0, -3.5), p2(6.0, -3.5)),
            line(p2(6.0, -3.5), p2(6.0, 3.5)),
            line(p2(6.0, 3.5), p2(1.0, 3.5)),
            line(p2(1.0, 3.5), p2(1.0, -3.5)),
        ],
        1.0,
    );
    let w = shift(&mut m, w, [0.0, 3.0, 0.0]);
    let gusset = |m: &mut Model, x: f64| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.0, 1.0), p2(3.0, 1.0)),
                line(p2(3.0, 1.0), p2(3.0, 6.0)),
                line(p2(3.0, 6.0), p2(0.0, 1.0)),
            ],
            1.0,
        );
        shift(m, g, [x, 0.0, 0.0])
    };
    let g1 = gusset(&mut m, 1.4);
    let g2 = gusset(&mut m, -2.4);
    m.rebuild_adjacency();
    let wg = boolean(&mut m, BoolKind::Fuse, w, g1).expect("wall + gusset")[0];
    m.rebuild_adjacency();
    let report = arrangement::operand_vertex_audit(&m, wg, g2).expect("the audit runs");
    let four: Vec<_> = report.iter().filter(|r| r.topo.len() >= 4).collect();
    for r in &four {
        eprintln!(
            "four-plane operand vertex {:?} {:?}: topo {:?} geom {:?} names {:?} dependent {:?} folded {:?} on vertex {:?}",
            r.vertex, r.point, r.topo, r.geom, r.names, r.dependent, r.folded, r.folded_on_vertex
        );
    }
    assert_eq!(four.len(), 2, "the apex edge's two ends: {four:?}");
    for r in &four {
        assert_eq!(r.side, 0, "the fused operand carries them");
        assert_eq!(
            r.geom.len(),
            4,
            "exactly four planes through the point: {:?}",
            r.geom
        );
        assert_eq!(
            r.topo, r.geom,
            "the topology already knows every plane through it"
        );
        let mut distinct: Vec<NodeId> = r.names.iter().map(|(_, n)| *n).collect();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 1, "one name for the vertex: {:?}", r.names);
        assert!(
            r.dependent.iter().all(|&d| !d),
            "no name shares a line: {:?}",
            r.names
        );
        let mut folded = r.folded.clone();
        folded.sort_unstable();
        folded.dedup();
        assert_eq!(folded.len(), 1, "one representative: {:?}", r.folded);
        assert_eq!(
            folded[0], distinct[0],
            "the fold is the name itself: {:?}",
            r.folded
        );
        assert!(
            r.folded_on_vertex.iter().all(|&on| on),
            "every fold lands on the vertex: {:?}",
            r.folded
        );
    }
    // Every other vertex of both operands is a plain three-plane corner, named once.
    for r in report.iter().filter(|r| r.topo.len() == 3) {
        let mut distinct: Vec<NodeId> = r.names.iter().map(|(_, n)| *n).collect();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 1, "{:?}: {:?}", r.point, r.names);
        assert!(r.dependent.iter().all(|&d| !d), "{:?}", r.point);
    }
}

/// ★ **The four-plane operand vertex fuses, in either operand order.** The fused wall
/// and gusset carry two vertices where four faces meet; fusing the second gusset onto that solid
/// must not depend on the order (a symptom met first — `FourPlane`, `DegenerateWitness` — in one,
/// a body in the other). Both operand orders give one valid body of volume `35 + 2 · 7.5`, every
/// result vertex has its own
/// coordinate (two names for one point would show here as two vertices at one place), and the
/// apex keeps its four faces.
#[test]
fn a_four_plane_operand_vertex_fuses_in_either_order() {
    for swapped in [false, true] {
        let (mut m, wg, g2) = wall_and_gusset_operand();
        let (a, b) = if swapped { (g2, wg) } else { (wg, g2) };
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the four-plane operand fuses");
        assert_eq!(out.len(), 1, "one body (swapped {swapped})");
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - 50.0).abs() < 1e-9, "volume {v} (swapped {swapped})");
        // Every result vertex has its own coordinate.
        let sol = m.solid(out[0]);
        let mut points: Vec<[u64; 3]> = Vec::new();
        let mut apex_faces = 0usize;
        for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
            for fh in m.shell(sh).faces.clone() {
                let face = m.face(fh);
                let mut seen: Vec<Handle<Vertex>> = Vec::new();
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        for &vh in m.edge(he.edge).vertices.iter() {
                            if seen.contains(&vh) {
                                continue;
                            }
                            seen.push(vh);
                            let p = m.vertex_point(vh);
                            if (p[0] - 1.4).abs() < 1e-9
                                && (p[1] - 3.0).abs() < 1e-9
                                && (p[2] - 6.0).abs() < 1e-9
                            {
                                apex_faces += 1;
                            }
                            let bits = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
                            if !points.contains(&bits) {
                                points.push(bits);
                            }
                        }
                    }
                }
            }
        }
        let mut vertices: Vec<Handle<Vertex>> = Vec::new();
        for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
            for fh in m.shell(sh).faces.clone() {
                let face = m.face(fh);
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        for &vh in m.edge(he.edge).vertices.iter() {
                            if !vertices.contains(&vh) {
                                vertices.push(vh);
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(
            points.len(),
            vertices.len(),
            "one vertex per coordinate (swapped {swapped})"
        );
        assert_eq!(
            apex_faces, 4,
            "the apex keeps its four faces (swapped {swapped})"
        );
    }
}

/// The user's four parts — the filleted, bored plate, the slot plate standing on it, and two
/// gussets at `gx` and `−(gx + 1)` — fused in both orders (slot plate second, slot plate last):
/// one valid body of volume `48 + (34 − π/4) + 15` each time, and the two results equal to the
/// bit. Locked at `gx = 1.4` and at the script's own `1.5`.
fn users_four_part_fold_in_either_order(gx: f64) {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let prism = |m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
        let profile = stated(edges).unwrap().remove(0);
        let frame = SketchFrame::world(m, axis);
        let OpOutput::Extrude { solid, .. } = apply(
            m,
            &Operation::Extrude {
                frame,
                profile,
                dist,
            },
        )
        .expect("the profile extrudes") else {
            unreachable!()
        };
        solid
    };
    let shift = |m: &mut Model, s: Handle<Solid>, t: [f64; 3]| -> Handle<Solid> {
        let r = |x: f64| nacre_exact::Rat::from_decimal(x).unwrap();
        let OpOutput::Transform { solid } = apply(
            m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::translation([r(t[0]), r(t[1]), r(t[2])]),
            },
        )
        .expect("the translation applies") else {
            unreachable!()
        };
        solid
    };
    let plate = |m: &mut Model| {
        prism(
            m,
            Axis::Z,
            vec![
                line(p2(-1.5, -4.0), p2(1.5, -4.0)),
                arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                line(p2(3.5, -2.0), p2(3.5, 4.0)),
                line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                line(p2(-3.5, 4.0), p2(-3.5, -2.0)),
                arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1),
                stated::circle(p2(-1.5, -2.0), 1.0),
                stated::circle(p2(1.5, -2.0), 1.0),
            ],
            1.0,
        )
    };
    let slot_plate = |m: &mut Model| {
        let s = prism(
            m,
            Axis::Y,
            vec![
                line(p2(1.0, -3.5), p2(6.0, -3.5)),
                line(p2(6.0, -3.5), p2(6.0, 3.5)),
                line(p2(6.0, 3.5), p2(1.0, 3.5)),
                line(p2(1.0, 3.5), p2(1.0, -3.5)),
                line(p2(3.0, -0.5), p2(4.0, -0.5)),
                arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2),
                line(p2(4.0, 0.5), p2(3.0, 0.5)),
                arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2),
            ],
            1.0,
        );
        shift(m, s, [0.0, 3.0, 0.0])
    };
    let gusset = |m: &mut Model, x: f64| {
        let g = prism(
            m,
            Axis::X,
            vec![
                line(p2(0.0, 1.0), p2(3.0, 1.0)),
                line(p2(3.0, 1.0), p2(3.0, 6.0)),
                line(p2(3.0, 6.0), p2(0.0, 1.0)),
            ],
            1.0,
        );
        shift(m, g, [x, 0.0, 0.0])
    };
    let fuse = |m: &mut Model, a: Handle<Solid>, b: Handle<Solid>| -> Handle<Solid> {
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Fuse, a, b).expect("the fold's fuse builds");
        assert_eq!(out.len(), 1, "one body");
        out[0]
    };
    let pi = std::f64::consts::PI;
    let want = 48.0 + (34.0 - pi / 4.0) + 15.0;
    let mut digests = Vec::new();
    for slot_plate_last in [false, true] {
        let mut m = Model::new();
        let p1 = plate(&mut m);
        let p2s = slot_plate(&mut m);
        let g1 = gusset(&mut m, gx);
        let g2 = gusset(&mut m, -(gx + 1.0));
        let part = if slot_plate_last {
            let a = fuse(&mut m, p1, g1);
            let a = fuse(&mut m, a, g2);
            fuse(&mut m, a, p2s)
        } else {
            let a = fuse(&mut m, p1, p2s);
            let a = fuse(&mut m, a, g1);
            fuse(&mut m, a, g2)
        };
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, part).expect("props").volume;
        assert!(
            (v - want).abs() < 1e-9,
            "{v} vs {want} (slot plate last {slot_plate_last})"
        );
        digests.push(brep_digest(&m, part));
    }
    let (d0, d1) = (&digests[0], &digests[1]);
    assert_eq!(
        d0.vertex_bits, d1.vertex_bits,
        "the two fold orders realize the same vertices"
    );
    assert_eq!(d0.faces, d1.faces, "and the same faces");
    assert_eq!(d0.edges, d1.edges, "and the same edges");
    assert_eq!(
        d0.volume_bits, d1.volume_bits,
        "and the same volume, to the bit"
    );
    assert_eq!(d0.triangle_bits, d1.triangle_bits, "and the same mesh");
}

/// ★ **The user's four-part fold builds in script order**, and the two fold orders
/// agree to the bit. The gussets at `1.4` and `−2.4`: one order leaves a four-plane vertex for
/// the fourth fuse; the census measured the two results' digests identical,
/// and this locks that measurement.
#[test]
fn the_users_four_part_fold_builds_in_either_order() {
    users_four_part_fold_in_either_order(1.4);
}

/// ★ **The user's fold builds at the script's own dimensions.** The gussets at `1.5`
/// and `−2.5` stand with a side plane through each fillet's axis, so each such plane holds the
/// fillet's tangent ruling with the plate's bottom wall — one line shared by two planes and a
/// cylinder. Same volume, both
/// orders, to the bit.
#[test]
fn the_users_fold_builds_at_the_scripts_own_dimensions() {
    users_four_part_fold_in_either_order(1.5);
}

/// ★ **One name under every rigid motion.** The four-plane operand vertex is named from
/// its incident classes, which a motion permutes and renumbers; under each motion of the oracle's
/// group the moved operands must still name each apex once, never by a dependent triple, with the
/// alias fold landing on the vertex — and the moved fuse must build to the same volume. (Why this
/// is not a rigid-motion oracle row: see [`wall_and_gusset_operand`].)
#[test]
fn a_four_plane_operand_vertex_has_one_name_under_rigid_motion() {
    for (mn, iso, _) in motion_group() {
        let (mut m, wg, g2) = wall_and_gusset_operand();
        let a = transform(&mut m, wg, &iso).expect("the operand moves");
        let b = transform(&mut m, g2, &iso).expect("the gusset moves");
        m.rebuild_adjacency();
        let report = arrangement::operand_vertex_audit(&m, a, b).expect("the audit runs");
        let four: Vec<_> = report.iter().filter(|r| r.topo.len() >= 4).collect();
        assert_eq!(four.len(), 2, "{mn}: the apex edge's two ends: {four:?}");
        for r in &four {
            let mut distinct: Vec<NodeId> = r.names.iter().map(|(_, n)| *n).collect();
            distinct.sort_unstable();
            distinct.dedup();
            assert_eq!(distinct.len(), 1, "{mn}: one name: {:?}", r.names);
            assert!(r.dependent.iter().all(|&d| !d), "{mn}: {:?}", r.names);
            assert!(
                r.folded_on_vertex.iter().all(|&on| on),
                "{mn}: every fold on the vertex: {:?}",
                r.folded
            );
        }
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved operand fuses");
        assert_eq!(out.len(), 1, "{mn}: one body");
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{mn}: {:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - 50.0).abs() < 1e-9, "{mn}: volume {v}");
    }
}

/// The cell-⑫ operand: a plate with one fillet (axis at `(1.5, −2)`, tangent to the bottom wall
/// along `x = 1.5`) and a slab standing on it whose face plane `x = 1.5` holds that axis.
fn fillet_plate_and_axis_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let profile = stated(vec![
        line(p2(-3.5, -4.0), p2(1.5, -4.0)),
        arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
        line(p2(3.5, -2.0), p2(3.5, 4.0)),
        line(p2(3.5, 4.0), p2(-3.5, 4.0)),
        line(p2(-3.5, 4.0), p2(-3.5, -4.0)),
    ])
    .unwrap()
    .remove(0);
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: plate, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 1.0,
        },
    )
    .expect("the plate extrudes") else {
        unreachable!()
    };
    let slab = m.add_cuboid(
        Point3::from_array([1.5, -3.0, 1.0]),
        Point3::from_array([2.5, 0.0, 2.0]),
    );
    m.rebuild_adjacency();
    (m, plate, slab)
}

/// ★ **A plane through a fillet's axis shares the fillet's tangent ruling.** The slab's
/// face plane `x = 1.5` holds the fillet's axis and so contains the fillet's tangent ruling with
/// the plate's bottom wall; the slab stands on the plate's top cap and overlaps it nowhere. Fuse
/// is the sum, `(52 + π) + 3`; cut is the plate, `52 + π`, or the slab, `3`; common is empty —
/// in both operand orders, every result valid.
#[test]
fn a_plane_through_the_fillet_axis_shares_the_tangent_ruling() {
    let pi = std::f64::consts::PI;
    let plate_v = 52.0 + pi;
    let slab_v = 3.0;
    for swapped in [false, true] {
        for (kind, want) in [
            (BoolKind::Fuse, Some(plate_v + slab_v)),
            (BoolKind::Cut, Some(if swapped { slab_v } else { plate_v })),
            (BoolKind::Common, None),
        ] {
            let (mut m, plate, slab) = fillet_plate_and_axis_slab();
            let (a, b) = if swapped {
                (slab, plate)
            } else {
                (plate, slab)
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
                None => assert!(out.is_empty(), "{kind:?} swapped {swapped}: empty: {out:?}"),
                Some(want) => {
                    assert_eq!(out.len(), 1, "{kind:?} swapped {swapped}: one body");
                    let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                    assert!(
                        (v - want).abs() < 1e-9,
                        "{kind:?} swapped {swapped}: volume {v} vs {want}"
                    );
                }
            }
        }
    }
}

/// ★ **One name under every rigid motion.** The tangent corner's identity with the
/// class's ruling crossing is read from the moved operands' own topology and the moved class
/// table (`seed_from_operands`, `side_of`), so under each motion of the oracle's group the two
/// corners on the axis plane must still fold onto one pierce representative, no class may
/// decline, and the moved fuse must build to the same volume.
#[test]
fn a_tangent_corner_has_one_name_under_rigid_motion() {
    for (mn, iso, _) in motion_group() {
        let (mut m, plate, slab) = fillet_plate_and_axis_slab();
        let a = transform(&mut m, plate, &iso).expect("the plate moves");
        let b = transform(&mut m, slab, &iso).expect("the slab moves");
        m.rebuild_adjacency();
        let corners = arrangement::pierce_corner_audit(&m, a, b).expect("the audit runs");
        assert_eq!(
            corners.len(),
            2,
            "{mn}: the two tangent corners on the axis plane: {corners:?}"
        );
        for c in &corners {
            assert_eq!(c.candidates.len(), 4, "{mn}: {:?}", c.candidates);
            let f = &c.folded;
            assert!(
                f[0] == f[1] && (f[2] == f[0]) != (f[3] == f[0]),
                "{mn}: one representative, the other root apart: {f:?}"
            );
            assert!(
                crate::combinatorics::pierce_name(f[0]).is_some(),
                "{mn}: represented on the cylinder: {:?}",
                f[0]
            );
        }
        let declines = arrangement::trace_declines(&m, a, b).expect("the traces run");
        assert!(declines.is_empty(), "{mn}: no class declines: {declines:?}");
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved operands fuse");
        assert_eq!(out.len(), 1, "{mn}: one body");
        m.rebuild_adjacency();
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{mn}: {:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 55.0 + std::f64::consts::PI;
        assert!((v - want).abs() < 1e-9, "{mn}: volume {v} vs {want}");
    }
}

/// ★ **A point on a cylinder that a foreign class passes through has one name.** The
/// class plane `x = 1.5` holds the fillet's axis, so one of its two rulings on the fillet is the
/// fillet's own tangent ruling with the plate's bottom wall, and the tangent corner at each cap
/// is named three ways — its own `Pierce … Double`, the three-plane name of its planes with the
/// class, and the class's ruling crossing at the same point. The seed (`seed_from_operands` →
/// `Aliases::record_on_cylinder`) joins them before any trace, the representative is a pierce
/// name (the point is represented on its cylinder), and the class's *other* root stays another
/// point. Unjoined, the four names would leave the lateral's sweep declining (`Ruling`, from
/// `theta_between`'s coincidence).
#[test]
fn a_tangent_corner_on_a_plane_through_the_axis_has_one_name() {
    let (m, plate, slab) = fillet_plate_and_axis_slab();
    let corners = arrangement::pierce_corner_audit(&m, plate, slab).expect("the audit runs");
    // The two tangent corners (z = 0 and z = 1) of the ruling `x = 1.5, y = −4`, each on the
    // slab's class `x = 1.5`; the fillet's other tangent corners (`y = −2` with the right wall)
    // are on no foreign class.
    let on_axis_plane: Vec<_> = corners
        .iter()
        .filter(|c| (c.point[0] - 1.5).abs() < 1e-9 && (c.point[1] + 4.0).abs() < 1e-9)
        .collect();
    assert_eq!(on_axis_plane.len(), 2, "{corners:?}");
    // Two vertices of the plate (operand A), both tangent corners, on the slab's one class.
    assert_ne!(on_axis_plane[0].vertex, on_axis_plane[1].vertex);
    assert_eq!(on_axis_plane[0].class, on_axis_plane[1].class);
    for c in on_axis_plane {
        assert_eq!(c.side, 0, "the plate's corner");
        assert!(
            matches!(
                crate::combinatorics::pierce_name(c.corner),
                Some((_, _, nacre_topo::QuadRoot::Double))
            ),
            "a tangent corner: {:?}",
            c.corner
        );
        // Candidates: the corner, the three-plane name, and the class's two ruling crossings.
        assert_eq!(c.candidates.len(), 4, "{:?}", c.candidates);
        // The corner, its three-plane name and the class's root *at* it fold onto one
        // representative; the class's other root is another point. Which of the two roots is
        // the corner's is a matter of the class's stored normal (`Lo`/`Hi`), so it is not fixed.
        let f = &c.folded;
        assert!(
            f[0] == f[1] && (f[2] == f[0]) != (f[3] == f[0]),
            "one representative, the other root apart: {f:?}"
        );
        assert!(
            crate::combinatorics::pierce_name(f[0]).is_some(),
            "a point on a cylinder is represented on the cylinder: {:?}",
            f[0]
        );
    }
    let declines = arrangement::trace_declines(&m, plate, slab).expect("the traces run");
    assert!(declines.is_empty(), "no class declines: {declines:?}");
}
