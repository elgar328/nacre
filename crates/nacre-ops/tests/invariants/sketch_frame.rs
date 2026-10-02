//! Locks for the public [`SketchFrame`]: a `Named` placement is a *claim*, checked exactly at
//! construction and rejected by name — never silently replaced by `Canonical` — and the value a
//! constructor accepts is letter-identical to what the operation road builds for the same words.

use nacre_exact::{Axis, Rat};
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{OpError, OpOutput, Operation, Profile2d, SketchFrame, SketchPlane, apply};
use nacre_topo::{FramePlacement, Model, Motion};

use crate::fixtures::datum_frame;

fn square(a: f64, b: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, a]),
        Point2::from_array([b, a]),
        Point2::from_array([b, b]),
        Point2::from_array([a, b]),
    ])
    .expect("a square is a fair profile")
}

/// The `Frame` motions the model's live faces sketch in — a wall's truth points at the node,
/// and the node's `plane` is the sketch plane's handle (not the wall's own surface).
fn frame_nodes(m: &Model) -> Vec<Motion> {
    m.live_solids()
        .iter()
        .flat_map(|&s| m.shell(m.solid(s).outer).faces.iter())
        .filter_map(|&fh| {
            let su = m.face(fh).surface;
            let motion = match m.surface(su) {
                nacre_topo::Surface::Plane {
                    motion: Some(h), ..
                } => *h,
                _ => return None,
            };
            match m.motion(motion).motion {
                node @ Motion::Frame { .. } => Some(node),
                _ => None,
            }
        })
        .collect()
}

/// ★★★ **Each bad claim gets its own name** — and the good claim is accepted with the caller's
/// words lifted, not normalized away.
#[test]
fn a_named_frame_rejects_each_bad_claim_by_name() {
    let m = Model::new();
    let z0 = m.world_plane(Axis::Z);
    let o = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
    let v = |x: f64, y: f64, z: f64| Vector3::from_array([x, y, z]);

    // An origin the decimal window cannot hold: the claim has no exact statement.
    assert_eq!(
        SketchFrame::named(&m, z0, o(1e300, 0.0, 0.0), v(1.0, 0.0, 0.0)),
        Err(OpError::FrameOutsideDecimalWindow)
    );
    // A stated origin off the stated plane.
    assert_eq!(
        SketchFrame::named(&m, z0, o(0.0, 0.0, 1.0), v(1.0, 0.0, 0.0)),
        Err(OpError::OriginNotOnPlane)
    );
    // A `ref_dir` with no in-plane part — parallel to the normal, or zero.
    assert_eq!(
        SketchFrame::named(&m, z0, o(0.0, 0.0, 0.0), v(0.0, 0.0, 3.0)),
        Err(OpError::RefDirParallelToNormal)
    );
    assert_eq!(
        SketchFrame::named(&m, z0, o(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0)),
        Err(OpError::RefDirParallelToNormal)
    );

    // The good claim: accepted, and the placement holds the caller's decimals lifted — with an
    // out-of-plane `ref_dir` component kept (its in-plane part names `+u`; projection happens
    // where the frame is built, not here).
    let f = SketchFrame::named(&m, z0, o(0.5, -0.25, 0.0), v(1.0, 1.0, 0.5)).expect("a fair claim");
    assert_eq!(f.plane(), z0);
    assert!(!f.flip(), "flip is the consuming operation's to measure");
    let d = |x: f64| Rat::from_decimal(x).unwrap();
    assert_eq!(
        *f.placement(),
        FramePlacement::Named {
            origin: [d(0.5), d(-0.25), d(0.0)],
            ref_dir: [d(1.0), d(1.0), d(0.5)],
        }
    );

    // And the frame it builds is that in-plane part: a unit square swept by `0.5` in it is the
    // unit-square prism, standing on `z = 0` with `+u` along `(1, 1, 0)`.
    let mut m = m;
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: f,
            profile: square(0.0, 1.0),
            dist: 0.5,
        },
    )
    .expect("an extrude in the accepted frame") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty(), "a valid prism");
    let volume = nacre_props::mass_props(&m, solid).expect("props").volume;
    assert!(
        (volume - 0.5).abs() < 1e-12,
        "the unit-square prism, got {volume}"
    );
    let u = nacre_ops::frame_plane(&m, &f)
        .expect("realizes")
        .x_axis()
        .as_array();
    let h = std::f64::consts::FRAC_1_SQRT_2;
    assert!(
        (u[0] - h).abs() < 1e-12 && (u[1] - h).abs() < 1e-12 && u[2] == 0.0,
        "+u is the in-plane part of (1, 1, 0.5), got {u:?}"
    );
}

/// ★★ **The constructor and the extrude road speak one value** — restating the axes-only plane's
/// own origin and `+u` through [`SketchFrame::named`] reproduces the operation's placement
/// letter for letter. The fixture's name is *wide* (full-width decimal axes: the cross-product
/// denominators pass `i128`), so the accepted origin also proves the residual check's
/// arbitrary-precision arm answers `0` on the plane — and the off-plane probe proves it answers
/// nonzero off it.
#[test]
fn a_named_frame_states_what_the_extrude_road_builds() {
    let mut m = Model::new();
    // The fixture of `a_prism_on_an_axes_only_tilted_frame_takes_the_exact_road` — chosen there
    // so the canonical name is genuinely Wide, asserted again below.
    let (origin, u, v) = (
        Point3::from_array([0.2547863291057384, -0.5123456789012345, 1.5432109876543211]),
        Vector3::from_array([0.7123456789012345, 0.5876543210987654, 0.4098765432101234]),
        Vector3::from_array([-0.5876543210987654, 0.7123456789012345, 0.1234567890123456]),
    );
    let plane = SketchPlane::from_axes(origin, u, v);
    let __f107 = datum_frame(&mut m, plane);
    let out = apply(
        &mut m,
        &Operation::Extrude {
            frame: __f107,
            profile: square(0.1, 1.6),
            dist: 0.7,
        },
    )
    .expect("the tilted extrude");
    let OpOutput::Extrude { .. } = out else {
        unreachable!()
    };
    let nodes = frame_nodes(&m);
    let node = nodes.first().expect("the prism sketched in a frame");
    let Motion::Frame {
        plane: plane_h,
        placement,
        ..
    } = node
    else {
        unreachable!()
    };
    assert!(
        m.surface_name.get(plane_h).unwrap().narrow().is_none(),
        "the fixture must actually be wide, or the arbitrary-precision arm went untested"
    );

    let f = SketchFrame::named(&m, *plane_h, origin, u).expect("the plane's own words");
    assert_eq!(
        f.placement(),
        placement,
        "one claim, one value — the constructor and the road must not lift differently"
    );

    // The same wide name must also *reject* exactly: the plane's own origin shifted 1e-9 in z —
    // far below any geometric feature, visible only to an exact residual. (Not one decimal digit:
    // near 1.5 the f64 ulp is ~2.2e-16, so a 1e-16 shift rounds back to the same double and the
    // probe would claim the very point it means to miss.)
    let off = Point3::from_array([0.2547863291057384, -0.5123456789012345, 1.543210988654321]);
    assert_ne!(off.as_array()[2], origin.as_array()[2], "a real shift");
    assert_eq!(
        SketchFrame::named(&m, *plane_h, off, u),
        Err(OpError::OriginNotOnPlane)
    );
}

/// ★★ **`face_sketch_frame` reports the very frame the pad sketches in** — same plane handle,
/// same placement, same measured flip as the `Frame` node the pad leaves behind. And on an
/// axis-aligned face, where the operation elides the node, the canonical frame it reports still
/// tells the truth about direction: the top face (outward `+ẑ`, canonical `ŵ = +ẑ`) is unflipped,
/// the seeded bottom (outward `−ẑ`) is flipped.
#[test]
fn face_sketch_frame_reports_the_frame_the_pad_uses() {
    let mut m = Model::new();
    let s = nacre_ops::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 2.0, 1.0]),
    );
    m.rebuild_adjacency();

    // Axis-aligned faces first: canonical placement, hand-checkable flips.
    let face_with_surface = |m: &Model, su| {
        m.shell(m.solid(s).outer)
            .faces
            .iter()
            .copied()
            .find(|&fh| m.face(fh).surface == su)
            .expect("the cuboid face")
    };
    let bottom = face_with_surface(&m, m.world_plane(Axis::Z));
    let fb = nacre_ops::face_sketch_frame(&m, bottom).expect("the seeded bottom");
    assert_eq!(fb.plane(), m.world_plane(Axis::Z));
    assert!(fb.flip(), "outward −ẑ against canonical ŵ = +ẑ");
    // ★ Not `placement == Canonical`: which placement spells the frame is an implementation
    // detail (the normal form), and it can hold while the frame is wrong — the contract is the
    // realization: the frame this returns must land, bit for bit, on the plane `face_plane`
    // reports — which is what the pad reads.
    // `tests/invariants/sketch_frame_contract.rs`
    // sweeps this same proposition over four placements × six faces.
    let realized = nacre_ops::frame_plane(&m, &fb).expect("the returned frame realizes");
    let reported = nacre_ops::face_plane(&m, bottom).expect("planar");
    let bits = |p: &nacre_ops::SketchPlane| {
        [
            p.origin().as_array().map(f64::to_bits),
            p.x_axis().as_array().map(f64::to_bits),
            p.y_axis().as_array().map(f64::to_bits),
        ]
    };
    assert_eq!(
        bits(&realized),
        bits(&reported),
        "the returned frame realizes somewhere the pad does not sketch"
    );

    // A tilted face: turn the box, pad on its (rotated) top, then compare the reported frame
    // with the node the pad actually sketched in.
    let turned = {
        let OpOutput::Transform { solid } = apply(
            &mut m,
            &Operation::Transform {
                solid: s,
                isometry: nacre_exact::Isometry::rotation(nacre_exact::Rotation {
                    axis: Axis::Y,
                    pivot: [Rat::from_int(0); 3],
                    angle: nacre_exact::Angle::from_deg(Rat::from_int(30)).unwrap(),
                }),
            },
        )
        .expect("rotate") else {
            unreachable!()
        };
        solid
    };
    m.rebuild_adjacency();
    let top = m
        .shell(m.solid(turned).outer)
        .faces
        .iter()
        .copied()
        .max_by(|&a, &b| {
            let z = |fh| {
                let f = m.face(fh);
                let p = match m.surface_cache(f.surface) {
                    nacre_geom::Surface::Plane(p) => *p,
                    _ => unreachable!(),
                };
                let sgn = f64::from(f.orientation.sign());
                p.normal().as_array()[2] * sgn
            };
            z(a).partial_cmp(&z(b)).unwrap()
        })
        .expect("a most-upward face");
    let reported = nacre_ops::face_sketch_frame(&m, top).expect("a planar live face");
    nacre_ops::fixtures::pad(&mut m, top, square(0.2, 0.8), 0.3).expect("pad on the tilted face");
    let used: Vec<_> = frame_nodes(&m)
        .into_iter()
        .filter(|n| matches!(n, Motion::Frame { plane, .. } if *plane == reported.plane()))
        .collect();
    let node = used.first().expect("the pad sketched in the face's frame");
    let Motion::Frame {
        plane,
        placement,
        flip,
    } = node
    else {
        panic!("a Frame node")
    };
    assert_eq!(*plane, reported.plane());
    assert_eq!(placement, reported.placement());
    assert_eq!(*flip, reported.flip());
}

// ─────────────────────────────────────────────────────────────────────────────
// A face's frame is chosen from the truth — the road, the origin and the axes
// ─────────────────────────────────────────────────────────────────────────────

/// A square `2·half` across, centred on `face`'s centroid in the frame the pad sketches in.
fn centred_square(m: &Model, face: nacre_store::Handle<nacre_topo::Face>, half: f64) -> Profile2d {
    let at = nacre_ops::face_plane(m, face).expect("planar");
    let d = nacre_props::face_props(m, face).expect("props").centroid - at.origin();
    let (cu, cv) = (d.dot(at.x_axis()), d.dot(at.y_axis()));
    Profile2d::polygon(vec![
        Point2::from_array([cu - half, cv - half]),
        Point2::from_array([cu + half, cv - half]),
        Point2::from_array([cu + half, cv + half]),
        Point2::from_array([cu - half, cv + half]),
    ])
    .expect("a square")
}

/// Pad `face` with a centred square and return the pad's top face.
fn pad_centred(
    m: &mut Model,
    face: nacre_store::Handle<nacre_topo::Face>,
    dist: f64,
) -> nacre_store::Handle<nacre_topo::Face> {
    let profile = centred_square(m, face, 0.2);
    let top_face = nacre_ops::fixtures::pad(m, face, profile, dist)
        .expect("the pad")
        .cap_face(m);
    m.rebuild_adjacency();
    top_face
}

/// Extrude `profile` by `dist` in `frame` and return the far cap's surface.
fn far_cap(
    m: &mut Model,
    frame: SketchFrame,
    profile: Profile2d,
    dist: f64,
) -> nacre_store::Handle<nacre_topo::Surface> {
    let OpOutput::Extrude { faces, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("the extrude") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    m.face(faces[1]).surface
}

/// A way to lift a block, and whether the lift is recorded as a node.
type Lift = fn(&mut Model, f64, Rat) -> nacre_store::Handle<nacre_topo::Solid>;

/// ★★ **A pad on a lifted face lands its top on the exact plane.**
///
/// A block `h` tall lifted by `1/3` has its top on `z = h + 1/3`, a rational with no short
/// decimal. The pad's frame origin used to be that plane's cache realized and lifted back with
/// `Rat::from_decimal` — the decimal the `f64` prints — so the pad's far cap stood on
/// `28666666666666667/2·10¹⁶` instead of `43/30` (measured, 6/6 lifts). Two roads to that top:
/// the plain block, which carries the lift into its statements, and the `Through` block, whose
/// every motion is a node. The oracle is the transform's own: a block `h + 1` tall lifted the
/// same way.
#[test]
fn a_pad_on_a_lifted_face_lands_on_the_exact_plane() {
    let roads: [(&str, Lift, bool); 2] = [
        ("plain", crate::fixtures::lifted_block, false),
        ("through", crate::fixtures::through_lifted_block, true),
    ];
    for (road, lifted, recorded) in roads {
        for (h, (p, q)) in [(0.1, (1, 3)), (0.3, (2, 7))] {
            let lift = Rat::new(p, q).expect("a lift");
            let mut m = Model::new();
            let block = lifted(&mut m, h, lift);
            let top = crate::fixtures::top_face(&m, block);
            assert_eq!(
                m.plane_motion(m.face(top).surface).is_some(),
                recorded,
                "{road}: whether the lift is recorded"
            );
            let padded = pad_centred(&mut m, top, 1.0);
            let mut want = Model::new();
            let want_block = lifted(&mut want, h + 1.0, lift);
            assert_eq!(
                m.world_plane_name(m.face(padded).surface),
                want.world_plane_name(crate::fixtures::top_surface(&want, want_block)),
                "{road}, h = {h}, lift = {p}/{q}: the pad's top is the lifted top of a block one \
                 taller"
            );
        }
    }
}

/// ★★ **A pad on a rational tilted face stays on the world road, so its top is the plane every
/// other statement of it names.**
///
/// The far cap of a prism on a `(3, 4, 0)` datum has rational axes, but its cached normal
/// realizes to `0.6000000000000001`, and the pad used to ask whether those realized axes lift back
/// to orthonormal rationals — they do not, so it sketched inside a frame node and its top landed
/// under a `(name, node)` key: one plane, two handles. The frame is asked of the name now.
#[test]
fn a_pad_on_a_rational_tilted_face_shares_the_plane_with_an_extrude() {
    let mut m = Model::new();
    let frame = datum_frame(
        &mut m,
        SketchPlane::from_origin_normal(
            Point3::from_array([0.0; 3]),
            Vector3::from_array([3.0, 4.0, 0.0]),
        )
        .expect("a plane"),
    );
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 1.0),
            dist: 1.0,
        },
    )
    .expect("the prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let padded = pad_centred(&mut m, faces[1], 1.0);
    let far = far_cap(&mut m, frame, square(0.0, 1.0), 2.0);
    assert_eq!(
        m.face(padded).surface,
        far,
        "the pad's top and the two-high prism's far cap are one plane, so one handle"
    );
}

/// ★★ **A plane whose origin has no short decimal is framed exactly.**
///
/// `3x + 4y + 12z = 12` through three corners of a `4 × 3 × 1` box has rational axes
/// (`(−4, 3, 0)/5`, a `v̂` over `65`) and its world origin `12/169·(3, 4, 12)` has no short
/// decimal — the realization the pad used to lift back could not be the plane's own point. A pad
/// on a prism raised off it must land where the prism twice as high does.
#[test]
fn a_pad_on_a_plane_with_a_long_origin_shares_the_plane_with_an_extrude() {
    let mut m = Model::new();
    let rect = Profile2d::polygon(vec![
        Point2::from_array([0.0, 0.0]),
        Point2::from_array([4.0, 0.0]),
        Point2::from_array([4.0, 3.0]),
        Point2::from_array([0.0, 3.0]),
    ])
    .expect("a rectangle");
    let world = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: world,
            profile: rect,
            dist: 1.0,
        },
    )
    .expect("a box") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let corner = |want: [f64; 3]| {
        m.shell(m.solid(solid).outer)
            .faces
            .iter()
            .flat_map(|&f| m.face(f).outer.half_edges.clone())
            .flat_map(|he| m.edge(he.edge).vertices)
            .find(|&v| m.vertex_point(v).as_array() == want)
            .expect("a box corner")
    };
    let vs = [
        corner([4.0, 0.0, 0.0]),
        corner([0.0, 3.0, 0.0]),
        corner([0.0, 0.0, 1.0]),
    ];
    let OpOutput::DatumPlane { frame, .. } = apply(
        &mut m,
        &Operation::DatumPlane {
            def: nacre_ops::DatumDef::ThroughVertices(vs),
        },
    )
    .expect("the datum") else {
        unreachable!()
    };
    let OpOutput::Extrude { faces, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square(0.0, 1.0),
            dist: 1.0,
        },
    )
    .expect("the prism") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    let padded = pad_centred(&mut m, faces[1], 1.0);
    let far = far_cap(&mut m, frame, square(0.0, 1.0), 2.0);
    assert_eq!(m.face(padded).surface, far, "one plane, one handle");
}

/// ★★ **The frame a face reports is the frame the pad builds in — on a lifted face and on a
/// mirrored one.**
///
/// `face_sketch_frame` hands a caller the frame to extrude in; an extrude there must land where
/// the pad lands, as the same handle. On the `Through` block, whose every motion is a node, a
/// lifted face's frame is carried out through the face's own chain, and a mirrored one's chain
/// has a reflection in it, so the carried frame is left-handed; the plain block carries its lift
/// and mirror into its statements and takes the world's frame. The pad and the extrude must agree
/// on every one, and the solid stay right side out.
#[test]
fn the_reported_frame_extrudes_where_the_pad_lands() {
    let lift = Rat::new(1, 3).expect("a lift");
    let blocks: [(&str, Lift); 4] = [
        ("plain", crate::fixtures::lifted_block),
        ("plain mirrored", crate::fixtures::mirrored_lifted_block),
        ("through", crate::fixtures::through_lifted_block),
        ("through mirrored", crate::fixtures::through_mirrored_block),
    ];
    for (what, make) in blocks {
        let mut m = Model::new();
        let block = make(&mut m, 0.1, lift);
        let top = crate::fixtures::top_face(&m, block);
        let frame = nacre_ops::face_sketch_frame(&m, top).expect("a lifted face has a frame");
        let profile = centred_square(&m, top, 0.2);
        let feat = nacre_ops::fixtures::pad(&mut m, top, profile.clone(), 1.0).expect("the pad");
        let (solid, top_face) = (feat.solid(), feat.cap_face(&m));
        m.rebuild_adjacency();
        let bad = nacre_validate::validate(&m);
        assert!(bad.is_empty(), "{what}: {bad:?}");
        let v = crate::fixtures::volume(&m, solid);
        assert!((v - (4.0 * 0.1 + 0.16)).abs() < 1e-12, "{what}: volume {v}");
        let far = far_cap(&mut m, frame, profile, 1.0);
        assert_eq!(
            m.face(top_face).surface,
            far,
            "{what}: the pad's top and the extrude's far cap are one handle"
        );
    }
}
