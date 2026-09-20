//! Locks for the public [`SketchFrame`]: a `Named` placement is a *claim*, checked exactly at
//! construction and rejected by name — never silently replaced by `Canonical` — and the value a
//! constructor accepts is letter-identical to what the operation road builds for the same words.

use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::{OpError, OpOutput, Operation, Profile2d, SketchFrame, SketchPlane, apply};
use nacre_scalar::{Axis, Rat};
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
    // The S6a terminal lock's fixture (`a_prism_on_an_axes_only_tilted_frame_takes_the_exact_road`)
    // — chosen there so the canonical name is genuinely Wide, asserted again below.
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
    let s = m.add_cuboid(
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
    // ★ Not `placement == Canonical` any more — that was an implementation detail standing in
    // for the actual contract, and it was false in substance: the canonical frame *flipped*
    // realizes point-symmetric to the axes the pad sketches this face in (the node is elided on
    // an axis-aligned face, and the pad uses the world axes derived from the outward normal).
    // The contract is the realization: the frame this returns must land, bit for bit, on the
    // plane `face_plane` reports — which is what the pad reads.
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
                isometry: nacre_scalar::Isometry::rotation(nacre_scalar::Rotation {
                    axis: Axis::Y,
                    pivot: [Rat::from_int(0); 3],
                    angle: nacre_scalar::Angle::from_deg(Rat::from_int(30)).unwrap(),
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
    apply(
        &mut m,
        &Operation::PadOnFace {
            face: top,
            profile: square(0.2, 0.8),
            dist: 0.3,
        },
    )
    .expect("pad on the tilted face");
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
