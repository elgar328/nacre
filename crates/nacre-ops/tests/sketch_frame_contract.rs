//! **`face_sketch_frame` speaks the pad's frame, or it declines by name — never a third thing.**
//!
//! The returned `SketchFrame` is not a report to read and discard: it is a value the kernel
//! accepts back (`Operation::Extrude { frame }`, `DatumDef::Offset { frame }`), so a wrong one is
//! worse than a wrong number — a sketch built in it lands somewhere the caller did not ask for.
//! Measured before the fix this file locks: on a plain axis-aligned block the returned frame put
//! the same footprint point-symmetric to the pad's, 2.55 apart.
//!
//! The contract asserted here is the realization itself, bit for bit:
//!
//! > For every face, `face_sketch_frame` either returns a frame whose realization
//! > (`frame_plane`) is **bit-identical** to `face_plane` — the frame the pad actually
//! > sketches in — or it declines by name. There is no middle.
//!
//! ★ Bits, not a tolerance: on the faces this sweeps, both sides come from the same exact
//! roads ({0,±1} axes, rational projections, or the same flip-baked chain realization), so
//! agreement is exact when it holds and a threshold would only paper over a third derivation.
//!
//! ★★ The refusal population is pinned per flavour below. It is defined by *verification
//! failure*, not by a condition list — and the invariant-plane restatement (2026-08-17)
//! delivered the intended news exactly as this header once predicted: a plane the rotation
//! maps onto itself keeps its world statement instead of carrying the motion, those cells
//! now verify bit-for-bit, and the rotated pins measured 2 → 0.

use nacre_math::Point2;
use nacre_ops::{
    OpError, OpOutput, Operation, Profile2d, SketchFrame, apply, face_plane, face_sketch_frame,
    frame_plane,
};
use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_store::Handle;
use nacre_topo::{Model, Solid};

fn rect(a: f64, b: f64, c: f64, d: f64) -> Profile2d {
    Profile2d::polygon(vec![
        Point2::from_array([a, b]),
        Point2::from_array([c, b]),
        Point2::from_array([c, d]),
        Point2::from_array([a, d]),
    ])
    .expect("a rectangle is a fair profile")
}

fn block(m: &mut Model) -> Handle<Solid> {
    let world = SketchFrame::world(m, Axis::Z);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: world,
            profile: rect(0.0, 0.0, 2.0, 2.0),
            dist: 1.0,
        },
    )
    .expect("the block") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

fn xf(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    )
    .expect("transform") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

fn rot(axis: Axis, deg: i128) -> Isometry {
    Isometry::rotation(Rotation {
        axis,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
    })
}

/// The four flavours: how the block is placed, and how many faces must decline.
///
/// Since the invariant-plane restatement, a plane the motion fixes (the z-caps here — the
/// rotation is about their own normal, the translation slides within them) keeps its world
/// statement and no motion node, so its canonical frame verifies and **nothing declines** in
/// any flavour. The declining population that remains for `FrameNotRepresentable` is the one
/// stage 1 leaves recorded: exactly-statable-but-shifted images (a z-translation after the
/// turn), mirror chains, and second-generation moved sources — none of which these flavours
/// build. Everything must agree to the bit — including the flip=true axis-aligned faces the
/// old fallback reported point-symmetric.
fn flavours() -> Vec<(&'static str, Model, Handle<Solid>, usize)> {
    let mut out = Vec::new();

    let mut m = Model::new();
    let s = block(&mut m);
    out.push(("plain", m, s, 0));

    // 90° is exact and restated (no motion recorded), so nothing may decline here either.
    let mut m = Model::new();
    let s = block(&mut m);
    let s = xf(&mut m, s, rot(Axis::Z, 90));
    out.push(("rot90", m, s, 0));

    // 30°: the two rotation-invariant planes (z = 0, z = 1) are restated — no motion, no
    // decline; the four walls carry the chain and transcribe.
    let mut m = Model::new();
    let s = block(&mut m);
    let s = xf(&mut m, s, rot(Axis::Z, 30));
    out.push(("rot30", m, s, 0));

    // 30° then an exact in-plane translation — the walls' history records it, while the caps
    // are fixed by both motions and stay restated. (Before the restatement the caps declined
    // here, and before the transcription fix their returned frame was 1.0 off in origin.)
    let mut m = Model::new();
    let s = block(&mut m);
    let s = xf(&mut m, s, rot(Axis::Z, 30));
    let s = xf(
        &mut m,
        s,
        Isometry::translation([Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)]),
    );
    out.push(("rot30+t", m, s, 0));

    out
}

#[test]
fn the_returned_frame_realizes_to_the_pads_or_declines_by_name() {
    let mut wrong: Vec<String> = Vec::new();

    for (tag, m, s, want_declined) in flavours() {
        let mut declined = 0usize;
        for fi in 0..6 {
            let f = m.shell(m.solid(s).outer).faces[fi];
            let a = face_plane(&m, f).expect("planar");
            match face_sketch_frame(&m, f) {
                Ok(sf) => {
                    let b = frame_plane(&m, &sf).expect("a returned frame realizes");
                    let bits = |p: &nacre_ops::SketchPlane| {
                        [
                            p.origin().as_array().map(f64::to_bits),
                            p.x_axis().as_array().map(f64::to_bits),
                            p.y_axis().as_array().map(f64::to_bits),
                        ]
                    };
                    if bits(&a) != bits(&b) {
                        wrong.push(format!(
                            "{tag} face[{fi}]: returned frame realizes to o={:?} x={:?} y={:?}, \
                             the pad sketches at o={:?} x={:?} y={:?}",
                            b.origin().as_array(),
                            b.x_axis().as_array(),
                            b.y_axis().as_array(),
                            a.origin().as_array(),
                            a.x_axis().as_array(),
                            a.y_axis().as_array(),
                        ));
                    }
                }
                // Declining is an honest answer; *which* faces decline is pinned below.
                Err(OpError::FrameNotRepresentable) => declined += 1,
                Err(e) => wrong.push(format!("{tag} face[{fi}]: unexpected error {e:?}")),
            }
        }
        if declined != want_declined {
            wrong.push(format!(
                "{tag}: {declined} faces declined, pinned {want_declined} — the representable \
                 population moved (restatement landing? update the pin with the measurement)"
            ));
        }
    }

    assert!(
        wrong.is_empty(),
        "faces holding a frame that is neither the pad's nor a refusal:\n{}",
        wrong.join("\n")
    );
}

/// **The value round-trips**: a sketch built in the returned frame lands exactly where the pad
/// lands. One flip=true face (the transcription population) and one flip=false face (the
/// canonical population) — the realization comparison above could in principle miss something an
/// actual operation reads, so the operation is the final witness.
#[test]
fn an_extrude_in_the_returned_frame_lands_with_the_pad() {
    for fi in [0usize, 1] {
        // face[0] = bottom (flip=true), face[1] = top (flip=false) on the plain block.
        let mut m = Model::new();
        let s = block(&mut m);
        let f = m.shell(m.solid(s).outer).faces[fi];
        let a = face_plane(&m, f).expect("planar");
        let d = nacre_props::face_props(&m, f).expect("props").centroid - a.origin();
        let (cu, cv) = (d.dot(a.x_axis()), d.dot(a.y_axis()));
        let profile = rect(cu + 0.05, cv - 0.15, cu + 0.25, cv + 0.05);
        let n = a.x_axis().cross(a.y_axis());

        let OpOutput::PadOnFace { top_face, .. } = apply(
            &mut m,
            &Operation::PadOnFace {
                face: f,
                profile: profile.clone(),
                dist: 0.5,
            },
        )
        .expect("the pad") else {
            unreachable!()
        };
        m.rebuild_adjacency();
        let pad_fp = nacre_props::face_props(&m, top_face)
            .expect("props")
            .centroid
            - n * 0.5;

        let mut m2 = Model::new();
        let s2 = block(&mut m2);
        let f2 = m2.shell(m2.solid(s2).outer).faces[fi];
        let sf = face_sketch_frame(&m2, f2).expect("a representable face");
        let OpOutput::Extrude { solid: prism, .. } = apply(
            &mut m2,
            &Operation::Extrude {
                frame: sf,
                profile,
                dist: 0.5,
            },
        )
        .expect("the extrude in the returned frame") else {
            unreachable!()
        };
        m2.rebuild_adjacency();
        let (lo, hi) = nacre_props::bounds(&m2, prism).expect("bounds");
        // The prism runs from its base along `n` for 0.5, so the base is centre − n·0.25.
        let prism_fp = (lo + (hi - lo) * 0.5) - n * 0.25;

        assert!(
            (pad_fp - prism_fp).norm() < 1e-12,
            "face[{fi}]: pad footprint at {:?}, extrude-in-returned-frame at {:?}",
            pad_fp.as_array(),
            prism_fp.as_array()
        );
    }
}
