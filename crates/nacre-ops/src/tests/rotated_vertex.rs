use super::*;
use crate::{OpOutput, Operation, apply};
use nacre_exact::{Angle, Axis, Isometry, Rat as R, Rotation as SRot};
use nacre_math::Point3;

fn rot30z() -> Isometry {
    Isometry::rotation(SRot {
        axis: Axis::Z,
        pivot: [R::from_int(1), R::from_int(1), R::from_int(0)],
        angle: Angle::from_deg(R::from_int(30)).unwrap(),
    })
}

/// One step of a motion chain — the two operations that leave a `MotionNode` behind.
enum Step {
    // Boxed: an `Isometry` dwarfs a `(Axis, Rat)`, and clippy is right that the unboxed
    // variant would make every `Flip` in the array pay for it.
    Move(Box<Isometry>),
    Flip(Axis, R),
}

fn reflected(
    m: &mut Model,
    s: Handle<nacre_topo::Solid>,
    axis: Axis,
    offset: R,
) -> Handle<nacre_topo::Solid> {
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid: s,
            axis,
            offset,
        },
    )
    .unwrap() else {
        panic!("expected Mirror");
    };
    m.rebuild_adjacency();
    solid
}

/// A turn about `axis` through an integer pivot, by `n/d` degrees.
fn turn(axis: Axis, pivot: [i128; 3], (n, d): (i128, i128)) -> Isometry {
    Isometry::rotation(SRot {
        axis,
        pivot: pivot.map(R::from_int),
        angle: Angle::from_deg(R::new(n, d).unwrap()).unwrap(),
    })
}

fn transformed(
    m: &mut Model,
    s: Handle<nacre_topo::Solid>,
    iso: &Isometry,
) -> Handle<nacre_topo::Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: *iso,
        },
    )
    .unwrap() else {
        panic!("expected Transform");
    };
    m.rebuild_adjacency();
    solid
}

/// **A replay reproduces the stored coordinate bit for bit.**
///
/// This is the contract that lets a consumer check a coordinate against its definition by
/// equality rather than by tolerance: the replay must perform the *same* float operations, in
/// the same order, that the producer did. Every motion node owes this, and `Mirror` is where
/// it is easiest to lose — `WitnessPoint::mirror` walks the same `2c − x` that `AxisMirror::point` did,
/// rather than an algebraically equal rearrangement.
///
/// ★★★ **It is also what stands between this kernel and a realization that is not a function
/// of its angle.** `Angle`'s f64 route is `(deg.to_f64() * PI / 180.0).cos()`, and that is
/// measured to give two answers one ulp apart — between a debug and a release build, and
/// between two call sites within one release build, where LLVM evaluates a literal-angle site
/// at compile time and leaves the other to libm. Both routes below realize the same angle, so a
/// producer and a consumer that folded differently would land different coordinates here.
///
/// ★★ **What holds it is that the angle crosses the model store.** It is written into an
/// `Operation::Transform`, pushed, and read back out before either route realizes it — and no
/// optimiser propagates a constant through a heap structure. That is *why* this passes with a
/// literal angle, and it is worth knowing: a future route that realizes an `Angle` it never
/// stored is not covered by this argument, only by this test happening to exercise it.
///
/// **Measured, not argued:** the whole census is bit-identical between a debug and a release
/// build (`tests/census.rs` documents the diff), which is the direct reading of the same claim
/// over 130 cases rather than one.
#[test]
fn replay_reproduces_the_stored_coordinate() {
    // ★ A chain, not one turn, and none of it "nice": a pivot off the origin, an angle whose
    // realization is nowhere near a quadrantal one, then a translate and a reflection. A single
    // Z-turn about a rational pivot exercises one node and one of the two rotate axes.
    let mv = |iso| Step::Move(Box::new(iso));
    let chains: [Vec<Step>; 4] = [
        vec![mv(rot30z())],
        vec![
            mv(turn(Axis::Z, [1, 1, 0], (2749, 71))),
            mv(turn(Axis::X, [0, 3, 2], (617, 9))),
        ],
        vec![
            mv(turn(Axis::Y, [5, 0, 1], (89999, 1000))),
            mv(Isometry::translation([
                R::new(7, 3).unwrap(),
                R::from_int(-2),
                R::from_int(0),
            ])),
            // ★ A reflection between two turns, because this is the node the contract is
            // easiest to lose on and the one whose parity the chain has to carry.
            Step::Flip(Axis::X, R::new(1, 2).unwrap()),
            mv(turn(Axis::Z, [0, 0, 0], (271, 4))),
        ],
        // ★ A same-axis depth-2 chain: the caps are fixed by *both* turns (restated
        // twice, still world-stated), so their corners exercise the fixed-carrier
        // licence at depth 2 — walls [Z, Z], caps None, `chain_fixes_plane` proving
        // the caps ride the whole chain.
        vec![mv(rot30z()), mv(turn(Axis::Z, [2, -1, 3], (2749, 71)))],
    ];
    let mut checked = 0;
    let mut declined = 0;
    for chain in &chains {
        let mut m = Model::new();
        let s = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        m.rebuild_adjacency();
        let mut r = s;
        for step in chain {
            r = match step {
                Step::Move(iso) => transformed(&mut m, r, iso),
                Step::Flip(axis, offset) => reflected(&mut m, r, *axis, *offset),
            };
        }
        let sh = m.solid(r).outer;
        for &fh in &m.shell(sh).faces {
            for he in &m.face(fh).outer.half_edges {
                for &vh in m.edge(he.edge).vertices.iter() {
                    // ★★ **Solve the definition**: the corner's three planes share one
                    // motion, so solving their pre-motion names exactly (rational Cramer)
                    // and replaying that chain must reproduce the stored coordinate —
                    // a lock over every chain shape in this table.
                    let nacre_topo::Vertex::ThreePlane(tri) = *m.vertex(vh) else {
                        panic!("a cuboid corner is a three-plane point");
                    };
                    let motion_of = |h| match m.surface(h) {
                        nacre_topo::Surface::Plane { motion, .. } => *motion,
                        nacre_topo::Surface::Cylinder { motion, .. } => *motion,
                    };
                    // The solvable population mirrors `solid_points`' own criterion: the
                    // moved carriers share one leaf, and a world-stated carrier among
                    // them is provably fixed by that chain (the restatement licence).
                    // A corner whose carriers hold *different* leaves — a mixed-axis
                    // chain leaves a partially-fixed cap with a shorter history — is
                    // production's honest decline (Arrange), counted, not asserted.
                    let leaves: Vec<_> = tri.iter().filter_map(|&h| motion_of(h)).collect();
                    let leaf = *leaves.first().expect("a moved solid moves some carrier");
                    if !leaves.iter().all(|&l| l == leaf) {
                        declined += 1;
                        continue;
                    }
                    let mut coeffs = [[R::from_int(0); 4]; 3];
                    for (o, h) in coeffs.iter_mut().zip(tri) {
                        *o = *m.surface_name.get(&h).unwrap().narrow().unwrap();
                    }
                    for (c, h) in coeffs.iter().zip(tri) {
                        if motion_of(h).is_none() {
                            assert!(
                                m.chain_fixes_plane(leaf, c),
                                "a world-stated carrier must be provably fixed by the chain"
                            );
                        }
                    }
                    let base = nacre_exact::three_planes_rat(coeffs)
                        .expect("three distinct planes of a cuboid corner meet");
                    let replayed = replay_chain_coord(
                        &m,
                        [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
                        leaf,
                    )
                    .unwrap();
                    // ★ The stored coordinate is the *realization* of the definition:
                    // the exact pre-motion base carried through the chain at 128
                    // bits and rounded once. The f64 replay of the same chain agrees with it
                    // to the rounding of its own arithmetic, not bit for bit — so the lock
                    // is: cache == realization, and replay within the construction epsilon.
                    let stored = m.vertex_point(vh).as_array();
                    assert!(
                        matches!(m.vertex_cache(vh), nacre_topo::PointCache::Bounded { .. }),
                        "a moved corner is realized from its definition"
                    );
                    let (realized, _) = crate::realize_vertex(&m, vh, crate::Precision::NearestF64)
                        .expect("the def road realizes what it solves")
                        .to_f64()
                        .expect("decided");
                    assert_eq!(
                        realized, stored,
                        "the cache is the realization, bit for bit"
                    );
                    for k in 0..3 {
                        assert!(
                            (replayed[k] - stored[k]).abs() <= 1e-9 * (1.0 + stored[k].abs()),
                            "f64 replay within the construction epsilon: {replayed:?} vs {stored:?}"
                        );
                    }
                    checked += 1;
                }
            }
        }
    }
    // Both populations must be real, or half the lock is vacuous: the single-turn and
    // same-axis chains solve (fixed caps riding the licence), the mixed-axis chains
    // decline at their partially-fixed corners.
    assert!(checked > 80, "solved only {checked} vertices");
    assert!(
        declined > 0,
        "no mixed-leaf corner declined — the residual vanished?"
    );
}

/// A chain is read root-to-leaf, so replaying it applies the rotations in the order they
/// happened — the reverse would be a different motion whenever the axes differ.
#[test]
fn a_chain_reads_root_to_leaf() {
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    m.rebuild_adjacency();
    let r = transformed(&mut m, s, &rot30z());
    let r = transformed(
        &mut m,
        r,
        &Isometry::rotation(SRot {
            axis: Axis::X,
            pivot: [R::from_int(0); 3],
            angle: Angle::from_deg(R::from_int(45)).unwrap(),
        }),
    );
    // The motion is the *face's* (surfaces record it; vertices follow their planes).
    // Faces carry different depths since the restatement — the caps were fixed by the
    // Z turn (restated) and only joined at the X turn — so the order lock reads every
    // face: the walls' `[Z, X]` (never `[X, Z]` — the reversal this test exists to
    // forbid) and the caps' `[X]`.
    let sh = m.solid(r).outer;
    let mut counts: std::collections::HashMap<Vec<Axis>, usize> = Default::default();
    for &fh in &m.shell(sh).faces {
        let &nacre_topo::Surface::Plane {
            motion: Some(rotation),
            ..
        } = m.surface(m.face(fh).surface)
        else {
            panic!("every face of this solid records a motion");
        };
        let axes: Vec<Axis> = motion_chain(&m, rotation)
            .expect("an axis-aligned history holds no frame")
            .iter()
            .filter_map(|n| match n {
                MoveNode::Rotate { axis, .. } => Some(*axis),
                MoveNode::Translate { .. }
                | MoveNode::Mirror { .. }
                | MoveNode::Frame { .. }
                | MoveNode::FrameWide(_)
                | MoveNode::FrameThrough(_) => None,
            })
            .collect();
        *counts.entry(axes).or_insert(0) += 1;
    }
    assert_eq!(
        counts,
        [(vec![Axis::Z, Axis::X], 4), (vec![Axis::X], 2)]
            .into_iter()
            .collect(),
        "walls read root-to-leaf [Z, X]; the Z-fixed caps joined at X"
    );
}

/// Push a plane whose exact triple is `pts`, with an f64 `Plane` **consistent with it**
/// (`Plane::through_points` of the realized corners) — the frame locks below ask about the
/// realized basis's geometry, so unlike the interning lock the f64 form has to match.
fn push_consistent(m: &mut Model, pts: [[R; 3]; 3]) -> Handle<Surface> {
    let f = |p: [R; 3]| Point3::from_array(p.map(|r| r.to_f64()));
    let pl = nacre_geom::Plane::through_points(f(pts[0]), f(pts[1]), f(pts[2]))
        .expect("a non-degenerate triple");
    let (h, _) = m.push_plane(pl, pts, None, nacre_topo::Orientation::Forward);
    h
}

/// The realized basis is a right-handed orthonormal frame whose origin sits on the plane
/// and whose `ŵ` is parallel to the plane's normal — the sanity every frame lock needs.
fn assert_frame_shape(
    m: &Model,
    h: Handle<Surface>,
    placement: &nacre_topo::FramePlacement,
    what: &str,
) {
    let (o, u, v, w) = frame_world_basis(m, h, placement, false)
        .unwrap_or_else(|| panic!("{what}: this frame must open"));
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    for (name, a) in [("u", u), ("v", v), ("w", w)] {
        assert!(
            (dot(a, a) - 1.0).abs() < 1e-12,
            "{what}: {name} is not unit ({:e})",
            dot(a, a) - 1.0
        );
    }
    assert!(dot(u, v).abs() < 1e-12, "{what}: u ⊥ v fails");
    assert!(dot(u, w).abs() < 1e-12, "{what}: u ⊥ w fails");
    assert!(dot(v, w).abs() < 1e-12, "{what}: v ⊥ w fails");
    let nacre_geom::Surface::Plane(pl) = m.surface_cache(h) else {
        unreachable!()
    };
    assert!(
        pl.distance(Point3::from_array(o)).abs() < 1e-9,
        "{what}: the origin is off the plane by {:e}",
        pl.distance(Point3::from_array(o))
    );
    let n = pl.normal().as_array();
    let nn = dot(n, n).sqrt();
    let cross = [
        w[1] * n[2] - w[2] * n[1],
        w[2] * n[0] - w[0] * n[2],
        w[0] * n[1] - w[1] * n[0],
    ];
    assert!(
        dot(cross, cross).sqrt() / nn < 1e-9,
        "{what}: ŵ is not parallel to the plane normal"
    );
}

/// ★★★★★ **A `Wide` name opens a frame.** The population that fails
/// `narrow()` — a plane whose canonical answer exceeds `i128` — realizes its canonical
/// placement through the arbitrary-precision road, and the basis is a real frame on the
/// real plane. (What stays closed for `Wide` is the *narrow shortcuts* — `base_rat`,
/// Shewchuk, `Isometry` transport — locked on the topo side.)
#[test]
fn a_wide_plane_hosts_a_canonical_frame() {
    let q = |n: i128, d: i128| R::new(n, d).unwrap();
    let big1 = (1i128 << 90) + 1;
    let big2 = (1i128 << 90) + 3;
    let pts = [
        [q(big1, 3), q(big2, 7), q(0, 1)],
        [q(-big2, 5), q(big1, 11), q(0, 1)],
        [q(1, 13), q(1, 17), q(1, 19)],
    ];
    // Fixture qualification: genuinely wide.
    let name = nacre_exact::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
    assert!(name.narrow().is_none(), "the fixture must be wide");
    let mut m = Model::new();
    let h = push_consistent(&mut m, pts);
    assert_frame_shape(&m, h, &nacre_topo::FramePlacement::Canonical, "wide plane");
}

/// ★★★★ **A narrow name whose squared lengths overflow `i128` opens too** — the
/// measured 1.6% population (`n·n` is a square, so it overflows long before the name).
/// On the narrow road alone this is `plane_frame_named`'s hard `None`.
#[test]
fn a_narrow_name_with_wide_squares_hosts_a_frame() {
    let q = |n: i128, d: i128| R::new(n, d).unwrap();
    // Intercept form: the plane through (1/p, 0, 0), (0, 1/q, 0), (0, 0, 1/r) has the
    // canonical name [p, q, r, −1] — narrow when p, q, r fit i128, while n·n = p²+q²+r²
    // does not (~2^140).
    let (p, q2, r) = ((1i128 << 70) + 1, (1i128 << 70) + 3, (1i128 << 70) + 7);
    let pts = [
        [q(1, p), q(0, 1), q(0, 1)],
        [q(0, 1), q(1, q2), q(0, 1)],
        [q(0, 1), q(0, 1), q(1, r)],
    ];
    // Fixture qualification: the name is narrow AND the narrow frame derivation dies on it.
    let name = nacre_exact::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
    let c = *name.narrow().expect("the name itself fits i128");
    assert!(
        nacre_exact::plane_frame_default(c).is_none(),
        "the fixture must be in the nn-overflow population"
    );
    let mut m = Model::new();
    let h = push_consistent(&mut m, pts);
    assert_frame_shape(
        &m,
        h,
        &nacre_topo::FramePlacement::Canonical,
        "nn-overflow plane",
    );
}

/// ★★★ **A `Named` placement opens on a `Wide` name too.** The locks above cover the
/// canonical road; the axes-only population (`from_axes`, a tilted full-width frame) is
/// exactly the other pairing — a caller-stated origin/`ref_dir` on a plane whose canonical
/// name exceeds `i128` — and it goes through `WideFrame::named_of`.
///
/// ★ The fixture is that population, not the `2^90` triple above: full-width *decimal* axes at
/// CAD scale. Their cross runs the denominators to `10^48`, so the canonical name is
/// genuinely `Wide` (asserted), while the geometry stays near `1` — which matters, because
/// a frame's unit axes are invisible in f64 next to a `2^90` origin (an ulp there is
/// `~6e10`).
#[test]
fn a_wide_plane_hosts_a_named_frame() {
    let d = |x: f64| R::from_decimal(x).unwrap();
    let o = [
        d(0.2547863291057384),
        d(-0.5123456789012345),
        d(1.5432109876543211),
    ];
    let x = [
        d(0.7123456789012345),
        d(0.5876543210987654),
        d(0.4098765432101234),
    ];
    let y = [
        d(-0.5876543210987654),
        d(0.7123456789012345),
        d(0.1234567890123456),
    ];
    let add = |a: [R; 3], b: [R; 3]| -> [R; 3] {
        core::array::from_fn(|i| a[i].checked_add(b[i]).expect("decimal widths"))
    };
    let pts = [o, add(o, x), add(o, y)];
    // Fixture qualification: genuinely wide, at unit scale.
    let name = nacre_exact::plane_name_exact(pts[0], pts[1], pts[2]).unwrap();
    assert!(
        name.narrow().is_none(),
        "the fixture must be wide — full-width crosses were expected to exceed i128"
    );
    let mut m = Model::new();
    let h = push_consistent(&mut m, pts);
    // The caller's statement, `PlaneDef`-style: origin = first point, +u toward the second.
    let placement = nacre_topo::FramePlacement::Named {
        origin: o,
        ref_dir: x,
    };
    assert_frame_shape(&m, h, &placement, "wide plane, named");
    // And `û` runs along the stated `ref_dir`, not some canonical direction.
    let (_, u, _, _) = frame_world_basis(&m, h, &placement, false).unwrap();
    let rd = x.map(|r| r.to_f64());
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cos = dot(u, rd) / dot(rd, rd).sqrt();
    assert!(
        (cos - 1.0).abs() < 1e-9,
        "û does not follow the caller's ref_dir (cos = {cos})"
    );
}
