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

/// **Every moved corner is its realization, and a corner on one chain replays near it.**
///
/// The cache is the definition's realization — read out to its nearest `f64` when the corner is
/// pushed (`push_vertex_realized`) — so `realize_vertex` gives it back bit for bit, whatever chains
/// the corner's carriers hold. `tests/invariants/realize_vertex.rs`'
/// `an_operations_cache_is_the_realization` holds that for boxes, booleans and a tilted frame; this
/// table holds it for motion chains.
///
/// Where the three carriers share one chain, the corner is also solved from their pre-motion names
/// and replayed through that chain in `f64` (`replay_chain_coord`), the road the judging witnesses
/// take. That replay is not the cache — it lands within the construction epsilon of it, not on it
/// (measured on corners drawn on a frame node: 88 of 800 bit-equal, the worst 1,591 ulps off) — and
/// what it checks is that the chain is read root to leaf and that each node's `f64` step is the
/// motion the node names. A corner whose carriers hold *different* chains — a cap a first turn fixed
/// and a later turn moved — has no one chain to replay, and the reuse road gives it up too
/// (`reuse::solid_points`): counted, per chain.
#[test]
fn a_moved_corner_is_its_realization_and_replays_near_it() {
    // ★ Chains, not one turn, and none of them "nice": pivots off the origin, angles whose
    // realization is nowhere near a quadrantal one, translations and reflections between turns.
    let mv = |iso| Step::Move(Box::new(iso));
    let chains: [Vec<Step>; 5] = [
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
            // ★ A reflection between two turns, the node whose parity the chain has to carry.
            Step::Flip(Axis::X, R::new(1, 2).unwrap()),
            mv(turn(Axis::Z, [0, 0, 0], (271, 4))),
        ],
        // ★ A same-axis depth-2 chain: the caps are fixed by *both* turns (restated
        // twice, still world-stated), so their corners exercise the fixed-carrier
        // licence at depth 2 — walls [Z, Z], caps None, `chain_fixes_plane` proving
        // the caps ride the whole chain.
        vec![mv(rot30z()), mv(turn(Axis::Z, [2, -1, 3], (2749, 71)))],
        // ★ A reflection **last** — the node a definition and its coordinate are easiest to part
        // on — at an offset that is not zero, where the mirror's `2c` term would vanish. The caps
        // are fixed by the Z turn and by an X mirror alike, so every corner rides the walls' one
        // [Z, M] chain through the licence.
        vec![mv(rot30z()), Step::Flip(Axis::X, R::new(1, 3).unwrap())],
    ];
    // Per chain, `(replayed, declined)` edge-end visits — a cuboid's 48 each. Measured: the two
    // mixed-axis chains leave every corner on two chains; the others put every corner on one.
    let want: [(usize, usize); 5] = [(48, 0), (0, 48), (0, 48), (48, 0), (48, 0)];
    for (chain, want) in chains.iter().zip(want) {
        let mut m = Model::new();
        let s = crate::fixtures::cuboid(
            &mut m,
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
        let (mut replayed, mut declined) = (0, 0);
        let sh = m.solid(r).outer;
        for &fh in &m.shell(sh).faces {
            for he in &m.face(fh).outer.half_edges {
                for &vh in m.edge(he.edge).vertices.iter() {
                    // The cache is the realization, on every corner.
                    let stored = m.vertex_point(vh).as_array();
                    assert!(
                        matches!(m.vertex_cache(vh), nacre_topo::PointCache::Realized { .. }),
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
                    let nacre_topo::Vertex::ThreePlane(tri) = *m.vertex(vh) else {
                        panic!("a cuboid corner is a three-plane point");
                    };
                    let motion_of = |h| match m.surface(h) {
                        nacre_topo::Surface::Plane { motion, .. } => *motion,
                        nacre_topo::Surface::Cylinder { motion, .. } => *motion,
                    };
                    // The replayable population mirrors `solid_points`' own criterion: the moved
                    // carriers share one leaf, and a world-stated carrier among them is provably
                    // fixed by that chain (the restatement licence).
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
                    let replay = replay_chain_coord(
                        &m,
                        [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
                        leaf,
                    )
                    .unwrap();
                    for k in 0..3 {
                        assert!(
                            (replay[k] - stored[k]).abs() <= 1e-9 * (1.0 + stored[k].abs()),
                            "f64 replay within the construction epsilon: {replay:?} vs {stored:?}"
                        );
                    }
                    replayed += 1;
                }
            }
        }
        assert_eq!(
            (replayed, declined),
            want,
            "the corners this chain replays and declines"
        );
    }
}

/// A chain is read root-to-leaf, so replaying it applies the rotations in the order they
/// happened — the reverse would be a different motion whenever the axes differ.
#[test]
fn a_chain_reads_root_to_leaf() {
    let mut m = Model::new();
    let s = crate::fixtures::cuboid(
        &mut m,
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
    let (h, _) = m.push_plane(
        pl,
        nacre_topo::CacheStanding::Unrealized,
        pts,
        None,
        nacre_topo::Orientation::Forward,
    );
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

/// A frame node chain realized the way [`frame_world_basis`] realizes one — `(origin, û, v̂, ŵ)`
/// as replayed points — so a node built by hand is read exactly as the narrow road's is.
fn chain_basis(chain: &[MoveNode]) -> WorldBasis {
    let at = |p: [i128; 3]| -> [f64; 3] {
        replay(WitnessPoint::at(p.map(R::from_int)), chain)
            .expect("the frame replays")
            .coord_approx()
    };
    let o = at([0, 0, 0]);
    let axis = |p: [i128; 3]| {
        let q = at(p);
        [q[0] - o[0], q[1] - o[1], q[2] - o[2]]
    };
    (o, axis([1, 0, 0]), axis([0, 1, 0]), axis([0, 0, 1]))
}

/// ★★★ **One name, one placement, one `flip` — one frame, whichever road builds it.**
///
/// The narrow road (`frame_chain`), the wide road (`WideFrame`, forced here on a name the narrow
/// road also takes) and the judged road (`FrameThrough` over the plane's own three points —
/// `Canonical` only, it has no `Named`) are three spellings of one convention, and a sketch must
/// not move when a name widens or a datum turns nameless. So the roads are compared with each
/// other over every placement × `flip`, rather than each against a property of its own: a new
/// disagreement between any two is red here whichever cell it lands in.
///
/// ★ The anchor that keeps the three from agreeing on the wrong thing: a `Canonical` frame's `û`
/// is the arbitrary axis of **its own** `ŵ` — `ẑ × ŵ`, or `ŷ × ŵ` on a horizontal plane — in
/// either `flip`, compared with its sign. That is the convention every road can spell: the
/// judged road proves no coefficient's sign, so a `û` tied to the canonical name's sign (first
/// nonzero component positive) has no spelling there.
#[test]
fn every_road_frames_a_flipped_plane_alike() {
    let q = |n: i128, d: i128| R::new(n, d).unwrap();
    let z = q(0, 1);
    // A tilted plane (canonical name `[15, 10, 6, −30]`) and a horizontal one, which takes the
    // arbitrary axis's vertical branch. Each triple spans the canonical name's own direction, so
    // the judged road's unflipped normal and the narrow road's face the same way.
    let planes = [
        (
            "tilted",
            [[q(2, 1), z, z], [z, q(3, 1), z], [z, z, q(5, 1)]],
            [q(2, 1), q(-3, 1), z],
            [q(2, 1), q(-3, 1), q(1, 1)],
        ),
        (
            "horizontal",
            [
                [z, z, q(1, 2)],
                [q(1, 1), z, q(1, 2)],
                [z, q(1, 1), q(1, 2)],
            ],
            [q(1, 1), q(1, 1), z],
            [q(1, 1), q(1, 1), q(1, 2)],
        ),
    ];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let near = |a: &WorldBasis, b: &WorldBasis| {
        [(a.0, b.0), (a.1, b.1), (a.2, b.2), (a.3, b.3)]
            .iter()
            .all(|(x, y)| (0..3).all(|k| (x[k] - y[k]).abs() < 1e-12))
    };
    let mut wrong: Vec<String> = Vec::new();
    for (what, pts, ref_dir, off_plane) in planes {
        let mut m = Model::new();
        let h = push_consistent(&mut m, pts);
        let name = m.surface_name.get(&h).expect("a named plane").clone();
        assert!(
            name.narrow().is_some(),
            "{what}: the narrow road must take it"
        );
        for flip in [false, true] {
            let narrow = |placement: &nacre_topo::FramePlacement| {
                frame_world_basis(&m, h, placement, flip).expect("narrow")
            };
            let wide = |wf: Option<nacre_judge::WideFrame>| {
                chain_basis(&[MoveNode::FrameWide(wf.expect("the wide road opens"))])
            };

            let canonical = narrow(&nacre_topo::FramePlacement::Canonical);
            let roads = [
                (
                    "wide",
                    wide(nacre_judge::WideFrame::canonical_of(&name, flip)),
                ),
                (
                    "judged",
                    chain_basis(&[MoveNode::FrameThrough(Box::new(
                        nacre_judge::FrameThrough::of(
                            pts.map(|p| nacre_judge::JudgedPoint::Pure(WitnessPoint::at(p))),
                            flip,
                        )
                        .expect("a healthy plane frames"),
                    ))]),
                ),
            ];
            for (road, b) in roads {
                if !near(&canonical, &b) {
                    wrong.push(format!(
                        "{what} Canonical flip={flip}: narrow {canonical:?} vs {road} {b:?}"
                    ));
                }
            }
            let (_, u, _, w) = canonical;
            let horizontal = w[0].abs() < 1e-12 && w[1].abs() < 1e-12;
            let raw = if horizontal {
                [w[2], 0.0, -w[0]] // ŷ × ŵ
            } else {
                [-w[1], w[0], 0.0] // ẑ × ŵ
            };
            let cos = dot(u, raw) / dot(raw, raw).sqrt();
            if (cos - 1.0).abs() > 1e-12 {
                wrong.push(format!(
                    "{what} Canonical flip={flip}: û {u:?} is not the arbitrary axis of ŵ {w:?} \
                     (cos {cos})"
                ));
            }

            // A `Named` `ref_dir` in the plane, and one with a component along the normal —
            // `SketchFrame::named` accepts it, and its in-plane part is what names `+u`.
            for (side, rd) in [("in-plane", ref_dir), ("off-plane", off_plane)] {
                let placement = nacre_topo::FramePlacement::Named {
                    origin: pts[0],
                    ref_dir: rd,
                };
                let n = narrow(&placement);
                let b = wide(nacre_judge::WideFrame::named_of(&name, &pts[0], &rd, flip));
                if !near(&n, &b) {
                    wrong.push(format!(
                        "{what} Named {side} flip={flip}: narrow {n:?} vs wide {b:?}"
                    ));
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "roads that frame one statement differently:\n{}",
        wrong.join("\n")
    );
}
