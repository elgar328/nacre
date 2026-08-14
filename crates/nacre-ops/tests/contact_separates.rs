//! Two bodies that only touch come back as two bodies — and one that cannot part does not.
//!
//! A `Fuse` whose operands meet along a line or at a point has an answer: the pieces themselves.
//! The kernel already returns several solids (disjoint operands do, `knife_edge.rs`), Cut and
//! Common already accept these contacts, and OCCT answers `SOLID 2` for all of them. What it must
//! *not* do is separate a body from itself: where the material loops around the contact, cutting
//! there leaves one piece and no 2-manifold contains it, so the reject stays.
//!
//! The criterion is connectivity through **manifold contacts only** — an edge two faces use. It is
//! the 3D reading of the rule `unify_coplanar_faces` already applies to coplanar faces one
//! dimension down ("faces that merely lie on the same plane without touching must each survive on
//! their own").

use nacre_math::{Point2, Point3};
use nacre_ops::{
    BoolError, BoolKind, DatumDef, OpOutput, Operation, Profile2d, RejectReason, SketchPlane,
    apply, boolean,
};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, VertexDef};

fn cube(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
    let s = m.add_cuboid(Point3::from_array(lo), Point3::from_array(hi));
    m.rebuild_adjacency();
    s
}

fn volume(m: &Model, s: Handle<Solid>) -> f64 {
    nacre_props::mass_props(m, s).expect("props").volume
}

/// Every vertex handle a solid reaches, through its outer shell and its cavities.
fn vertices_of(m: &Model, s: Handle<Solid>) -> Vec<Handle<nacre_topo::Vertex>> {
    let solid = m.solids.get(s);
    let mut out = Vec::new();
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            let f = m.faces.get(fh);
            for l in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &l.half_edges {
                    out.extend(m.edges.get(he.edge).vertices);
                }
            }
        }
    }
    out.sort_unstable_by_key(|h| h.index());
    out.dedup();
    out
}

/// The surfaces a solid actually has faces on — what a vertex of it is allowed to be named by.
fn foreign_definitions(m: &Model, s: Handle<Solid>) -> Vec<String> {
    let solid = m.solids.get(s);
    let mut own = Vec::new();
    for &sh in std::iter::once(&solid.outer).chain(solid.cavities.iter()) {
        for &fh in &m.shells.get(sh).faces {
            own.push(m.faces.get(fh).surface);
        }
    }
    let mut bad = Vec::new();
    for v in vertices_of(m, s) {
        // ★ A `match`, not a `let`: `VertexDef` has a second variant (`OnSeam`), and a `let`
        // binding would need an else-arm that silently passes everything it cannot read.
        match m.vertices.get(v).def {
            VertexDef::ThreePlane(tri) => {
                for t in tri {
                    if !own.contains(&t) {
                        bad.push(format!("vertex {} names surface {}", v.index(), t.index()));
                    }
                }
            }
            VertexDef::OnSeam(pair) => {
                for t in pair {
                    if !own.contains(&t) {
                        bad.push(format!("vertex {} names surface {}", v.index(), t.index()));
                    }
                }
            }
        }
    }
    bad
}

/// ① Two cubes sharing exactly the vertical line `x = 1, y = 1`.
///
/// ★ **`validate` is the strong check here.** It counts edge uses over the *whole model*, so a
/// contact edge left welded between the two bodies would come back `NonManifoldEdge` — its silence
/// is the evidence that each body got its own copy. The disjoint handle sets say the same thing in
/// one specific proposition, so a future change cannot make both vacuous at once.
#[test]
fn two_cubes_touching_along_an_edge_come_back_as_two_bodies() {
    let mut m = Model::new();
    let a = cube(&mut m, [0.0; 3], [1.0; 3]);
    let b = cube(&mut m, [1.0, 1.0, 0.0], [2.0, 2.0, 1.0]);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a line contact separates");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2, "two bodies that only touch are two bodies");
    let total: f64 = out.iter().map(|&s| volume(&m, s)).sum();
    assert!((total - 2.0).abs() < 1e-12, "volumes sum to 2, got {total}");
    assert_eq!(
        nacre_validate::validate(&m),
        Vec::new(),
        "a welded contact edge would read as NonManifoldEdge here"
    );
    let (u, v) = (vertices_of(&m, out[0]), vertices_of(&m, out[1]));
    assert!(
        u.iter().all(|x| !v.contains(x)),
        "the two bodies share a vertex handle: {u:?} / {v:?}"
    );
}

/// ② Two cubes sharing exactly the point `(1,1,1)`.
#[test]
fn two_cubes_touching_at_a_corner_come_back_as_two_bodies() {
    let mut m = Model::new();
    let a = cube(&mut m, [0.0; 3], [1.0; 3]);
    let b = cube(&mut m, [1.0; 3], [2.0; 3]);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a point contact separates");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 2);
    let total: f64 = out.iter().map(|&s| volume(&m, s)).sum();
    assert!((total - 2.0).abs() < 1e-12, "volumes sum to 2, got {total}");
    assert_eq!(nacre_validate::validate(&m), Vec::new());
}

/// ③ ★★ **The negative control, along a line.** A and B meet only along `x = 2, y = 2`, but a
/// bridge overlaps both, so the material loops around the contact: cut it there and one piece
/// remains. There is no pair of solids to return, and the reject stays.
///
/// Without this the change is indistinguishable from deleting the reject.
#[test]
fn a_body_pinched_along_a_line_is_still_a_reject() {
    let mut m = Model::new();
    let a = cube(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let bridge = cube(&mut m, [1.0, 0.3, 0.2], [3.0, 3.0, 0.8]);
    let b = cube(&mut m, [2.0, 2.0, 0.0], [4.0, 4.0, 1.0]);
    let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
    m.rebuild_adjacency();
    assert_eq!(
        boolean(&mut m, BoolKind::Fuse, ab[0], b),
        Err(BoolError::Unsupported {
            reason: RejectReason::NonManifoldResultEdge
        })
    );
}

/// ④ ★★ **The negative control, at a point.** The same loop, closed through two bridge pieces that
/// go around the contact rather than through it, so `(2,2,1)` stays a pinch.
#[test]
fn a_body_pinched_at_a_point_is_still_a_reject() {
    let mut m = Model::new();
    let a = cube(&mut m, [0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
    let g1 = cube(&mut m, [1.0, 0.3, 0.2], [4.5, 1.3, 0.8]);
    let g2 = cube(&mut m, [3.5, 0.3, 0.2], [4.5, 3.8, 1.6]);
    let b = cube(&mut m, [2.0, 2.0, 1.0], [4.0, 4.0, 2.0]);
    let t1 = boolean(&mut m, BoolKind::Fuse, a, g1).expect("a and g1 overlap");
    m.rebuild_adjacency();
    let t2 = boolean(&mut m, BoolKind::Fuse, t1[0], g2).expect("g1 and g2 overlap");
    m.rebuild_adjacency();
    assert_eq!(
        boolean(&mut m, BoolKind::Fuse, t2[0], b),
        Err(BoolError::Unsupported {
            reason: RejectReason::NonManifoldVertex
        })
    );
}

/// A prism on the world XY plane raised from `z` by `dist`.
fn prism(m: &mut Model, pts: &[[f64; 2]], z: f64, dist: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("poly");
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(
                SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
            ),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    };
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("extrude") else {
        unreachable!()
    };
    m.rebuild_adjacency();
    solid
}

/// The block every cavity below is cut from: `[0,4]² × [0,3]`.
fn block(m: &mut Model) -> Handle<Solid> {
    prism(
        m,
        &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]],
        0.0,
        3.0,
    )
}

/// ⑤ ★★★ **The control that tells the *unit* apart.** Two square voids inside one block, meeting
/// along the vertical line `x = 2, y = 2`. Four faces use that segment — two from each void — and
/// the count is what says no 2-manifold contains it.
///
/// Handles are minted per **output solid** (a material component and the cavities nested in it),
/// not per connected component. Scope them per component and each void gets its own copies of
/// those vertices, the shared segment becomes two edges of two uses each, **the count stops seeing
/// anything**, and a solid whose material has zero thickness along that line comes back clean.
/// Nothing else in this file can catch that: ⑤b is caught by the self-touch test and ⑤c before the
/// pinch is reached, and both group a void with its host under either rule.
#[test]
fn two_cavities_meeting_along_a_line_are_still_a_reject() {
    let mut m = Model::new();
    let b = block(&mut m);
    let v1 = prism(
        &mut m,
        &[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]],
        1.0,
        1.0,
    );
    let hollow = boolean(&mut m, BoolKind::Cut, b, v1).expect("the first void is ordinary");
    m.rebuild_adjacency();
    let v2 = prism(
        &mut m,
        &[[2.0, 2.0], [3.0, 2.0], [3.0, 3.0], [2.0, 3.0]],
        1.0,
        1.0,
    );
    assert_eq!(
        boolean(&mut m, BoolKind::Cut, hollow[0], v2),
        Err(BoolError::Unsupported {
            reason: RejectReason::NonManifoldResultEdge
        }),
        "two voids pinching the material between them is not a solid"
    );
}

/// ⑤b **A void touching a wall from inside** — its corner on the face, not on an edge, so the
/// contact is a line in the *interior* of the block's wall. Every edge is still used twice and the
/// count sees nothing; the self-touch test is what names it. Kept because it is the cavity twin of
/// the self-touching solids, and it must not start passing when contacts learn to separate.
#[test]
fn a_cavity_touching_its_host_s_wall_is_still_a_reject() {
    let mut m = Model::new();
    let b = block(&mut m);
    let void = prism(&mut m, &[[0.0, 2.0], [2.0, 1.0], [2.0, 3.0]], 1.0, 1.0);
    assert_eq!(
        boolean(&mut m, BoolKind::Cut, b, void),
        Err(BoolError::Unsupported {
            reason: RejectReason::SelfTouchingResult
        })
    );
}

/// ⑤c ★ **Every candidate node grazes.** A diamond void whose four corners sit on the four walls:
/// a component's nesting depth is decided by casting from one of its own nodes, and here every one
/// of them lies on the block's boundary. The kernel says that rather than guessing.
///
/// A new population, not a new defect — until contacts separated, a touching void and its host
/// were one component and the question was never asked. Locked because it is the one reject this
/// change adds that a caller can actually meet, and its name should not drift.
#[test]
fn a_void_whose_every_corner_grazes_declines_by_name() {
    let mut m = Model::new();
    let b = block(&mut m);
    let diamond = prism(
        &mut m,
        &[[0.0, 2.0], [2.0, 0.0], [4.0, 2.0], [2.0, 4.0]],
        1.0,
        1.0,
    );
    assert_eq!(
        boolean(&mut m, BoolKind::Cut, b, diamond),
        Err(BoolError::Unsupported {
            reason: RejectReason::NoClearRay
        })
    );
}

/// ⑥ ★ Each body names itself by **its own** planes. The definition triple is derived per solid, so
/// a contact vertex of one body can never be named by a plane only the other body has a face on —
/// the defect `reconstruct`'s `debug_assert` catches in debug builds and this catches in release.
#[test]
fn each_body_names_itself_by_its_own_planes() {
    let mut m = Model::new();
    let a = cube(&mut m, [0.0; 3], [1.0; 3]);
    let b = cube(&mut m, [1.0, 1.0, 0.0], [2.0, 2.0, 1.0]);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("a line contact separates");
    m.rebuild_adjacency();
    for &s in &out {
        assert_eq!(
            foreign_definitions(&m, s),
            Vec::<String>::new(),
            "a body named by a plane it has no face on"
        );
    }
}

/// ⑦ Cut and Common were already accepting these contacts and must not move.
#[test]
fn cut_and_common_are_unchanged_by_a_contact() {
    for (lo, hi) in [
        ([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]), // line
        ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]), // point
    ] {
        let mut m = Model::new();
        let a = cube(&mut m, [0.0; 3], [1.0; 3]);
        let b = cube(&mut m, lo, hi);
        let cut = boolean(&mut m, BoolKind::Cut, a, b).expect("cut");
        assert_eq!(cut.len(), 1, "a contact removes nothing");
        m.rebuild_adjacency();
        assert!((volume(&m, cut[0]) - 1.0).abs() < 1e-12);

        let mut m = Model::new();
        let a = cube(&mut m, [0.0; 3], [1.0; 3]);
        let b = cube(&mut m, lo, hi);
        assert!(
            boolean(&mut m, BoolKind::Common, a, b)
                .expect("common")
                .is_empty(),
            "a contact has no volume in common"
        );
    }
}
