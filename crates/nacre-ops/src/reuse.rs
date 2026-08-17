//! **A plane the other operand cannot reach carries no news.**
//!
//! A boolean asks "what is the arrangement of all these planes?" when the question it is really
//! answering is "how does `b` change `a`?". Those differ by everything the second operand cannot
//! touch — and in an incremental fold that is nearly all of it: folding eighty fins onto a hub,
//! the eightieth boolean arranges 169 plane classes of which the incoming fin reaches 27.
//!
//! So this module answers the other 142 without arranging them. A class whose plane is **proved**
//! disjoint from the other solid contributes exactly the faces its own solid already has there,
//! and those faces are restated in the arrangement's own vocabulary — plane-class triples — so
//! nothing downstream can tell which route a face came by.
//!
//! **The proof is the whole design.** Separation is established with
//! [`nacre_cip::orient3d_filter`], which returns a sign only when its error bound clears zero and
//! `None` otherwise; `None` means "arrange it after all", which is what the engine did before this
//! module existed. **So a missed proof costs time and nothing else, and no tolerance enters** —
//! the same reason the kernel is allowed to have a fast path at all.

use crate::boolean::{LocalFace, Node};
use crate::planes::{FaceInfo, SolidSide, WorkingPlane};
use crate::{BoolKind, he_start};
use nacre_cip::{WitnessPoint, orient3d_filter};
use nacre_scalar::Orient;
use nacre_store::Handle;
use nacre_topo::{Model, Solid, Vertex};
use std::collections::HashMap;

/// Whether a boolean may take the shortcut this module exists for.
///
/// **`Off` is not a fallback, it is the reference.** Debug builds run the boolean both ways and
/// require the same answer, which is the only check that can see the one thing skipping the trace
/// pass changes: the alias table. Keeping the switch permanent means the shortcut stays falsifiable
/// after the engine around it moves on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ClassReuse {
    Off,
    Proved,
}

/// What a plane class contributes to the result, once its relation to the other operand is known.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ClassPlan {
    /// Arrange it — the ordinary path, and the answer whenever nothing was proved.
    Arrange,
    /// The other solid cannot reach this plane, so `side`'s faces here survive unchanged.
    PassThrough(SolidSide),
    /// The other solid cannot reach this plane, and the operation keeps nothing without it.
    Empty,
}

/// Every vertex of a solid as an exact [`WitnessPoint`], from its **definition** — or `None` when one
/// declines, and the class falls back to [`ClassPlan::Arrange`] (slower, never wrong).
///
/// Three branches (measured against the pre-S7 `Origin` road, bit-identical where both
/// answered — the C2 differential):
/// * measured vertex (`vertex_tol` Some) — decline, an implicit point has no rational base;
/// * all three defining surfaces world-stated (`motion: None`) — the coordinate is the
///   statement, `WitnessPoint::exact`, letter-identical to the old `Constructed` arm;
/// * the moved carriers sharing one motion, any world-stated carrier **fixed by that
///   chain** — solve the triple's **narrow names in the shared pre-motion frame**
///   ([`nacre_scalar::three_planes_rat`]) and replay the chain: the very computation
///   measured bit-identical to the stored base-and-replay road (8/8). (A fixed plane's
///   world equation *is* its pre-motion equation, which is what admits a restated cap.)
///
/// ★ **Mixed frames decline — unless the odd carrier is provably fixed.** Since the
/// invariant-plane restatement, a turned block's corner is a restated world cap × two
/// chained walls; the cap's world equation is *also* its equation in the walls' pre-motion
/// frame precisely when the chain fixes the plane ([`crate::rotated_vertex::chain_fixes_plane`]
/// — the consumer-side twin of the producer's own gate), so the triple solves in that
/// frame and the chain replays as before. A mixed corner the chain does **not** fix — a
/// prism's base ring under a caller-stated world plane, the frame realization being
/// irrational — still has no rational pullback and declines honestly (Arrange — slower,
/// never wrong); C2's differential counts that population (open item 14).
fn solid_points(model: &Model, s: Handle<Solid>) -> Option<Vec<WitnessPoint>> {
    use nacre_topo::VertexDef;

    let sol = model.solids.get(s);
    let mut seen: std::collections::HashSet<Handle<Vertex>> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
        for &fh in &model.shells.get(sh).faces {
            let f = model.faces.get(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for &he in &lp.half_edges {
                    let vh = he_start(model, he);
                    if !seen.insert(vh) {
                        continue;
                    }
                    if model.vertex_tol(vh).is_some() {
                        return None; // measured — an implicit point has no rational base
                    }
                    let tri = match model.vertices.get(vh).def {
                        VertexDef::ThreePlane(tri) => tri,
                        // No rational base point: OnSeam and Branch coordinates are not
                        // rational, so reuse declines and the boolean takes the slower road.
                        // ★ The old comment's premise ("a cylinder never reaches a boolean")
                        // expires at M6-2 — the decline stays correct then, the premise does
                        // not.
                        VertexDef::OnSeam(_) | VertexDef::Branch { .. } => return None,
                    };
                    let motions = tri.map(|h| model.plane_motion(h));
                    out.push(if motions.iter().all(Option::is_none) {
                        WitnessPoint::exact(model.vertex_point(vh).as_array())?
                    } else {
                        // One shared leaf among the moved carriers, or a decline.
                        let mut leaf = None;
                        for m in motions.iter().flatten() {
                            match leaf {
                                None => leaf = Some(*m),
                                Some(l) if l == *m => {}
                                Some(_) => return None, // two histories — no shared frame
                            }
                        }
                        let leaf = leaf.expect("not the all-world arm");
                        let mut coeffs = [[nacre_scalar::Rat::from_int(0); 4]; 3];
                        for (o, h) in coeffs.iter_mut().zip(tri) {
                            *o = *model.surface_name.get(&h)?.narrow()?;
                        }
                        // A world-stated carrier among chained ones is admissible iff the
                        // chain fixes it — then its equation holds in the pre-motion frame
                        // too. Otherwise: mixed frames, decline (see the doc).
                        for (c, m) in coeffs.iter().zip(&motions) {
                            if m.is_none()
                                && !crate::rotated_vertex::chain_fixes_plane(model, leaf, c)
                            {
                                return None;
                            }
                        }
                        let base = nacre_scalar::three_planes_rat(coeffs)?;
                        let chain = crate::rotated_vertex::motion_chain(model, leaf)?;
                        crate::rotated_vertex::replay(WitnessPoint::at(base), &chain)?
                    });
                }
            }
        }
    }
    Some(out)
}

/// Is every point of `q` **provably** strictly on one side of plane class `wc`?
///
/// A signed distance to a plane is affine, so all vertices on one side puts their whole convex
/// hull there; a planar polygon lies inside the hull of its own vertices even when it is not
/// convex, so the solid's whole boundary is on that side, and a bounded solid whose boundary a
/// plane misses is one the plane misses entirely.
///
/// **That argument needs every face to be planar.** It is: `collect_planes` rejects a non-planar
/// face before any of this runs. A curved face would bulge outside its vertices' hull, so if
/// curved surfaces ever arrive here this must test their bounds, not their corners.
fn plane_misses(geom: &WorkingPlane, q: &[WitnessPoint]) -> bool {
    let t = &geom.tri_pt3;
    let mut side: Option<Orient> = None;
    for v in q {
        // Whichever way the witness triangle is wound, "all the same sign" reads the same — so no
        // orientation convention is consulted here, and none can be got wrong.
        match orient3d_filter(&t[0], &t[1], &t[2], v) {
            Some(o) if o != Orient::Zero => match side {
                None => side = Some(o),
                Some(s) if s == o => {}
                _ => return false, // the plane has points on both sides: it cuts the solid
            },
            _ => return false, // on the plane, or not proved off it
        }
    }
    side.is_some()
}

/// The per-class plan for one boolean.
///
/// The rule table has no geometry left in it, because a plane the other solid misses is a plane
/// every point of which is *outside* that solid — so `keep` collapses:
///
/// | kind     | class is `a`'s | class is `b`'s |
/// |----------|----------------|----------------|
/// | `Fuse`   | pass through   | pass through   |
/// | `Cut`    | pass through   | empty          |
/// | `Common` | empty          | empty          |
///
/// A class both operands bound (`class_owner == None`) is a coplanar contact and is always
/// arranged: the two solids meet there by definition.
pub(crate) fn class_plans(
    model: &Model,
    reuse: ClassReuse,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    geom: &[WorkingPlane],
    class_owner: &[Option<SolidSide>],
) -> Vec<ClassPlan> {
    let mut plans = vec![ClassPlan::Arrange; geom.len()];
    if reuse == ClassReuse::Off || !class_owner.iter().any(|o| o.is_some()) {
        return plans;
    }
    // Built lazily and once: the query set is the *other* operand's, so each is needed only if
    // some class belongs to the one it is not.
    let mut q: [Option<Option<Vec<WitnessPoint>>>; 2] = [None, None];
    for (wc, plan) in plans.iter_mut().enumerate() {
        let Some(side) = class_owner[wc] else {
            continue;
        };
        let other = match side {
            SolidSide::A => 1,
            SolidSide::B => 0,
        };
        let pts =
            q[other].get_or_insert_with(|| solid_points(model, if other == 0 { a } else { b }));
        let Some(pts) = pts.as_deref() else { continue };
        if !plane_misses(&geom[wc], pts) {
            continue;
        }
        *plan = match (kind, side) {
            (BoolKind::Fuse, s) => ClassPlan::PassThrough(s),
            (BoolKind::Cut, SolidSide::A) => ClassPlan::PassThrough(SolidSide::A),
            (BoolKind::Cut, SolidSide::B) | (BoolKind::Common, _) => ClassPlan::Empty,
        };
    }
    plans
}

/// Every plane class incident to each vertex of one operand — the name an arrangement vertex
/// carries, read off the solid that already has it.
///
/// The arrangement names a vertex `sorted3([w, f, third])`: the three plane classes that meet
/// there. A vertex of a solid meets exactly the classes of the faces around it, so the same name
/// is available without arranging anything.
pub(crate) struct VertexClasses {
    vertices: HashMap<Handle<Vertex>, Vec<usize>>,
    /// Each edge's two plane classes — **the wall a ring edge rides**, from the incidence rather
    /// than from the endpoint names (see `boolean::Ring`). Built in the same walk.
    edges: HashMap<Handle<nacre_topo::Edge>, Vec<usize>>,
}

impl VertexClasses {
    /// Over the faces of one operand — `faces[range]` — collect each vertex's plane classes.
    pub(crate) fn of(
        model: &Model,
        faces: &[FaceInfo],
        plane_ix: &[usize],
        range: std::ops::Range<usize>,
    ) -> VertexClasses {
        let mut vertices: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
        let mut edges: HashMap<Handle<nacre_topo::Edge>, Vec<usize>> = HashMap::new();
        for fi in range {
            let wc = plane_ix[fi];
            let f = model
                .faces
                .get(faces[fi].face.expect("reuse only sees real faces"));
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for &he in &lp.half_edges {
                    let e = vertices.entry(he_start(model, he)).or_default();
                    if !e.contains(&wc) {
                        e.push(wc);
                    }
                    let e = edges.entry(he.edge).or_default();
                    if !e.contains(&wc) {
                        e.push(wc);
                    }
                }
            }
        }
        VertexClasses { vertices, edges }
    }

    /// The class on the other side of `e` from `wc` — the wall a ring edge on `wc` rides.
    fn wall(&self, e: Handle<nacre_topo::Edge>, wc: usize) -> Option<usize> {
        let cs = self.edges.get(&e)?;
        let [a, b] = cs[..] else { return None };
        Some(if a == wc { b } else { a })
    }

    /// The vertex's triple, or `None` when it is not three planes.
    ///
    /// **Four or more is not a failure, it is a different question.** That is a concurrency, which
    /// the arrangement resolves through its alias table by picking one triple as the name — a
    /// choice this cannot reproduce from the solid alone. So the class is arranged instead, which
    /// is always right and merely slower.
    fn triple(&self, v: Handle<Vertex>) -> Option<[usize; 3]> {
        let c = self.vertices.get(&v)?;
        let [a, b, d] = c[..] else { return None };
        let mut t = [a, b, d];
        t.sort_unstable();
        Some(t)
    }
}

/// A face set as something two routes to it can be compared by.
///
/// **Neither the ring's starting vertex nor the order of the faces is part of the answer.** The
/// arrangement's rings start wherever its DCEL traversal did and a solid's start where the loop
/// was stored, so comparing them literally would report a difference that is not one — and the
/// first instinct on seeing that would be to "fix" an engine that is right. So each ring is
/// rotated to begin at its least node and the faces are sorted.
///
/// `(plane, flip, outer ring, hole rings)` — every field of a [`LocalFace`] but those freedoms.
#[cfg(any(debug_assertions, test))]
pub(crate) type CanonFace = (usize, bool, Vec<[usize; 3]>, Vec<Vec<[usize; 3]>>);

#[cfg(any(debug_assertions, test))]
pub(crate) fn canonical(faces: &[LocalFace]) -> Vec<CanonFace> {
    let ring = |r: &[Node]| -> Vec<[usize; 3]> {
        let t: Vec<[usize; 3]> = r.iter().map(|Node::Seam(t)| *t).collect();
        match t.iter().enumerate().min_by_key(|(_, v)| **v) {
            Some((i, _)) => t[i..].iter().chain(&t[..i]).copied().collect(),
            None => t,
        }
    };
    let mut out: Vec<_> = faces
        .iter()
        .map(|f| {
            let mut inner: Vec<Vec<[usize; 3]>> = f.inner.iter().map(|r| ring(r)).collect();
            inner.sort();
            (f.plane_idx, f.flip, ring(&f.loop_nodes), inner)
        })
        .collect();
    out.sort();
    out
}

/// The faces `side` already has on class `wc`, restated as arrangement output.
///
/// `None` when any vertex cannot be named — the class is then arranged, which is what the engine
/// did before and is never wrong.
#[allow(clippy::too_many_arguments)]
pub(crate) fn pass_through(
    model: &Model,
    wc: usize,
    geom: &WorkingPlane,
    faces: &[FaceInfo],
    plane_ix: &[usize],
    range: std::ops::Range<usize>,
    vc: &VertexClasses,
    canon: impl Fn([usize; 3]) -> [usize; 3],
) -> Option<Vec<LocalFace>> {
    let mut out = Vec::new();
    for fi in range {
        if plane_ix[fi] != wc {
            continue;
        }
        let fa = &faces[fi];
        // **Which way this face faces, in the class's frame.** The two normals are parallel — the
        // faces are coplanar by the very judgement that put them in one class — so their dot is
        // ±1 and no rounding can move its sign. Anything near zero would mean the class itself is
        // wrong, so it is refused rather than rounded.
        let n = geom.tri_n_out();
        let d = fa.n_out.dot(n);
        // `n_out` is a unit vector but the class's is a raw cross product, so the test has to be
        // on the **cosine**, not on the dot. Comparing the dot to a constant instead reads a small
        // witness triangle as a near-perpendicular one — which rejected 94% of the classes this
        // was built for, silently and while still being correct.
        if d.abs() < 0.5 * n.norm() {
            return None;
        }
        let same = d > 0.0;
        let ring = |lp: &nacre_topo::Loop| -> Option<crate::boolean::Ring> {
            let mut r: Vec<Node> = lp
                .half_edges
                .iter()
                .map(|&he| Some(Node::Seam(canon(vc.triple(he_start(model, he))?))))
                .collect::<Option<_>>()?;
            // ★ The wall of the edge leaving vertex `i` is the plane of the face on the other side
            // of that edge — carried, not derived (see `boolean::Ring`).
            let mut walls: Vec<usize> = lp
                .half_edges
                .iter()
                .map(|he| vc.wall(he.edge, wc))
                .collect::<Option<_>>()?;
            // Rings are emitted CCW about the *class*'s outward normal (`emit_faces`), and `flip`
            // alone carries which chamber is material. A face wound against the class frame is
            // therefore reversed here rather than signalled downstream.
            if !same {
                r.reverse();
                // Edge `i` leaves node `i`; reversing the nodes re-pairs them, so the walls follow
                // the same permutation shifted by one.
                walls.reverse();
                walls.rotate_left(1);
            }
            Some(crate::boolean::Ring::new(r, walls))
        };
        let f = model
            .faces
            .get(fa.face.expect("reuse only sees real faces"));
        out.push(LocalFace {
            plane_idx: wc,
            loop_nodes: ring(&f.outer)?,
            inner: f.inner.iter().map(ring).collect::<Option<_>>()?,
            flip: !same,
        });
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
    /// when the plane is not one the model already holds (a seed, or a face's).
    fn datum_frame(m: &mut Model, plane: crate::SketchPlane) -> crate::SketchFrame {
        match crate::apply(
            m,
            &crate::Operation::DatumPlane {
                def: crate::DatumDef::Stated(plane),
            },
        ) {
            Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
            other => panic!("stating a plane: {other:?}"),
        }
    }

    use crate::{OpOutput, Operation, apply};
    use nacre_math::{Point2, Point3, Vector3};
    use nacre_scalar::{Angle, Axis, Isometry, Rat, Rotation};

    fn square(a: f64, b: f64) -> crate::Profile2d {
        crate::Profile2d::polygon(vec![
            Point2::from_array([a, a]),
            Point2::from_array([b, a]),
            Point2::from_array([b, b]),
            Point2::from_array([a, b]),
        ])
        .unwrap()
    }

    /// The def road's answer on one solid: `Some` with every point realized finite, or `None`.
    fn assert_answers(m: &Model, s: Handle<Solid>, want_some: bool, what: &str) {
        match solid_points(m, s) {
            Some(pts) => {
                assert!(
                    want_some,
                    "{what}: expected a decline, got {} points",
                    pts.len()
                );
                assert!(!pts.is_empty(), "{what}: answered with no points");
                for p in &pts {
                    assert!(
                        p.coord.iter().all(|c| c.is_finite()),
                        "{what}: a realized coordinate is not finite"
                    );
                }
            }
            None => assert!(!want_some, "{what}: expected an answer, got a decline"),
        }
    }

    /// ★★★ **Open item 0's lock on the population that reaches this module.** A decimal-framed
    /// prism's constructed corners solve from their carriers' in-frame names, whose rational
    /// Cramer overflowed on every one (measured 8/8 — while all eight points fit `Rat`), and
    /// [`solid_points`] gives the whole solid up on the first failure — so a framed operand
    /// used to cost a boolean its entire class reuse. The investigation probe, promoted.
    #[test]
    fn a_framed_prisms_corners_solve_for_reuse() {
        let mut m = Model::new();
        let plane = crate::SketchPlane::from_axes(
            Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
            Vector3::from_array([0.6, 0.8, 0.0]),
            Vector3::from_array([-0.48, 0.36, 0.8]),
        );
        let frame = datum_frame(&mut m, plane);
        let OpOutput::Extrude { solid: prism, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: crate::Profile2d::polygon(vec![
                    Point2::from_array([0.1111111111111111, 0.1234567890123456]),
                    Point2::from_array([4.123456789012345, 0.2345678901234567]),
                    Point2::from_array([3.9876543210987654, 3.1234567890123459]),
                    Point2::from_array([0.2222222222222222, 2.765432109876543]),
                ])
                .unwrap(),
                dist: 2.5,
            },
        )
        .expect("the framed prism") else {
            unreachable!()
        };
        m.rebuild_adjacency();

        // ① The solve: all eight corners — the narrow Cramer alone answered none of them.
        let pts = solid_points(&m, prism).expect("a framed prism's corners solve");
        assert_eq!(pts.len(), 8, "eight corners, none given up");

        // ② One point, one realization road: every replayed coordinate IS a stored vertex
        // coordinate, bit for bit. A ulp here would mean construction and replay realize one
        // rational through two roads — a finding, not a tolerance.
        let mut stored: std::collections::HashSet<[u64; 3]> = std::collections::HashSet::new();
        let sol = m.solids.get(prism);
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shells.get(sh).faces {
                for &he in &m.faces.get(fh).outer.half_edges {
                    stored.insert(
                        m.vertex_point(he_start(&m, he))
                            .as_array()
                            .map(f64::to_bits),
                    );
                }
            }
        }
        assert_eq!(stored.len(), 8, "a prism has eight distinct corners");
        for p in &pts {
            assert!(
                stored.contains(&p.coord.map(f64::to_bits)),
                "a replayed corner {:?} is not any stored coordinate",
                p.coord
            );
        }

        // ③ The gate bites, in the consultation direction the ceiling used to close: a class
        // owned by the world cuboid asks for the *prism's* points, so it can now leave
        // `Arrange`. (The other direction — prism-owned classes consulting the cuboid — was
        // always open, so it proves nothing here.)
        let cub = m.add_cuboid(
            Point3::from_array([-2.0, -2.0, -2.0]),
            Point3::from_array([10.0, 10.0, 10.0]),
        );
        m.rebuild_adjacency();
        let setup = crate::planes::plane_index_setup(&m, prism, cub).expect("plane setup");
        let plans = class_plans(
            &m,
            ClassReuse::Proved,
            BoolKind::Fuse,
            prism,
            cub,
            &setup.geom,
            &setup.class_owner,
        );
        let opened = plans
            .iter()
            .zip(&setup.class_owner)
            .filter(|(p, o)| **o == Some(SolidSide::B) && **p != ClassPlan::Arrange)
            .count();
        assert!(
            opened > 0,
            "no cuboid-owned class left Arrange — the reopened road was not exercised"
        );
    }

    /// ★★ S7: **which populations the def road answers for**, pinned per producer. Measured
    /// against the pre-S7 `Origin` road while both existed (C2): constructed and moved agree
    /// **bit for bit**; discovered declines on both. The one recorded difference is the
    /// mixed-frame population (④): the `Origin` road answered it through the sketch-frame base
    /// vertex S7 dissolves, and the def road declines honestly — reuse falls back to Arrange,
    /// which is slower and never wrong.
    #[test]
    fn the_def_road_answers_for_the_populations_it_can_name() {
        let mut m = Model::new();
        // ① Constructed, decimal-friendly and decimal-unfriendly corners.
        let a = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([1.0, 2.0, 3.0]),
        );
        let b = m.add_cuboid(
            Point3::from_array([0.1, 5.0, 0.3]),
            Point3::from_array([1.7, 6.9, 2.2]),
        );
        m.rebuild_adjacency();
        assert_answers(&m, a, true, "constructed (integer corners)");
        assert_answers(&m, b, true, "constructed (decimal corners)");

        // ② Moved: an inexact rotation records a chain — the 8/8 population, now suite-wide.
        let turned = {
            let OpOutput::Transform { solid } = apply(
                &mut m,
                &Operation::Transform {
                    solid: b,
                    isometry: Isometry::rotation(Rotation {
                        axis: Axis::Z,
                        point: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                    }),
                },
            )
            .unwrap() else {
                unreachable!()
            };
            solid
        };
        m.rebuild_adjacency();
        assert_answers(&m, turned, true, "moved (rotated cuboid)");

        // ③ Discovered: a fuse's result declines on both roads.
        let c = m.add_cuboid(
            Point3::from_array([0.5, 0.5, 0.5]),
            Point3::from_array([1.5, 2.5, 3.5]),
        );
        m.rebuild_adjacency();
        let fused = {
            let OpOutput::Boolean { solids } = apply(
                &mut m,
                &Operation::Boolean {
                    kind: crate::BoolKind::Fuse,
                    a,
                    b: c,
                },
            )
            .unwrap() else {
                unreachable!()
            };
            solids[0]
        };
        m.rebuild_adjacency();
        assert_answers(&m, fused, false, "discovered (a fused result)");

        // ④ The recorded difference, pinned: a tilted-frame prism's base ring sits under a
        //    world-stated cap (mixed frames), and no rational pullback exists — so the def road
        //    declines where the `Origin` road answered through the base vertex S7 dissolved.
        let tilted = crate::SketchPlane::from_origin_normal(
            Point3::from_array([0.25, -0.5, 1.5]),
            Vector3::from_array([0.3141592653589793, -0.2718281828459045, 1.0]),
        )
        .unwrap();
        let prism = {
            let __frame0 = datum_frame(&mut m, tilted);
            let OpOutput::Extrude { solid, .. } = apply(
                &mut m,
                &Operation::Extrude {
                    frame: __frame0,
                    profile: square(0.5, 2.5),
                    dist: 1.1,
                },
            )
            .unwrap() else {
                unreachable!()
            };
            solid
        };
        m.rebuild_adjacency();
        assert_answers(
            &m,
            prism,
            false,
            "mixed frames (a tilted prism's base ring)",
        );
    }
}
