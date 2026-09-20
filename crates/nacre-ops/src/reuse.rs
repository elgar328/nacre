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
//! [`nacre_judge::orient3d_filter`], which returns a sign only when its error bound clears zero and
//! `None` otherwise; `None` means "arrange it after all", which is what the engine did before this
//! module existed. **So a missed proof costs time and nothing else, and no tolerance enters** —
//! the same reason the kernel is allowed to have a fast path at all.

use crate::boolean::LocalFace;
use crate::combinatorics::{Canon3, NodeId};
use crate::planes::{ClassIx, FaceRow, SolidSide, WorkingPlane};
use crate::{BoolKind, he_start};
use nacre_judge::{WitnessPoint, orient3d_filter};
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
/// One road for every vertex: its base is **solved from its definition** — the three carriers'
/// narrow names, in the frame the planes are stated in ([`nacre_scalar::three_planes_rat`]) —
/// and the coordinate handed to the predicates is that base's own nearest `f64` with the
/// rounding it carries ([`WitnessPoint::at_nearest`] on the world arm; [`WitnessPoint::at`],
/// measured, on the motion arm, whose replay transports it). Two arms after that:
/// * all three carriers world-stated (`motion: None`) — the base is the world point, done;
/// * the moved carriers sharing one motion, any world-stated carrier **fixed by that
///   chain** — the base is the pre-motion point; replay the chain (measured bit-identical to
///   the stored base-and-replay road, 8/8). (A fixed plane's world equation *is* its
///   pre-motion equation, which is what admits a restated cap.)
///
/// ★ **The cache is never read here.** The world arm used to hand `vertex_point` to
/// `WitnessPoint::exact` — the rounded `f64` lifted back to a rational with tol 0 — which names
/// a different point for every coordinate `f64` cannot hold (measured: at most 316 of a census
/// corpus's 6,988 operand vertices, `0.3` among them). A boolean's own result vertices used to decline
/// wholesale on the same ground ("an implicit point has no rational base"); their definition
/// has one, so they take the road too.
///
/// ★ **Mixed frames decline — unless the odd carrier is provably fixed.** Since the
/// invariant-plane restatement, a turned block's corner is a restated world cap × two
/// chained walls; the cap's world equation is *also* its equation in the walls' pre-motion
/// frame precisely when the chain fixes the plane ([`nacre_topo::Model::chain_fixes_plane`]
/// — the consumer-side twin of the producer's own gate), so the triple solves in that
/// frame and the chain replays as before. A mixed corner the chain does **not** fix — a
/// prism's base ring under a caller-stated world plane, the frame realization being
/// irrational — still has no rational pullback and declines honestly (Arrange — slower,
/// never wrong).
fn solid_points(model: &Model, s: Handle<Solid>) -> Option<Vec<WitnessPoint>> {
    use nacre_topo::Vertex;

    let sol = model.solid(s);
    let mut seen: std::collections::HashSet<Handle<Vertex>> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
        for &fh in &model.shell(sh).faces {
            let f = model.face(fh);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for &he in &lp.half_edges {
                    let vh = he_start(model, he);
                    if !seen.insert(vh) {
                        continue;
                    }
                    let tri = match *model.vertex(vh) {
                        Vertex::ThreePlane(tri) => tri,
                        // No rational base point: OnSeam and Pierce coordinates are not
                        // rational, so reuse declines and the boolean takes the slower road.
                        Vertex::OnSeam(_) | Vertex::Pierce { .. } => return None,
                    };
                    // The base is the definition's: the three narrow names, solved in the frame
                    // the planes are stated in. A wide or missing name declines.
                    let mut coeffs = [[nacre_scalar::Rat::from_int(0); 4]; 3];
                    for (o, h) in coeffs.iter_mut().zip(tri) {
                        *o = *model.surface_name.get(&h)?.narrow()?;
                    }
                    let motions = tri.map(|h| model.plane_motion(h));
                    // One shared leaf among the moved carriers, or a decline.
                    let mut leaf = None;
                    for m in motions.iter().flatten() {
                        match leaf {
                            None => leaf = Some(*m),
                            Some(l) if l == *m => {}
                            Some(_) => return None, // two histories — no shared frame
                        }
                    }
                    let Some(leaf) = leaf else {
                        // World-stated throughout: the base is the world point itself.
                        out.push(WitnessPoint::at_nearest(nacre_scalar::three_planes_rat(
                            coeffs,
                        )?));
                        continue;
                    };
                    // A world-stated carrier among chained ones is admissible iff the
                    // chain fixes it — then its equation holds in the pre-motion frame
                    // too. Otherwise: mixed frames, decline (see the doc).
                    for (c, m) in coeffs.iter().zip(&motions) {
                        if m.is_none() && !model.chain_fixes_plane(leaf, c) {
                            return None;
                        }
                    }
                    let base = nacre_scalar::three_planes_rat(coeffs)?;
                    let chain = crate::rotated_vertex::motion_chain(model, leaf)?;
                    out.push(crate::rotated_vertex::replay(
                        WitnessPoint::at(base),
                        &chain,
                    )?);
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
///
/// ★ That is *a* name, not the only one — the arrangement also names pierce points
/// ([`crate::combinatorics::NodeId::Pierce`]). This table stays three-plane by construction: its
/// population is the vertices a **solid already has**, and a solid gains a pierce vertex only when
/// the arc split starts building them.
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
        faces: &[FaceRow],
        plane_ix: &[ClassIx],
        range: std::ops::Range<usize>,
    ) -> VertexClasses {
        let mut vertices: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
        let mut edges: HashMap<Handle<nacre_topo::Edge>, Vec<usize>> = HashMap::new();
        for fi in range {
            let wc = plane_ix[fi].plane();
            let f = model.face(faces[fi].plane().face.expect("reuse only sees real faces"));
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
    /// **Four or more is not a failure, it is a different question.** That is a concurrency.
    /// Its name *is* reproducible from the solid alone — `canonical_triple` of the incident
    /// classes, the rule the ring road and the alias table share — but this pass has no judge to
    /// ask (`class_plans` runs before one is made), and passing a class through around such a
    /// vertex also has to weld its faces to arranged ones by that name, which is unmeasured. So
    /// the class is arranged instead, which is always right and merely slower; how many classes
    /// that costs is a number to read off the corpus before deciding otherwise.
    fn triple(&self, v: Handle<Vertex>) -> Option<NodeId> {
        let c = self.vertices.get(&v)?;
        let [a, b, d] = c[..] else { return None };
        Some(NodeId::three_planes(Canon3::three([a, b, d])))
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
/// One element of a canonicalized bound: a node's name, or the class a curved bound is.
///
/// ★ **The vessel is wide rather than refusing.** The differential's job is to *separate*
/// faces the two routes could disagree about, so a pierce node needs a spelling here — declining
/// would silently stop covering every face that contains one.
/// Being an enum also means **no `usize::MAX` sentinels** for a circle and a band:
/// those exist only when the vessel is a triple, which is the same
/// defect one layer down. It separates strictly more than sentinels would — a band with
/// `lo == hi` would collide with that class's circle.
/// A band's rim in the canonical key: its plane class, or a chain (whose nodes follow).
#[cfg(any(debug_assertions, test))]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) enum CanonRim {
    Circle(usize),
    Chain,
}

#[cfg(any(debug_assertions, test))]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) enum CanonNode {
    /// A three-plane vertex, by its sorted classes.
    Three([usize; 3]),
    /// A `plane ∩ plane ∩ cylinder` vertex, by everything that names it.
    Pierce {
        planes: [usize; 2],
        cyl: usize,
        root: nacre_topo::QuadRoot,
    },
    /// A whole circle, by its cylinder class.
    Circle(usize),
    /// A band, by its two rim classes.
    Band(CanonRim, CanonRim),
}

/// `(plane, flip, outer ring, hole rings)` — every field of a [`LocalFace`] but those freedoms.
#[cfg(any(debug_assertions, test))]
impl CanonNode {
    /// A node's spelling. ★ It reads the identity directly rather than through
    /// `three_plane_name`: this caller's answer for a pierce node is not "refused" but "the key
    /// says which one it is", so the door's single answer is the wrong one here.
    fn of(n: NodeId) -> CanonNode {
        match n {
            NodeId::ThreePlane(t) => CanonNode::Three(t),
            NodeId::Pierce { planes, cyl, root } => CanonNode::Pierce { planes, cyl, root },
        }
    }
}

#[cfg(any(debug_assertions, test))]
pub(crate) type CanonFace = (usize, bool, Vec<CanonNode>, Vec<Vec<CanonNode>>);

#[cfg(any(debug_assertions, test))]
pub(crate) fn canonical(faces: &[LocalFace]) -> Vec<CanonFace> {
    let ring = |r: &[NodeId]| -> Vec<CanonNode> {
        let t: Vec<CanonNode> = r.iter().map(|&n| CanonNode::of(n)).collect();
        match t.iter().enumerate().min_by_key(|(_, v)| **v) {
            Some((i, _)) => t[i..].iter().chain(&t[..i]).copied().collect(),
            None => t,
        }
    };
    let mut out: Vec<_> = faces
        .iter()
        .map(|f| {
            // ★ A **curved** bound canonicalizes by its classes rather than by nodes, because it
            // has none. The key still has to *separate* faces the two routes could disagree
            // about, so a circle contributes its cylinder class and a band its two rims —
            // spelled as one-element triples so they share the vessel with the node lists.
            let canon_rim = |r: &crate::boolean::Rim| match r {
                crate::boolean::Rim::Circle(c) => CanonRim::Circle(*c),
                crate::boolean::Rim::Chain(_) => CanonRim::Chain,
            };
            let ring_b = |b: &crate::boolean::Bound| match b {
                crate::boolean::Bound::Ring(r) => ring(r),
                crate::boolean::Bound::Circle { cyl } => vec![CanonNode::Circle(*cyl)],
                // A chain rim's nodes follow the marker, so two bands differing only in a
                // chain still separate.
                crate::boolean::Bound::Band { lo, hi } => {
                    let mut v = vec![CanonNode::Band(canon_rim(lo), canon_rim(hi))];
                    for r in [lo, hi].into_iter().filter_map(crate::boolean::Rim::ring) {
                        v.extend(ring(r));
                    }
                    v
                }
            };
            let mut inner: Vec<Vec<CanonNode>> = f.inner.iter().map(ring_b).collect();
            inner.sort();
            // A cylinder face has no plane class; `usize::MAX − k` keeps the key total and keeps
            // the two kinds apart (a plane index can never reach it).
            let cls = match f.surf {
                crate::planes::ClassIx::Plane(c) => c,
                crate::planes::ClassIx::Cyl(k) => usize::MAX - k,
            };
            (cls, f.flip, ring_b(&f.outer), inner)
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
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    range: std::ops::Range<usize>,
    vc: &VertexClasses,
    canon: impl Fn(NodeId) -> NodeId,
) -> Option<Vec<LocalFace>> {
    let mut out = Vec::new();
    for fi in range {
        if plane_ix[fi].plane() != wc {
            continue;
        }
        let fa = &faces[fi];
        // **Which way this face faces, in the class's frame.** The two normals are parallel — the
        // faces are coplanar by the very judgement that put them in one class — so their dot is
        // ±1 and no rounding can move its sign. Anything near zero would mean the class itself is
        // wrong, so it is refused rather than rounded.
        let n = geom.tri_n_out();
        let d = fa.plane().n_out.dot(n);
        // `n_out` is a unit vector but the class's is a raw cross product, so the test has to be
        // on the **cosine**, not on the dot. Comparing the dot to a constant instead reads a small
        // witness triangle as a near-perpendicular one — which rejected 94% of the classes this
        // was built for, silently and while still being correct.
        if d.abs() < 0.5 * n.norm() {
            return None;
        }
        let same = d > 0.0;
        let ring = |lp: &nacre_topo::Loop| -> Option<crate::boolean::Ring> {
            let mut r: Vec<NodeId> = lp
                .half_edges
                .iter()
                .map(|&he| Some(canon(vc.triple(he_start(model, he))?)))
                .collect::<Option<_>>()?;
            // ★ The wall of the edge leaving vertex `i` is the plane of the face on the other side
            // of that edge — carried, not derived (see `boolean::Ring`).
            let mut walls: Vec<crate::boolean::Wall> = lp
                .half_edges
                .iter()
                .map(|he| Some(crate::boolean::Wall::Plane(vc.wall(he.edge, wc)?)))
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
        let f = model.face(fa.face().expect("reuse only sees real faces"));
        out.push(LocalFace {
            surf: crate::planes::ClassIx::Plane(wc),
            outer: crate::boolean::Bound::Ring(ring(&f.outer)?),
            inner: f
                .inner
                .iter()
                .map(|l| Some(crate::boolean::Bound::Ring(ring(l)?)))
                .collect::<Option<_>>()?,
            flip: !same,
        });
    }
    Some(out)
}

#[cfg(test)]
#[path = "tests/reuse.rs"]
mod tests;
