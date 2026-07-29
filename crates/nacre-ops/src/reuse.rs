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
use crate::planes::{FaceInfo, PlaneGeom, SolidSide};
use crate::{BoolKind, he_start};
use nacre_cip::{Pt3, orient3d_filter};
use nacre_scalar::Orient;
use nacre_store::Handle;
use nacre_topo::{Model, Origin, Solid, Vertex};
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

/// Every vertex of a solid as an exact [`Pt3`], or `None` when one of them is `Discovered`.
///
/// **A discovered vertex is an implicit point** — the meet of three surfaces — so it has no
/// rational base to stand on. Handing [`Pt3`] the rounded coordinate as if it were the definition
/// would make the shared-motion exact path answer about a point that is not the one asked for,
/// which is a silent wrong answer rather than a slow one. So this declines instead, and the class
/// falls back to [`ClassPlan::Arrange`].
///
/// This costs nothing where it matters: in a fold it is the *incoming* operand that is queried,
/// and that one is freshly built.
fn solid_points(model: &Model, s: Handle<Solid>) -> Option<Vec<Pt3>> {
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
                    let v = model.vertices.get(vh);
                    out.push(match v.origin {
                        Origin::Constructed => Pt3::exact(v.point.as_array())?,
                        Origin::Discovered { .. } => return None,
                        Origin::Moved { base, motion } => {
                            let chain = crate::rotated_vertex::motion_chain(model, motion);
                            let bp = model.vertices.get(base).point.as_array();
                            let base = crate::rotated_vertex::coord_rat(bp).ok()?;
                            crate::rotated_vertex::replay(Pt3::at(base), &chain)
                        }
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
fn plane_misses(geom: &PlaneGeom, q: &[Pt3]) -> bool {
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
    geom: &[PlaneGeom],
    class_owner: &[Option<SolidSide>],
) -> Vec<ClassPlan> {
    let mut plans = vec![ClassPlan::Arrange; geom.len()];
    if reuse == ClassReuse::Off || !class_owner.iter().any(|o| o.is_some()) {
        return plans;
    }
    // Built lazily and once: the query set is the *other* operand's, so each is needed only if
    // some class belongs to the one it is not.
    let mut q: [Option<Option<Vec<Pt3>>>; 2] = [None, None];
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
pub(crate) struct VertexClasses(HashMap<Handle<Vertex>, Vec<usize>>);

impl VertexClasses {
    /// Over the faces of one operand — `faces[range]` — collect each vertex's plane classes.
    pub(crate) fn of(
        model: &Model,
        faces: &[FaceInfo],
        plane_ix: &[usize],
        range: std::ops::Range<usize>,
    ) -> VertexClasses {
        let mut map: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
        for fi in range {
            let wc = plane_ix[fi];
            let f = model.faces.get(faces[fi].face);
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for &he in &lp.half_edges {
                    let e = map.entry(he_start(model, he)).or_default();
                    if !e.contains(&wc) {
                        e.push(wc);
                    }
                }
            }
        }
        VertexClasses(map)
    }

    /// The vertex's triple, or `None` when it is not three planes.
    ///
    /// **Four or more is not a failure, it is a different question.** That is a concurrency, which
    /// the arrangement resolves through its alias table by picking one triple as the name — a
    /// choice this cannot reproduce from the solid alone. So the class is arranged instead, which
    /// is always right and merely slower.
    fn triple(&self, v: Handle<Vertex>) -> Option<[usize; 3]> {
        let c = self.0.get(&v)?;
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
#[cfg(debug_assertions)]
pub(crate) type CanonFace = (usize, bool, Vec<[usize; 3]>, Vec<Vec<[usize; 3]>>);

#[cfg(debug_assertions)]
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
    geom: &PlaneGeom,
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
        let ring = |lp: &nacre_topo::Loop| -> Option<Vec<Node>> {
            let mut r: Vec<Node> = lp
                .half_edges
                .iter()
                .map(|&he| Some(Node::Seam(canon(vc.triple(he_start(model, he))?))))
                .collect::<Option<_>>()?;
            // Rings are emitted CCW about the *class*'s outward normal (`emit_faces`), and `flip`
            // alone carries which chamber is material. A face wound against the class frame is
            // therefore reversed here rather than signalled downstream.
            if !same {
                r.reverse();
            }
            Some(r)
        };
        let f = model.faces.get(fa.face);
        out.push(LocalFace {
            plane_idx: wc,
            loop_nodes: ring(&f.outer)?,
            inner: f.inner.iter().map(ring).collect::<Option<_>>()?,
            flip: !same,
        });
    }
    Some(out)
}
