//! Combinatorial queries over the per-face planar arrangement (design §8 M5, sub-unit 3).
//!
//! This module answers combinatorial questions about how one solid's boundary cuts
//! a face of the other. It lives in `nacre-ops` and not in `nacre-geom` because it
//! needs `Model`/`Face`/`Edge`, and geom sits below topo (design §1: dependencies
//! flow upward only).
//!
//! Everything decided here is decided by an exact predicate. Coordinates that
//! appear (`three_planes`' cache) are never the basis of a decision — the truth of
//! a seam point is its plane triple, as it is for a measured vertex (design §4).
//!
//! # One `usize`, two meanings — now two tables
//!
//! An arrangement reasons **per plane**, but a solid gives you **faces**: an earlier boolean can
//! split one geometric plane between two faces (a base's exposed top and the cantilever underside
//! above it) whose outward normals **oppose**. Both were rows of one `planes` table, so the same
//! `usize` meant "face" here and "plane" there, and when the two meanings met in one comparison the
//! answer was wrong *silently* — four times on this branch, most recently as 117748 predicate calls
//! that read "different plane" for two faces of one plane and answered from rounding noise.
//!
//! That used to be held by a naming convention (`fp` / `fc`) and a debug-time net. It is now the
//! type: the predicates here take [`crate::planes::WorkingPlane`], which has no face geometry to offer, and
//! `plane_ix` is the one place a face index becomes a plane index (in [`loop_triples`]).
//!
//! **Exception:** the code that *defines* the classes (`crate::fill_classes` →
//! `crate::shares_or_coplanar` → `Judge::planes_coplanar`) necessarily runs before a
//! plane table exists, so it takes face indices — hence that predicate's generic `Witness` bound.

use crate::planes::{ClassIx, WorkingPlane, edge_incidence};
use crate::tolerant::Judge;
use crate::{BoolError, RejectReason, reject};
use nacre_store::Handle;
use nacre_topo::{Edge, Face, Loop, Model, Solid, Vertex};
use std::collections::HashMap;

/// A solid's edges, each with its endpoints and **the indices of its two faces**.
///
/// Named for what it holds: the pair is `planes`-table slots, i.e. *faces*, not plane classes —
/// `edge_faces` builds it from `surf_ix: HashMap<Handle<Face>, usize>`. It was `EdgePlanes`, and
/// that name is how a face index gets read as a plane one. The face→plane step is `plane_ix`, and
/// it happens in `loop_triples`, nowhere else.
pub(crate) type EdgeFaces = HashMap<Handle<Edge>, ([Handle<Vertex>; 2], [usize; 2])>;

/// Index a solid's edges by handle — [`edge_incidence`] keyed for lookup.
///
/// Hole-ring edges are in here too: `edge_incidence` walks every loop of every
/// face, and rejects an edge whose incidence is not a pair.
pub(crate) fn edge_faces(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<EdgeFaces, BoolError> {
    let mut out = EdgeFaces::new();
    for (eh, bounds, inc) in edge_incidence(model, solid, surf_ix)? {
        out.insert(eh, (bounds, inc));
    }
    Ok(out)
}

/// Order two crossings of the line `P ∩ Q` along that line: `-1` if `V_i` precedes
/// `V_j`, `+1` if it follows, `0` if they coincide.
///
/// With `V_k = P ∩ Q ∩ R_k` and `d = n_P × n_Q`, we want `sign((V_i − V_j)·d)`. Since
/// `V_j ∈ R_j`, [`three_plane_orient3d`](nacre_geom::intersect::three_plane_orient3d) gives `sign((V_i − V_j)·N_j)` for `N_j` the
/// right-hand normal of `R_j.tri`; multiplying by `sign(d·N_j)` recovers the order.
/// Both factors are exact predicates, so the comparator is a true total order.
///
/// The pair `(P, Q)` is a parameter, not the seam pair: sub-unit 3d orders two seam
/// crossings along an *edge* of `f` by calling this with `(P, R)`, the edge's own
/// two planes. No new predicate is needed for that.
pub(crate) fn order_along(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    i: usize,
    j: usize,
) -> i8 {
    jd.orient3d(p, q, i, j) * dir_sign(jd, p, q, j)
}

/// **The order of two points along `L = P ∩ Q`, whatever pins them** — the one place that rule
/// lives, and both its consumers (the tracer's node sort and [`edge_dir`]) call it.
///
/// ★★★★★ **Which direction is `order_along`'s, derived rather than matched.** [`Judge::orient3d`]
/// is *not* a four-plane determinant — it is the point `V = {p, q, i}` against the **triangle** `j`,
/// so it is symmetric in its first three. With `Vᵢ − Vⱼ = t·(n_p × n_q)` (both lie on `L`),
///
/// ```text
/// order_along = sign((Vᵢ−Vⱼ)·N_out(j)) · sign((n_p×n_q)·N_out(j)) = sign(t·(…)²) = sign(t)
/// ```
///
/// — so it orders along **`n_p × n_q` taken from the raw coefficients**, and `j`'s own orientation
/// cancels as a square. That is why the second road below asks
/// [`nacre_geom::intersect::plane_pair_dir_sign`] — the *same* primitive — for the sign of a
/// component of that direction, instead of forming a cross product of its own. (A cross product of
/// the classes' `world_rat` normals is **not** it: those may oppose the raw ones per plane, which
/// flips the direction. Measured: the derived road agrees with `order_along` 160/160 where a
/// `world_rat` cross agreed 80/160.)
///
/// **Two roads, one rule.** Plane-pinned pairs keep the integer predicates; anything a cylinder
/// pinned has no third plane to be ordered by, and goes through the `a + b√c` tower on one
/// coordinate axis. `None` is a missing description (a rational endpoint with no exact coordinates,
/// a cylinder with no world statement), never a shape this cannot order.
#[inline]
pub(crate) fn order_pinned(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    q: usize,
    a: (NodeId, EndPin),
    b: (NodeId, EndPin),
) -> Option<i8> {
    order_on(jd, cyls, p, q, PointOn::of(a), PointOn::of(b))
}

/// [`order_pinned`] for callers that already hold the points in [`PointOn`] form.
#[inline]
pub(crate) fn order_on(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    q: usize,
    a: PointOn,
    b: PointOn,
) -> Option<i8> {
    let (a, b) = (locate(jd, cyls, p, q, a)?, on_line(jd, cyls, p, q, b)?);
    order_located(jd, p, q, &a, &b)
}

/// **What a point on `L = P ∩ Q` is given as** — the input [`locate`] and [`on_line`] take.
///
/// ★★★★ **A plane pin brings no name, and that is the point.** `{p, q, pin}` *is* the name, so
/// asking a caller for one means building a `NodeId` — a sorted triple — that both constructors
/// then ignore. The arrangement's split-point sort does that twice per comparison, and it showed:
/// the sort phase ran 29% longer with the name in the signature than without it.
#[derive(Clone, Copy, Debug)]
pub(crate) enum PointOn {
    Class(usize),
    /// Must be a [`NodeId::Branch`]; the quadric road reads its cylinder and root.
    Branch(NodeId),
}

impl PointOn {
    /// ★★★★ **A cylinder pin arrives with a branch name — a producer's invariant, so this stops
    /// loudly rather than naming a refusal.** Every site that writes `EndPin::Cylinder` writes the
    /// `NodeId::Branch` beside it, so the two disagreeing is a defect in this kernel rather than a
    /// property of the model, and [`RejectReason::BranchVertexUnnamed`] would be a false sentence
    /// for it: that one says the point *is* exactly named and only the path lacks a way to carry
    /// it. The rule is `crate::planes::ClassIx::plane`'s — *"a loud panic beats a silently wrong
    /// plane."*
    /// ☑ Measured unexercised over the whole suite and the ignored sweep before it was made loud.
    #[inline]
    fn of(pt: (NodeId, EndPin)) -> PointOn {
        match pt.1 {
            EndPin::Class(r) => PointOn::Class(r),
            EndPin::Cylinder => match pt.0 {
                NodeId::Branch { .. } => PointOn::Branch(pt.0),
                n @ NodeId::ThreePlane(_) => {
                    unreachable!("a cylinder pin was written beside a three-plane name: {n:?}")
                }
            },
        }
    }
}

/// **A point on `L = P ∩ Q` with the parts a comparison needs already computed** — the form
/// [`order_located`] takes.
///
/// ★★★★★ **This exists so the rule can have one implementation *and* keep its hoists.** The
/// arrangement's interval overlay is the boolean's hottest named phase (12–24% of it, 1.4–3.2M
/// containment questions per 60-fin fold), and its speed comes from building a point's Cramer
/// parts **once per wall pair** and its `dir_sign` **once per segment endpoint**, then asking many
/// questions of them. A rule that takes bare `(name, pin)` throws both away on every call. The
/// alternative — a fast inlined copy beside the general one — is the defect shape this codebase
/// keeps finding in itself, so the hoist becomes an *argument* instead of a second spelling.
///
/// ★★ **Both sides are `Located`, and the symmetry is load-bearing.** A draft took one side bare
/// plus a separate `b_ds: Option<i8>`, which is a value that has to match arguments it cannot see
/// — the same shape as a bug this cell's review caught (passing the wrong `q`). Built here, `ds`
/// cannot disagree with the `(p, q)` it was built for. And the symmetry pays again: a segment
/// endpoint's `Located` is built once per segment, so a branch end's [`branch_meet`] is never
/// re-solved per (split point × segment).
pub(crate) enum OnLine {
    /// A point three planes name. `ds` is `dir_sign(p, q, pin)`.
    ///
    /// ★ `ds` is read when this point is the **second** argument — the formula is
    /// `at_a.orient3d(b.pin) × ds_b`, which is [`order_along`] with `a`'s half hoisted.
    /// ★ **No `name`.** A plane-pinned point on `L = P ∩ Q` *is* `{p, q, pin}`, so the name is
    /// derivable and storing it put a 32-byte `NodeId` in a `Vec` the containment loop reads
    /// millions of times. What it was for — the cold road's coordinates, and identity against the
    /// other side — needs neither: coordinates come from the derived triple, and a plane-pinned
    /// point can never *equal* a cylinder-pinned one, so cross-variant identity is `false` by type.
    Class { pin: usize, ds: i8 },
    /// A point a cylinder pins: the meet line of its two planes and the root along it.
    ///
    /// ★★★★ **Boxed, and the profile is why** — the same reason [`EdgeDir::Arc`] is: a `MeetLine`
    /// and a `QuadVal` are an order of magnitude wider than a plane pin's two words, and the
    /// arrangement stores one of these **per segment endpoint** while every question on the hot
    /// path is a plane pin. Unboxed it moved `split_at_crossings` from 23.2% of a rotated fold to
    /// 28.2%, and 12.1% of an axis-aligned one to 19.0% — with the *same trip counts*, so the cost
    /// was the width and nothing else. One allocation on the cold road buys that back.
    Branch {
        name: NodeId,
        meet: Box<(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal)>,
    },
}

/// A point in the form a comparison's **first** argument needs: an [`OnLine`] plus the Cramer
/// handle the plane road asks `orient3d` of.
///
/// ★★★★★ **The asymmetry is measured, not assumed.** A draft made both sides `Located` for the
/// symmetry, and the profile refused it: an endpoint's handle is **never read** — the plane road is
/// `at_a.orient3d(b.pin) × ds_b` — while `ImplicitPoint` is wide enough (a `OnceCell` of Cramer
/// parts) that storing two per segment moved `split_at_crossings` from 23.2% of a rotated fold to
/// 27.2%, and 12.1% of an axis-aligned one to 17.8%. First arguments are built a few times per
/// wall; second arguments are stored per segment and asked millions of times.
/// ★ What the symmetry was protecting against survives: `ds` is still built beside the `pin` it
/// belongs to, by one constructor, so it cannot disagree with the `(p, q)` it was made for.
pub(crate) enum Located<'a> {
    /// ★ **No `ds` here, and that is the point.** The plane road is
    /// `at_a.orient3d(b.pin) × ds_b` — a *first* argument's `dir_sign` is never read. Computing it
    /// anyway cost a `plane_pair_dir_sign` per wall pair and per sub-interval: measured, the cover
    /// loop went 48.6→90.8ms on an axis-aligned fold with the trip counts unchanged.
    Class {
        pin: usize,
        at: crate::tolerant::ImplicitPoint<'a, WorkingPlane>,
    },
    Branch {
        name: NodeId,
        meet: Box<(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal)>,
    },
}

impl Located<'_> {
    /// The name this point is known by, given the line it sits on — derived for a plane pin.
    pub(crate) fn name(&self, p: usize, q: usize) -> NodeId {
        match self {
            Located::Class { pin, .. } => NodeId::three_planes([p, q, *pin]),
            Located::Branch { name, .. } => *name,
        }
    }
}

impl OnLine {
    /// The name this point is known by, given the line it sits on — derived for a plane pin.
    pub(crate) fn name(&self, p: usize, q: usize) -> NodeId {
        match self {
            OnLine::Class { pin, .. } => NodeId::three_planes([p, q, *pin]),
            OnLine::Branch { name, .. } => *name,
        }
    }
}

/// Whether two located points are **the same point by identity**, without asking a predicate.
///
/// ★ Two plane-pinned points compare by **pin**, which is the integer test the arrangement's
/// overlay has always used: on one line `{p, q, pin}` is injective in `pin`, so equal pins are one
/// point. It is a *sufficient* condition and not the whole test — where four planes meet, one point
/// wears two pins — which is why callers still ask the order afterwards.
/// ★ Anything a cylinder pins has no pin payload to compare, so those compare by **name**.
#[inline]
pub(crate) fn same_point(a: &Located<'_>, b: &OnLine) -> bool {
    match (a, b) {
        (Located::Class { pin: x, .. }, OnLine::Class { pin: y, .. }) => x == y,
        (Located::Branch { name: x, .. }, OnLine::Branch { name: y, .. }) => x == y,
        // ★ A point three planes name and a point a cylinder pins are never the same *name* — the
        // two are different variants of `NodeId`. Whether they are the same **point** is a
        // different question, and the order predicate below is what asks it.
        _ => false,
    }
}

/// **Is `at` within the closed extent an edge's two ends mark out on `L = P ∩ Q`?** — endpoints
/// included.
///
/// On an endpoint it is contained, checked by **identity**, because `order_along(x, x)` is not
/// defined to return 0 (the old `strictly_inside` never compared a class with itself). Otherwise it
/// is contained iff it lies between the two — opposite order signs — or a `0` says it *is* one of
/// them by a second name.
///
/// ★ **Asked before the predicates, not after.** Identity is a *sufficient* condition for
/// containment, so answering it first skips both orientations. It is not the whole test: where four
/// planes meet, one point wears two pins, and `at` may be the group's representative while the edge
/// still remembers the other — which is what the `== 0` arms catch. Subsumed, not dropped; the
/// order between them is free.
///
/// ★★★★★ **`ends` must have been located for the very `(p, q)` handed in here.** Each [`OnLine`]
/// carries a `dir_sign` — or a `(line, s)` — made for one plane pair, and asking about a point
/// located for a *different* pair silently compares two different lines. The type does not say so:
/// `OnLine` records no pair, and giving it one widens the vector the overlay's cover loop reads
/// millions of times. **Callers hold the pair beside the ends.**
///
/// ★ **The point is a parameter, and that is what lets a caller hoist it.** The point asked about
/// is constant for every edge on one wall, so a caller locates it once per **wall pair** and every
/// edge there shares its Cramer parts. Taking a bare class instead would cap the sharing at one
/// edge's two endpoints, which is what a `_pair` predicate did.
///
/// ★ `None` is a description that could not be formed exactly, never a shape this cannot answer;
/// the callers turn it into [`crate::RejectReason::WitnessNotRational`].
#[inline]
pub(crate) fn closed_contains(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    at: &Located<'_>,
    ends: &[OnLine; 2],
) -> Option<bool> {
    let [l0, l1] = ends;
    if same_point(at, l0) || same_point(at, l1) {
        return Some(true);
    }
    let a = order_located(jd, p, q, at, l0)?;
    let b = order_located(jd, p, q, at, l1)?;
    Some(a == 0 || b == 0 || a != b)
}

/// Put a point on `L = P ∩ Q` into the form [`order_located`] takes.
///
/// ☑ **Infallible for a plane pin** — `Judge::point` and [`dir_sign`] are total — which is what
/// keeps [`order_pinned`] total on the plane/plane pair its two older callers ask about.
#[inline]
pub(crate) fn on_line(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    q: usize,
    pt: PointOn,
) -> Option<OnLine> {
    match pt {
        PointOn::Class(pin) => Some(OnLine::Class {
            pin,
            ds: dir_sign(jd, p, q, pin),
        }),
        PointOn::Branch(name) => {
            let NodeId::Branch { cyl, .. } = name else {
                return None;
            };
            let meet = branch_meet(jd, cyl, &cyls.get(cyl)?.def, name)?;
            Some(OnLine::Branch {
                name,
                meet: Box::new(meet),
            })
        }
    }
}

/// Put a point into the **first**-argument form — an [`on_line`] with the Cramer handle beside it.
///
/// ☑ **Infallible for a plane pin**, like [`on_line`], which is what keeps [`order_pinned`] total
/// on the plane/plane pair its two older callers ask about.
#[inline]
pub(crate) fn locate<'j>(
    jd: &'j Judge<'j, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    q: usize,
    pt: PointOn,
) -> Option<Located<'j>> {
    match pt {
        PointOn::Class(pin) => Some(Located::Class {
            pin,
            at: jd.point(p, q, pin),
        }),
        PointOn::Branch(name) => {
            let NodeId::Branch { cyl, .. } = name else {
                return None;
            };
            let meet = branch_meet(jd, cyl, &cyls.get(cyl)?.def, name)?;
            Some(Located::Branch {
                name,
                meet: Box::new(meet),
            })
        }
    }
}

/// **The order of two located points along `L = P ∩ Q`** — [`order_pinned`]'s body, with the
/// hoists handed in rather than rebuilt.
///
/// ★★★★★ The plane/plane test is the **first statement** on purpose: everything after it computes
/// the axis component `(k, dsign)` with three `plane_pair_dir_sign` calls, and the overlay asks
/// this millions of times on pairs that never reach there.
/// ★★★★★ **Hot and cold are separate functions, and the profile is what says so.** With the
/// quadratic road in the same body the plane road stopped being inlined, and the arrangement's two
/// containment loops — 1.4M and 3.2M questions per 60-fin fold — **doubled**: collect 78→153ms
/// rotated, cover 49→100ms axis-aligned, with the trip counts unchanged. The same shape
/// `Judge::orient3d` uses for its own cheap route.
#[inline(always)]
pub(crate) fn order_located(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    a: &Located<'_>,
    b: &OnLine,
) -> Option<i8> {
    // ☑ `order_along(jd, p, q, i, j)` exactly: `ImplicitPoint::orient3d(j)` is
    // `Judge::orient3d(p, q, pin, j)` on the same four planes with the Cramer parts cached, and
    // `ds` is that call's `dir_sign` factor.
    if let (Located::Class { at, .. }, OnLine::Class { pin: j, ds, .. }) = (a, b) {
        return Some(at.orient3d(*j) * ds);
    }
    order_located_quad(jd, p, q, a, b)
}

/// [`order_located`]'s other road: anything a cylinder pins, compared on one coordinate axis
/// through the `a + b√c` tower. Kept out of line so the plane road above stays small enough to be.
#[inline(never)]
fn order_located_quad(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    a: &Located<'_>,
    b: &OnLine,
) -> Option<i8> {
    use nacre_scalar::Orient;
    use nacre_scalar::quad::{cmp_coord_branch, cmp_coord_meet_branch};
    let axis = |k: usize| {
        let mut n = [0.0; 3];
        n[k] = 1.0;
        nacre_geom::Plane::from_point_normal(
            nacre_math::Point3::from_array([0.0; 3]),
            nacre_math::Vector3::from_array(n),
        )
    };
    let (k, dsign) = (0..3).find_map(|k| {
        let s = nacre_geom::intersect::plane_pair_dir_sign(
            &jd.planes[p].plane,
            &jd.planes[q].plane,
            &axis(k)?,
        );
        (s != 0).then_some((k, s))
    })?;
    let sign = |o: Orient| -> i8 {
        match o {
            Orient::Positive => 1,
            Orient::Negative => -1,
            Orient::Zero => 0,
        }
    };
    let rat = |n: NodeId| -> Option<nacre_scalar::MeetPoint> {
        node_coords_rat(jd, n).map(nacre_scalar::MeetPoint::Narrow)
    };
    let cmp = match (a, b) {
        // ☑ Unreachable — two plane-pinned points returned through the integer road above. Spelled
        // as a decline rather than a panic or a second rational comparison nobody would exercise.
        (Located::Class { .. }, OnLine::Class { .. }) => return None,
        (Located::Class { .. }, OnLine::Branch { meet, .. }) => sign(cmp_coord_meet_branch(
            &rat(a.name(p, q))?,
            &meet.0,
            &meet.1,
            k,
        )),
        (Located::Branch { meet, .. }, OnLine::Class { .. }) => -sign(cmp_coord_meet_branch(
            &rat(b.name(p, q))?,
            &meet.0,
            &meet.1,
            k,
        )),
        (Located::Branch { meet: m1, .. }, OnLine::Branch { meet: m2, .. }) => {
            sign(cmp_coord_branch((&m1.0, &m1.1), (&m2.0, &m2.1), k))
        }
    };
    Some(cmp * dsign)
}

/// **A vertex's identity** — the sorted plane triple that names it.
///
/// `Eq`/`Hash` give identity dedup so an A-piece and a B-piece that meet at a seam node share one
/// result vertex/edge; `Ord` gives the deterministic node order replay needs — and it is the bare
/// triple's lexicographic order, so every "smallest name wins" rule reads unchanged.
///
/// This was `boolean::Node`, spoken only by the assembler. It lives here because the arrangement
/// names the same vertices, and the variant is spelled like [`nacre_topo::VertexDef::ThreePlane`]
/// so the arrangement, the assembler and the topology store call the thing by one name.
///
/// Read it with `match`, never `let`-`else`: a new variant lights up the first and falls silently
/// into the second — a defect this repository has already had.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub(crate) enum NodeId {
    ThreePlane([usize; 3]), // sorted triple (key into the seam map)
    /// Where two plane classes' meet line crosses a cylinder's lateral surface — the point
    /// [`nacre_topo::VertexDef::Branch`] names, and the point the next rung splits a circle
    /// into arcs at.
    ///
    /// ★ **Two index spaces in one name.** `planes` are plane-class indices and `cyl` is a
    /// **cylinder**-class index (`crate::planes::ClassIx::Cyl`'s payload); they are separate
    /// numberings and a value from one is meaningless in the other. `reuse::canonical` already
    /// carries both kinds in one key, so the precedent is the file's, not this type's.
    Branch {
        /// The two cutting plane classes, ascending — the [`NodeId::ThreePlane`] precedent, and
        /// the order `root` is defined against.
        planes: [usize; 2],
        /// The cylinder class whose lateral surface the meet line crosses.
        cyl: usize,
        /// Which crossing, along the meet line of `planes` in stored order.
        root: nacre_topo::QuadRoot,
    },
}

impl NodeId {
    /// The canonical name of the point where three plane classes meet — **the only way one is
    /// made**, so "two spellings of one vertex are one name" holds by construction.
    ///
    /// ★★ **Sorting is *this variant's* canonicalization, not the definition of canonical** — see
    /// [`NodeId::branch`], whose pair carries a root that a re-sort has to restate.
    ///
    /// It does **not** check for a collapsed triple. Two names being equal is a real condition
    /// with *different answers at different callers* — `arrangement::plane_ring` declines
    /// (`CollapsedTriple`), [`loop_triples`] falls back to naming the vertex from every plane
    /// touching it — so making the constructor fallible would copy that fork to all eight minting
    /// sites.
    pub(crate) fn three_planes(mut t: [usize; 3]) -> NodeId {
        t.sort_unstable();
        NodeId::ThreePlane(t)
    }

    /// The canonical name of a `plane ∩ plane ∩ cylinder` point — **the only way one is made**, so
    /// "two spellings of one vertex are one name" holds by construction here too.
    ///
    /// ★★★ **It is not a sort.** Ordering the two classes can reverse the meet line, and the root
    /// is defined against that line, so the pair and the root move **together**. That rule is
    /// [`nacre_topo::QuadRoot::canonical`]'s and this reads it; spelling it again here is how the
    /// same point ends up with two names, one of them pointing at the other root.
    ///
    /// `first`/`second` are the two plane classes **in whatever order the caller solved them**;
    /// `root` is that solve's answer (`Double` for a tangency).
    pub(crate) fn branch(
        first: usize,
        second: usize,
        cyl: usize,
        root: nacre_topo::QuadRoot,
    ) -> NodeId {
        let (planes, root) = nacre_topo::QuadRoot::canonical([first, second], root);
        NodeId::Branch { planes, cyl, root }
    }
}

/// **The one door out of the identity and into the machinery that assumes every vertex has a
/// three-plane name** — the ray casts, the ring walks, the wall-and-handle derivations, and the
/// comparison keys. They take `[usize; 3]`, and rightly: an in-flight probe point like
/// `point_in_component`'s `{a, b, q}` is three planes without being any vertex of the arrangement,
/// so a name is the wrong type for their parameter.
///
/// **`None` is a branch point**, and the answer for every caller behind this door is the same one:
/// it has no three-plane name and none of those paths has another to give. What differs is what
/// each does about it, and that splits in two — see [`three_plane_probes`] for the split.
///
/// Sites that *dispatch* instead of declining stay out of it deliberately — [`node_coords_rat`],
/// [`branch_point`] and [`loop_winding`]'s lexicographic scan, each with its own answer for the
/// other variant. ★ And one is **not** in the ring-and-segment world at all: `reuse::canonical`
/// builds a comparison key, and its answer is that the key's vessel widens (`CanonNode`), not that
/// the question is refused.
///
/// ★★ **The gate that keeps this honest**, and the two files it exempts:
///
/// ```text
/// rg 'NodeId::(ThreePlane|Branch)' crates/ \
///   -g '!**/combinatorics.rs' -g '!**/reuse.rs' | grep -vE ':\s*//'
/// ```
///
/// It must be empty. Spelling a variant anywhere else means a site went around the door instead of
/// answering — the same gate found **fifteen** of those in the previous cell with the whole suite
/// already green, and it costs nothing and does not break on a rename.
///
/// ★★★ **Empty again since [`branch_name`] exists** (2026-08-22). Three sites in `arrangement`
/// used to open the `Branch` variant directly — `arc_split_witness`' `separates` closure and two
/// inside `split_circles` — because this door answers only the three-plane half and its twin did
/// not exist. The seam table becoming a fourth consumer is what finally paid for the twin; all
/// four go through it now.
///
/// ★ Not automated, and that is how it rotted: the three arrived a cell ago and nothing ran the
/// check, then a fourth was nearly added with the suite green. `tests/rotation_sweep.rs`'
/// `side_of` guard is the precedent for making a source scan a test.
pub(crate) fn three_plane_name(n: NodeId) -> Option<[usize; 3]> {
    match n {
        NodeId::ThreePlane(t) => Some(t),
        NodeId::Branch { .. } => None,
    }
}

/// [`three_plane_name`]'s twin — the payload of a **branch** name, `None` for a three-plane one.
///
/// The three questions its consumers ask are all payload reads: does this root *separate* (a
/// tangency's `Double` does not), which cylinder is the point on, where does it sort along the
/// meet line. Spelling the variant at those sites instead is what the gate above forbids — the
/// door is total over the enum, so a third variant becomes a compile error here rather than a
/// silent fall-through at four call sites.
pub(crate) fn branch_name(n: NodeId) -> Option<([usize; 2], usize, nacre_topo::QuadRoot)> {
    match n {
        NodeId::ThreePlane(_) => None,
        NodeId::Branch { planes, cyl, root } => Some((planes, cyl, root)),
    }
}

/// **The names of a list of *candidates*, branch points dropped.**
///
/// ★★★ **This is the licence, and its name is where the licence is stated.** Behind the door there
/// are two shapes and only one of them may lose a member:
///
/// - a **probe list** may — its consumers try each member until one decides, and an exhausted list
///   is already a named decline (`ring_in_ring` answers `NoClearRay`, `boolean`'s `first_deciding`
///   answers `Ok(None)` and leaves the rejection to its caller);
/// - a **ring** may not — dropping a node from a cyclic sign sequence produces a *different
///   polygon* and answers a different question, confidently. Those sites collect through
///   `Option<Vec<_>>` instead, which cannot drop a member even by accident.
///
/// Written as one named function rather than a `filter_map` at each site so the licence travels
/// with the call: typing this name at a ring walk is a visible category error.
pub(crate) fn three_plane_probes(nodes: impl IntoIterator<Item = NodeId>) -> Vec<[usize; 3]> {
    nodes.into_iter().filter_map(three_plane_name).collect()
}

/// **What pins an endpoint on its edge's line** — the third carrier, in the vocabulary that names
/// it.
///
/// ★ A traced segment's ends are always plane triples, so this had been a bare `usize` (the third
/// plane class) everywhere. The arc split puts a **cylinder** crossing in the middle of a segment,
/// and that point has no third *plane* — what pins it is the quadric, and its name is the
/// [`NodeId::Branch`] the edge already carries in its endpoint list. So the pin says **which kind**
/// and the name is read from beside it, rather than a second copy living here.
///
/// ★★ The two arms are two *orders*, not two spellings of one: [`order_along`] reads a class
/// through `orient3d` × `dir_sign` (integer predicates), and a branch point through the
/// `a + b√c` tower. Naming the kind is what makes the second reachable at all — a `usize` had
/// nowhere to say "not a plane".
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum EndPin {
    /// The third plane class: the point is `P ∩ wall ∩ this`.
    Class(usize),
    /// A cylinder crossing: the point is the [`NodeId::Branch`] this endpoint is named by.
    Cylinder,
}

impl EndPin {
    /// The plane class, for the sites that are still plane-only — and `None` is the honest answer
    /// where a cylinder pinned the point, never a stand-in index.
    pub(crate) fn class(self) -> Option<usize> {
        match self {
            EndPin::Class(c) => Some(c),
            EndPin::Cylinder => None,
        }
    }
}

/// **What carries a ring edge** — the plane whose meet with `P` the edge rides, or the circle an
/// arc rides.
///
/// ★ The two arms are not symmetric and should not be made so. A plane carrier is an *index*: the
/// class table answers everything about it, and the direction it gives is the same at both ends. An
/// arc carrier has to carry the cylinder itself, because the direction it gives depends on **where
/// on the circle** it is asked — which is what [`dir_at`]'s `node` is for.
#[derive(Clone, Debug)]
pub(crate) enum Carrier {
    Plane {
        /// The plane whose meet with `P` carries this edge.
        wall: usize,
        /// The travel sense along `n_P × n_wall`, when the **endpoints** cannot supply it.
        ///
        /// ★★ An arc split cuts a segment at branch points, and `order_along` speaks three-plane
        /// classes only — so a sub-segment with a cut end has nothing to derive its sense from.
        /// The split does know it (it sorted those points along the line), so it carries it here
        /// rather than leaving a hole for [`edge_dir`] to fall into. `None` on an edge whose two
        /// named ends still answer, which is every edge no split touched.
        sense: Option<i8>,
    },
    Arc(Box<ArcCarrier>),
    /// A straight edge on the **lateral surface** — see [`RulingCarrier`]. Like an arc it must
    /// carry the cylinder itself; unlike an arc its direction (`±m`) is the same at both ends.
    Ruling(Box<RulingCarrier>),
}

/// The circle an arc rides, and which way around it the arc runs.
#[derive(Clone, Debug)]
pub(crate) struct ArcCarrier {
    /// The cylinder class — the arc's **identity**, so "two arcs of one circle" is an index
    /// comparison and not a geometric one.
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    /// `true` when travel is counter-clockwise about the cylinder's axis — the sense
    /// `arrangement::split_circles` builds every `MergedArc` in, inverted for the twin half-edge.
    pub ccw: bool,
}

/// The **ruling** a straight lateral edge rides (the M6-2 rulings ladder): a wall plane parallel
/// to a cylinder's axis meets the lateral surface in up to two axis-parallel lines, and
/// `(cyl, side)` names which of the two this is.
///
/// `side` is the sign of `(x − o) · (m × n̂)` for any point `x` on the ruling, with `o`/`m` the
/// cylinder's origin/axis and `n̂` the class's **canonical** coefficients
/// ([`class_coeffs_rat`] — one spelling; the stored normal opposes the canonical one on half the
/// classes, which is the `stored_coeffs_rat` lesson).
#[derive(Clone, Debug)]
pub(crate) struct RulingCarrier {
    /// The cylinder class — the ruling's identity, with `side`.
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    /// Which of the two parallel rulings, by the sign convention above.
    pub side: i8,
    /// `true` when travel runs along `+m` — the sense the ruling split builds every
    /// `MergedRuling` in (`end[0] → end[1]` ascends the axis), inverted for the twin half-edge.
    pub up: bool,
}

impl Carrier {
    /// A plane carrier whose sense its endpoints still supply — every edge outside a split.
    pub(crate) fn plane(wall: usize) -> Carrier {
        Carrier::Plane { wall, sense: None }
    }

    /// The carrying plane class, `None` for an arc or a ruling. The named-road consumers (the ray
    /// casts, the on-ring test) speak plane classes and nothing else, so this is where they
    /// decline.
    pub(crate) fn wall(&self) -> Option<usize> {
        match self {
            Carrier::Plane { wall, .. } => Some(*wall),
            Carrier::Arc(_) | Carrier::Ruling(_) => None,
        }
    }
}

/// One edge of a ring on plane `P`, carrying **its own geometry** rather than leaving it to be
/// recovered from the two endpoint names.
///
/// ★ **Why this type exists.** A vertex name is a plane triple, and for a long time the engine read
/// an edge's supporting plane back out of its endpoints — "the class the two names share besides
/// `P`". That works only while every vertex lies on exactly three planes. It is an accident of the
/// corpus, not an invariant: let four planes meet at a point, give the point one canonical name, and
/// the shared class is **some other plane than the one the edge rides**, silently. So the walker
/// that knows the edge — the DCEL half-edge, which was told its wall — hands the geometry over
/// instead, and only rings whose provenance is *names alone* go through `ring_from_names` (test-only).
#[derive(Clone, Debug)]
pub(crate) struct RingEdge {
    /// Identity of the vertex this edge leaves.
    pub node: NodeId,
    /// Identity of the vertex it reaches.
    ///
    /// ★★ **The pins were already two and the names only one, and that asymmetry was the bug's
    /// hiding place.** A direction is a property of `(edge, node)` — for a straight edge the two
    /// ends give the same answer, so nothing ever had to say which end it meant, and a direction
    /// taken at the *wrong* node was unspellable-looking but perfectly legal. With both names here,
    /// [`dir_at`] can check that the node it is asked about is actually on this edge.
    ///
    /// ★ Measured before it was asserted: 153,798 rings in the suite, **0** where an edge's far end
    /// is not the next edge's start. "A ring is a chain" was a producer's promise until now.
    pub to: NodeId,
    /// What carries the edge — a plane's meet with `P`, or a circle ([`Carrier`]).
    pub carrier: Carrier,
    /// What pins each endpoint on the carrier — a third plane, or the cylinder an arc split put
    /// there ([`EndPin`]). Not a name: see [`RingEdge`]'s note.
    pub from_h: EndPin,
    pub to_h: EndPin,
}

/// Recover a ring's edges from its vertex names — the classic derivation, now in **one** place.
///
/// Sound exactly while each name lists all three of its planes and no more (see [`RingEdge`]).
/// ★ **Test-only since 2026-08-17**: the last production consumer (the tracer's crossed-edge
/// wall) reads the wall the producer carries (`NamedRing`) instead, so no production path
/// derives ring geometry from names any more. Kept for hand-built test rings, whose vertices
/// are clean three-plane points by construction.
/// ★ It takes bare triples, not [`NodeId`]s, because it is a **fixture constructor**: its callers
/// hold hand-written literals, so this is where those become names (`NodeId::three_planes`).
#[cfg(test)]
pub(crate) fn ring_from_names(p: usize, ring: &[[usize; 3]]) -> Result<Vec<RingEdge>, BoolError> {
    (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let shared: Vec<usize> = a.iter().copied().filter(|x| b.contains(x)).collect();
            if shared.len() != 2 || !shared.contains(&p) {
                return Err(reject(RejectReason::RingNaming));
            }
            let wall = shared[usize::from(shared[0] == p)];
            let third = |t: [usize; 3]| t.iter().copied().find(|&x| x != p && x != wall);
            let (Some(from_h), Some(to_h)) = (third(a), third(b)) else {
                return Err(reject(RejectReason::RingNaming));
            };
            Ok(RingEdge {
                node: NodeId::three_planes(a),
                to: NodeId::three_planes(b),
                carrier: Carrier::plane(wall),
                from_h: EndPin::Class(from_h),
                to_h: EndPin::Class(to_h),
            })
        })
        .collect()
}

/// A ring's edges from its nodes and the **carried** wall of each edge.
///
/// **Which plane of a point's name pins it on the line `a ∩ b`** — the one place that rule lives.
///
/// ★★★★★ **It was written five times before it was written once.** `ring_edges_with_walls`, the ray
/// caster's namer, and *both* arms of `arrangement::third_on_l` each spelled it, and the tracer's
/// seated road spelled a **reduced** version — no cut test, no smallest rule, a `panic` where this
/// returns `None`. That is this codebase's dominant defect shape: the correct rule inlined in a
/// sibling while a second site uses a smaller one. So it lives here now and they all call it.
///
/// **The rule, and why each half of it.** A handle only has to be *some* plane through the node
/// that **cuts** `a ∩ b` — that is all [`order_along`] asks of it, since it reads the handle through
/// `orient3d × dir_sign`. One *parallel* to the line names no point on it and would read `0` against
/// everything, fabricating a coincidence rather than missing one; hence
/// [`Judge::plane_pair_dir_sign`] `!= 0`. Under a concurrency several qualify and **any will do**,
/// so the **smallest** is taken and replay stays stable. `None` is "this name pins nothing here",
/// which each caller turns into its own vocabulary rather than sharing one label.
///
/// ★ `t` comes in sorted ([`NodeId::three_planes`] is the only constructor and it sorts), so
/// "smallest qualifying" is the first that qualifies — the two spellings the call sites used are
/// the same value.
///
/// ★★★★★ **Which half of this rule decides an answer, measured by breaking each.** Candidates
/// number `2` on 9,348 calls across the suite, so "smallest" is not vacuous *as a choice* — yet
/// taking the **largest** instead leaves the bit census **identical**. That is not a limp
/// instrument: it is the first measurement of the invariant this rule rests on, *"under a
/// concurrency several qualify and any will do"*, which until now was only asserted. What does
/// decide is the **cut test**: invert it and the census collapses from 269 rows to 5. So the
/// smallest rule buys replay stability, and the cut test buys correctness.
pub(crate) fn pin_on_line(
    jd: &Judge<'_, WorkingPlane>,
    a: usize,
    b: usize,
    t: [usize; 3],
) -> Option<usize> {
    t.iter()
        .copied()
        .find(|&c| c != a && c != b && jd.plane_pair_dir_sign(a, b, c) != 0)
}

/// ★ **This is the sound twin of `ring_from_names`.** That one reads an edge's supporting plane
/// back out of its two endpoint names — "the class they share besides `P`" — which works only while
/// every vertex lies on exactly three planes. Let four meet at a point, give it one canonical name,
/// and the name need not mention the plane the edge rides at all; the two names can even share
/// nothing but `P`. Measured (2026-07, `docs/dev-log.md`): a boolean whose fourth plane fell on an
/// arrangement vertex produced exactly that, and `ring_from_names` declined `RingNaming`.
///
/// The **handle** is still derived, and soundly — see [`pin_on_line`], where that rule now lives.
pub(crate) fn ring_edges_with_walls(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    nodes: &[NodeId],
    walls: &[usize],
) -> Result<Vec<RingEdge>, BoolError> {
    if nodes.len() != walls.len() {
        return Err(reject(RejectReason::RingNaming));
    }
    // ★ **The names are taken before the handles, and a branch node stops the ring here** — not
    // inside `handle`. Letting it fall into that `None` would report `RingNaming` ("these names do
    // not chain") for a vertex that has no three-plane name to chain in the first place: one label
    // for two causes. Collected through `Option<Vec<_>>` so a member cannot be dropped instead.
    let names: Vec<[usize; 3]> = nodes
        .iter()
        .map(|&n| three_plane_name(n))
        .collect::<Option<_>>()
        .ok_or_else(|| reject(RejectReason::BranchVertexUnnamed))?;
    let handle = |t: [usize; 3], wall: usize| pin_on_line(jd, p, wall, t);
    (0..nodes.len())
        .map(|i| {
            let j = (i + 1) % nodes.len();
            let wall = walls[i];
            let (Some(from_h), Some(to_h)) = (handle(names[i], wall), handle(names[j], wall))
            else {
                return Err(reject(RejectReason::RingNaming));
            };
            Ok(RingEdge {
                node: nodes[i],
                to: nodes[j],
                carrier: Carrier::plane(wall),
                from_h: EndPin::Class(from_h),
                to_h: EndPin::Class(to_h),
            })
        })
        .collect()
}

/// **A direction on a working plane**: the carrier that supplies it, and which way along it.
///
/// ★★★ **The fields are private, and that is the whole point.** Every rule that reads this
/// representation — [`turn`], [`antiparallel`], [`parallel_carriers`] via [`continuation`] — lives
/// in this module, and Rust's module privacy is what keeps it that way: a sibling module cannot
/// reach in even by accident. The previous cell had to catch fifteen such bypasses with a `rg`
/// gate because the identity it raised was an enum whose variants were nameable everywhere; here
/// the boundary is a compile error instead of a grep.
///
/// ★★ **A test that must read a direction is not a violation — it is an oracle.** It should read
/// **its own** inputs (the `(carrier, sense)` it built) rather than this type, both because that
/// keeps the boundary and because an oracle read back out of the value under test is an oracle
/// derived from its own answer. Do not add accessors for it.
///
/// ★★★ **`carrier`, not `wall` — and it grew when arcs arrived.** A straight edge's direction is a plane carrier and a sign;
/// an arc's is a **tangent at one end**, and the two ends differ. So the type is a sum, and the arc
/// arm carries *everything the comparison needs* — the end as `(line, s)`, the circle's centre, and
/// which way the axis points against the class's stored normal — rather than an index a reader
/// would have to resolve against a table. That is what keeps [`turn`] from needing one.
#[derive(Clone, Debug)]
pub(crate) enum EdgeDir {
    /// A straight edge: the plane whose meet with `P` carries it, `+1` when travel runs along
    /// `n_P × n_carrier`.
    Line { carrier: usize, sense: i8 },
    /// An arc **at one of its ends**. Boxed: the exact end is a `MeetLine` and a `QuadVal`, an
    /// order of magnitude wider than a line's two words, and every direction on the hot path is a
    /// line.
    Arc(Box<ArcDir>),
    /// A ruling: travel along `±m`, the cylinder's axis — the same at both ends, like a line, but
    /// with no plane-class carrier to name a `(carrier, sense)` pair by. Boxed for the axis
    /// vector's width.
    Ruling(Box<RulingDir>),
}

/// The payload of [`EdgeDir::Ruling`] — private fields like [`ArcDir`]'s, for the same reason.
#[derive(Clone, Debug)]
pub(crate) struct RulingDir {
    cyl: usize,
    /// Which of the two parallel rulings ([`RulingCarrier::side`]) — read by [`antiparallel`] and
    /// [`continuation`], where "one carrier" means one line, not one cylinder.
    side: i8,
    /// The axis direction `m`, rational — the one geometric fact [`turn`] needs.
    axis: [nacre_scalar::Rat; 3],
    /// Travel runs along `+m`.
    up: bool,
}

/// The payload of [`EdgeDir::Arc`] — private fields for the same reason the enum has them: every
/// rule that reads a direction lives in this module.
#[derive(Clone, Debug)]
pub(crate) struct ArcDir {
    cyl: usize,
    /// The end this direction is taken at, exactly.
    at: (nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal),
    /// The circle's centre on `P` — rational, because a circle bound's plane is ⊥ to the axis.
    centre: [nacre_scalar::Rat; 3],
    /// `n_P · m > 0`: the class's **stored** normal against the cylinder's axis.
    ///
    /// ★★ **It used to say "measured unexercised, `true` on every class the corpus reaches", and
    /// that went stale the moment a new population arrived.** ☑ Re-measured with the hole trace in:
    /// [`smooth_extremum_winding`] reads it `false` on half its firings, and dropping it there moves
    /// the chained fixtures' wall. What the old note still gets right is that no *fixture* forced
    /// the factor into being — the algebra did — and [`arc_side`] records which factors its own
    /// population locks.
    axis_up: bool,
    /// Travel runs counter-clockwise about the circle's own normal (the axis).
    ///
    /// ★ This one **is** locked: dropping it makes the two arcs at a crossing share a bucket in
    /// `arrangement::angular_order`, and the walk comes back `UnorderedEdges`.
    ccw: bool,
}

impl EdgeDir {
    /// A straight direction stated directly. Production makes them through [`dir_at`]; this exists
    /// for fixtures that state a `(carrier, sense)` pair by hand.
    pub(crate) fn new(carrier: usize, sense: i8) -> EdgeDir {
        EdgeDir::Line { carrier, sense }
    }

    /// The travel sense of a straight direction — `None` for an arc, whose direction is not a sign
    /// against a fixed carrier, and for a ruling, whose carrier is not a plane class (the one
    /// reader carries **plane** sub-segment senses forward).
    pub(crate) fn sense(&self) -> Option<i8> {
        match self {
            EdgeDir::Line { sense, .. } => Some(*sense),
            EdgeDir::Arc(_) | EdgeDir::Ruling(_) => None,
        }
    }
}

/// **An edge's direction on plane `p`** — its carrier is `wall` and its sense is `+1` when the edge
/// runs along `d = n_p × n_wall`, `-1` against it. The one place a direction is made.
///
/// ★ It returns the sense **paired with the carrier it belongs to** ([`EdgeDir`]) rather than a
/// bare `i8`, so the callers that used to do that pairing by hand no longer can. (`EdgeDir::new`
/// still states a pair directly — that is for fixtures, and it is the one place a wrong pairing is
/// still spellable.)
///
/// ★★ **It used to be made twice, two different ways.** `order_along` is `sign((V_i − V_j)·d)`, so
/// the direction of travel is either `−order_along(from, to)` (invert the result) or
/// `order_along(to, from)` (swap the arguments) — the same value by the antisymmetry of a
/// difference, and this file spelled it the first way while `arrangement`'s `angular_order`
/// spelled it the second, inline. They agreed only because two independent inversions happened to
/// cancel; a change to `order_along`'s convention would have moved one and not the other.
///
/// ★ **The zero policy lives here, and that is not a matter of taste**: both callers answered a
/// coincidence the same way (`CoincidentNodes`). Where two consumers want *different* answers —
/// the turn's zero, which is a rejection to one and a bucket to the other — the policy stays with
/// them and only the sign is shared (see [`turn`]).
///
/// ★★ **Both ends arrive as a name *and* a pin, because a cylinder-pinned one needs both.** The
/// pin says which kind of thing holds the point on `L`; the name says which point. A plane pin
/// carries its own name in its payload, so the two used to be one argument — but
/// [`EndPin::Cylinder`] has no payload by design ("the name is read from beside it"), and beside it
/// is here. [`order_pinned`] then picks the road.
pub(crate) fn edge_dir(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    wall: usize,
    from: (NodeId, EndPin),
    to: (NodeId, EndPin),
) -> Result<EdgeDir, BoolError> {
    // ★★ **`None` is `WitnessNotRational`, and the change of name is the change of proposition.**
    // This arm used to be `RingNaming` because the one thing it caught was a cylinder-pinned end —
    // a node with no third plane, which is genuinely a naming fact. That case is answered now, and
    // what is left is [`order_pinned`]'s own `None`: a class with no narrow rational description,
    // a coordinate past `Rat`. The names chain perfectly; the *value* could not be formed. That is
    // the sentence `arc_at` next door already uses for the same cause.
    let sense = match order_pinned(jd, cyls, p, wall, from, to) {
        Some(-1) => 1,
        Some(1) => -1,
        Some(_) => return Err(reject(RejectReason::CoincidentNodes)), // two nodes coincide
        None => return Err(reject(RejectReason::WitnessNotRational)),
    };
    Ok(EdgeDir::new(wall, sense))
}

/// **An edge's direction of travel, read at one of its ends** — and the only way a direction is
/// made.
///
/// ★★★★ **`node` is not decoration.** A straight edge's tangent is the same at both ends, so
/// "the edge's direction" and "the direction *at* this end" have always been one thing and no
/// caller had to say which it meant. An arc's two ends differ — and worse, a direction taken at a
/// node the edge does not even touch is a perfectly legal call that returns a confident, wrong
/// sign. That is not a hypothetical: it is what the first attempt at the arc walk did, and the
/// `turn == 0` it produced took a measurement to explain.
///
/// So the node comes in and is **checked here, at the one place a direction is born** — not at
/// `turn`, which would have to become fallible and drag `Result` through two sorts, and not at the
/// callers, who would each have to remember. `debug_assert` because "a ring is a chain" is a
/// *producer's* invariant, not an input's: the census runs in both profiles, so the debug one is
/// the instrument.
///
/// ★ **The hoist survives.** `angular_order` builds one direction per outgoing half-edge and sorts
/// on those values; building them inside the comparison instead would re-run `order_along` per
/// comparison — the very mistake this file measured and fixed once before (`split_at_crossings`'
/// `end_ds`: "1.5M `dir_sign` calls where 113k are distinct").
pub(crate) fn dir_at(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    e: &RingEdge,
    node: NodeId,
) -> Result<EdgeDir, BoolError> {
    debug_assert!(
        node == e.node || node == e.to,
        "a direction was asked at a node this edge does not touch — the ring is not a chain, \
         or the caller paired two edges that do not meet"
    );
    match &e.carrier {
        // ★ The carried sense first, and only then the endpoints. This used to say the endpoints
        // *cannot* answer a cut end — that was true while a branch point had no order, and
        // [`order_pinned`] has since given it one. What is left is the better reason: the sense is
        // the **splitter's own statement** about a piece it made, taken once for the whole segment
        // and handed to every sub-segment, rather than re-derived per piece from two names.
        // ★ The coincident-endpoint check `edge_dir` makes is not lost with it: the split refuses
        // two crossings at one parameter (`CoincidentNodes`) before any sub-segment is built, which
        // is the only way a carried sense can exist at all.
        Carrier::Plane {
            wall,
            sense: Some(s),
        } => Ok(EdgeDir::new(*wall, *s)),
        Carrier::Plane { wall, sense: None } => {
            edge_dir(jd, cyls, p, *wall, (e.node, e.from_h), (e.to, e.to_h))
        }
        Carrier::Arc(a) => arc_at(jd, p, a, node),
        // A ruling's direction is `±m` at both ends — the carrier states the travel, no endpoint
        // order is asked (its ends are branch points, which have no third plane to order by).
        Carrier::Ruling(r) => Ok(EdgeDir::Ruling(Box::new(RulingDir {
            cyl: r.cyl,
            side: r.side,
            axis: r.def.dir(),
            up: r.up,
        }))),
    }
}

/// **An arc's direction of travel at one of its ends** — the curved half of [`dir_at`].
///
/// Everything here is a *cache of the node's name*: [`branch_meet`] re-solves the point from the
/// name rather than taking a producer's coordinate, and the circle's centre is the axis point at
/// this plane's parameter. The direction itself is never materialized — [`turn`] reads it through
/// one `a + b√c` sign, and the three booleans below are what orient that sign.
fn arc_at(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &ArcCarrier,
    node: NodeId,
) -> Result<EdgeDir, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    // ★★ **Two causes, two names.** `branch_meet` folds four `None`s into one, and they are not the
    // same fact: a three-plane node on an arc is a *naming* failure (the split shipped the
    // two-names case it exists to refuse), while a class with no rational description, a meet that
    // does not solve, and a root the meet does not have are all "the exact route declined". Asking
    // the kind first is what keeps `RingNaming`'s sentence true where it is raised.
    //
    // ★ Only the first is reachable from here: the other three would have stopped the split that
    // built this arc, since it re-solves the same pair of classes for the same cylinder.
    let NodeId::Branch { .. } = node else {
        return Err(reject(RejectReason::RingNaming));
    };
    let at = branch_meet(jd, a.cyl, &a.def, node).ok_or_else(undecided)?;
    let coeffs = class_coeffs_rat(jd, p).ok_or_else(undecided)?;
    // The circle's centre: where the axis pierces this plane. ★ The canonical sign the coefficients
    // carry cancels in the parameter (numerator and denominator both flip), so this one does not
    // need the stored-frame turn that `arc_side` does.
    let t = crate::planes::axis_param_of_plane(&coeffs, &a.def).ok_or_else(undecided)?;
    let (o, m) = (a.def.origin(), a.def.dir());
    let centre: [nacre_scalar::Rat; 3] = (|| {
        let mut c = [nacre_scalar::Rat::from_int(0); 3];
        for k in 0..3 {
            c[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
        }
        Some(c)
    })()
    .ok_or_else(undecided)?;
    Ok(EdgeDir::Arc(Box::new(ArcDir {
        cyl: a.cyl,
        at,
        centre,
        axis_up: crate::planes::plus_t_is_above(&jd.planes[p], &a.def),
        ccw: a.ccw,
    })))
}

/// A class's exact description **turned to face its stored normal**.
///
/// ★★★ [`class_coeffs_rat`] hands back the class's *canonical* name — first nonzero component
/// positive — which points the **other way** from the stored normal on half the classes. Every
/// direction sign in this file is written in the stored frame ([`turn`]'s `frame_sign` bridge, the
/// cell labels' "above"), so a rule spelled against the canonical normal reads backwards on exactly
/// those classes and nowhere else. `planes::plus_t_is_above` carries the same warning and the
/// thirty-six tests that went red at once when it was not heeded.
///
/// The turn is decided by an `f64` dot of two **parallel** vectors — the class's own normal, twice
/// over — so the product is `±|a||b|`, a full magnitude from the sign boundary rather than a
/// near-zero comparison. That is `plus_t_is_above`'s argument verbatim, and it is why this is total
/// where the integer fold below is not.
///
/// ★ **And it is cross-checked against that integer fold.** `nacre_cip::predicate::name_stored_ints`
/// does the same turn exactly, in integers, on the classes it can (`None` for a wide name or a
/// witness that speaks another frame) — so where it answers, the two must agree.
pub(crate) fn stored_coeffs_rat(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
) -> Option<[nacre_scalar::Rat; 4]> {
    let coeffs = class_coeffs_rat(jd, c)?;
    let n = nacre_math::Vector3::from_array([
        coeffs[0].to_f64(),
        coeffs[1].to_f64(),
        coeffs[2].to_f64(),
    ]);
    let agrees = jd.planes[c].plane.normal().dot(n) > 0.0;
    debug_assert!(
        jd.planes[c].name_ints.as_ref().is_none_or(|ni| {
            let k = (0..4).find(|&k| coeffs[k] != nacre_scalar::Rat::from_int(0));
            // ★ Through `NameInts`' own accessor, not by touching the integers: production code
            // in this crate reaches `BigInt` only through `nacre-cip`'s types (the `num-bigint`
            // dependency is dev-only, and an assertion is not a reason to promote it).
            k.is_none_or(|k| ((ni.coeff_sign(k) > 0) == (coeffs[k].numer() > 0)) == agrees)
        }),
        "the canonical-to-stored turn disagrees with the class's integer fold"
    );
    if agrees {
        return Some(coeffs);
    }
    let zero = nacre_scalar::Rat::from_int(0);
    let mut out = [zero; 4];
    for k in 0..4 {
        out[k] = zero.checked_sub(coeffs[k])?;
    }
    Some(out)
}

/// The turn at ring node `i`, about the face's **outward** normal: `+1` left, `-1` right.
///
/// **No point is materialized, and no coordinate is read** — the algebra that makes that true
/// lives with the atom, [`turn`], and is not repeated here.
///
/// ★ **It is `0` exactly when the two walls name the same plane.** The node is `P ∩ A ∩ B`, and a
/// point exists there only if the three normals are independent — so with two genuinely different
/// planes the determinant cannot vanish. Two edges through one node whose lines are parallel are
/// two edges on **one** line, and that is the straight stretch a ring can arrive as (the
/// arrangement names a point wherever a feature crosses an edge); [`loop_winding`] walks back past
/// exactly that.
///
/// ★★ **"Same plane" is not "same class".** Aliasing lets two classes name one plane, which is why
/// `angular_order` folds them (`Aliases::union_wall`) and still keeps an honest `UnorderedEdges`
/// for the pair aliasing did not reach. (The older wording here said "never `0`" — a claim about
/// the predicate, where the truth is a claim about the two walls.)
///
/// A ring is not convex, so this is **not** the winding — at a reflex node it is its
/// opposite. [`loop_winding`] asks it at a hull vertex, where the two agree.
// Used by `loop_winding`'s tests and by the winding goldens in `lib.rs`; production reads the
// turn through `turn_between`, which lets the caller skip a straight stretch.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn turn_at(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
    i: usize,
) -> Result<i8, BoolError> {
    let n = ring.len();
    // ★ Both read **at this node** — the one they share. For a straight edge that is today's
    // answer either way; naming it is what makes an arc's two ends tellable apart.
    let arriving = dir_at(jd, cyls, p, &ring[(i + n - 1) % n], ring[i].node)?;
    let leaving = dir_at(jd, cyls, p, &ring[i], ring[i].node)?;
    turn_between(jd, p, &arriving, &leaving)
}

/// **Do two carriers give one direction on `P`?** — `P ∩ a` and `P ∩ b` are parallel.
///
/// ★ **`parallel`, not `collinear`, and the difference is the caller's.** Two parallel lines are
/// the same line only when they share a point. [`loop_winding`]'s walk-back has that extra premise
/// (its two edges are chained through a ring), so there parallel *does* mean collinear;
/// `arrangement`'s wall direction families do not — its own comment says so: *"Same family ⇒ the
/// two lines are parallel and meet in no point."* One predicate, two premises, and the premise
/// belongs to whoever has it.
///
/// ★★ **The same primitive answers a different question elsewhere and that is left alone.**
/// `plane_pair_dir_sign(p, wall, c) != 0` also spells "does this third plane *cut* the line" — a
/// fact about whether a **point** can be named on it, not about an edge's direction. Same
/// arithmetic, different sentence; folding them would put one name on two propositions.
pub(crate) fn parallel_carriers(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: usize,
    b: usize,
) -> bool {
    jd.plane_pair_dir_sign(p, a, b) == 0
}

/// **Two directions along one carrier, running opposite ways** — the `π` end of the angular order.
///
/// ★ **It is not "the angle is π"**, and the guard that makes it one stays with the caller:
/// `arrangement`'s `angular_order` asks this only inside the bucket where the turn already came
/// back `0`. Pulled out of that guard the predicate answers a different question, because two
/// directions on *different* carriers can also be collinear (aliasing that `union_wall` did not
/// reach), and those belong in the `0` bucket rather than the `π` one.
pub(crate) fn antiparallel(a: &EdgeDir, b: &EdgeDir) -> bool {
    match (a, b) {
        (
            EdgeDir::Line {
                carrier: ca,
                sense: sa,
            },
            EdgeDir::Line {
                carrier: cb,
                sense: sb,
            },
        ) => ca == cb && sa != sb,
        // ★ The two arcs of one circle meeting at one node are tangent by construction, so the
        // only question left is which way each travels. **They are at one node by construction**
        // too — a fan is built from a single vertex — which is why nothing here has to compare
        // positions (an earlier attempt used the sign of `s` as a stand-in for "same node", and
        // two different nodes can share it).
        (EdgeDir::Arc(x), EdgeDir::Arc(y)) => x.cyl == y.cyl && x.ccw != y.ccw,
        // One ruling, opposite travel — the straight reading of the arc arm above. Two *different*
        // rulings are parallel lines and share no node, so `(cyl, side)` identity is the "one
        // carrier" premise, same as `ca == cb` for lines.
        (EdgeDir::Ruling(x), EdgeDir::Ruling(y)) => {
            x.cyl == y.cyl && x.side == y.side && x.up != y.up
        }
        // ★★ **A line and an arc are never the π pole here, and the reason is upstream**: this is
        // asked only where the turn already came back `0`, and a `0` turn against an arc means the
        // segment is *tangent* at that node — while the split cuts only at **transversal**
        // crossings. So `false` is not a shrug: it sends the pair to the `0` bucket, where
        // `angular_order`'s `UnorderedEdges` names the surprise rather than ranking it. The same
        // sentence covers a line against a ruling (their `0` is a parallel-concurrency
        // degeneracy) and a ruling against an arc (they cannot share a node at all).
        _ => false,
    }
}

/// What an earlier ring edge does relative to `later` — the three answers [`loop_winding`]'s
/// walk-back needs, as one word each instead of two inline tests.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Continuation {
    /// Not parallel: the loop turns here, so this is the edge the turn is read against.
    Turns,
    /// Parallel and the same way: one straight run, keep walking back.
    Straight,
    /// Parallel and the *opposite* way: the ring doubles back along the line it came in on — an
    /// antenna, whose tip has no turn and whose neighbours' turn belongs to a different vertex.
    DoublesBack,
}

/// [`Continuation`] of `earlier` with respect to `later`, both on `P`.
///
/// ★★★ **`later` is the *fixed* reference, never the previous candidate.** The walk-back compares
/// every candidate with the edge the turn will be read at, not with its neighbour — chaining it
/// would weaken "the whole run goes one way" into "each neighbouring pair does", and those differ
/// on a run that reverses twice.
///
/// ★★ **Measured 2026-08-21: nothing in the suite reaches `DoublesBack`** — `straight_angle` is
/// raised nowhere at all (`--features reject-trace` over the workspace: 13 reasons, 48 raises, this
/// one zero), and it is not in the reject census's frozen corpus either. It is an unfired backstop
/// like `angular_order`'s `UnorderedEdges`, kept because upstream is *supposed* to make it
/// impossible (`merge_coincident`, `split_at_crossings`) and that is an argument rather than a
/// check. So this extraction is defended by derivation, not by a test — recorded, not hidden.
pub(crate) fn continuation(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    earlier: &EdgeDir,
    later: &EdgeDir,
) -> Continuation {
    match (earlier, later) {
        (
            EdgeDir::Line {
                carrier: ce,
                sense: se,
            },
            EdgeDir::Line {
                carrier: cl,
                sense: sl,
            },
        ) => {
            if !parallel_carriers(jd, p, *ce, *cl) {
                Continuation::Turns
            } else if se == sl {
                Continuation::Straight
            } else {
                Continuation::DoublesBack
            }
        }
        // Two arcs of one circle share a tangent at the node they meet at — the curved reading of
        // "parallel carriers". Same travel is a straight stretch; opposite is the antenna the
        // caller refuses.
        (EdgeDir::Arc(e), EdgeDir::Arc(l)) if e.cyl == l.cyl => {
            if e.ccw == l.ccw {
                Continuation::Straight
            } else {
                Continuation::DoublesBack
            }
        }
        // One ruling continuing through a node — the straight reading again, keyed like
        // `antiparallel`'s arm: `(cyl, side)` is the line's identity.
        (EdgeDir::Ruling(e), EdgeDir::Ruling(l)) if e.cyl == l.cyl && e.side == l.side => {
            if e.up == l.up {
                Continuation::Straight
            } else {
                Continuation::DoublesBack
            }
        }
        // A line and an arc at a **transversal** crossing turn — that is what transversal means,
        // and the split makes no other kind of node.
        _ => Continuation::Turns,
    }
}

/// The turn from the direction a loop **arrives on** to the one it **leaves on** — [`turn`] with
/// the straight-angle policy its ring callers share.
///
/// ★★ **It takes two directions, not two edges, and the point they meet at is the caller's.**
/// [`loop_winding`] reads this at `ring[lo].node`, and the arriving direction may come from an edge
/// that is *not* adjacent to it: the walk-back hands the edge on the far side of a straight run,
/// licensed by [`Continuation::Straight`] — the run is parallel and travelled the same way, so its
/// direction *is* the one the loop arrives on. That premise belongs to the caller that established
/// it, the same way [`antiparallel`]'s guard does.
///
/// ★ **Which end each direction was read at is [`dir_at`]'s to say** — a line's tangent does not
/// change along it, an arc's does, and that is why the direction is built from `(edge, node)`
/// rather than from an edge alone.
///
/// ★★ **The `0` refused here is the *straight* angle, and that is a plane sentence.** Two arcs of
/// one circle meeting at a node are tangent-continuous, so [`turn`] answers `0` there without
/// anything being wrong — the caller's walk-back reads that node as
/// [`Continuation::Straight`] and steps past it, which is why the `0` that reaches here is the
/// one where a loop doubles back on itself.
fn turn_between(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    arriving: &EdgeDir,
    leaving: &EdgeDir,
) -> Result<i8, BoolError> {
    match turn(jd, p, arriving, leaving)? {
        0 => Err(reject(RejectReason::StraightAngle)),
        t => Ok(t),
    }
}

/// **The signed turn between two directions at a point of plane `p`**, about that face's outward
/// normal — `+1` left, `-1` right, `0` collinear. The one place the sign is made.
///
/// Each direction is given as `(wall, dir)`: the plane it rides beside `p`, and [`edge_dir`]'s
/// sign along `n_p × n_wall`. **No point is materialized and no coordinate is read** — the two
/// directions are `s_a·(n_p × n_a)` and `s_b·(n_p × n_b)`, and
///
/// ```text
///   (n_p × n_a) × (n_p × n_b) = n_p · det[n_p, n_a, n_b]      (a×b)×(a×c) = a·det(a,b,c)
///     ⇒  turn = s_a · s_b · sign(det[n_p, n_a, n_b]) · orient_sign(p)
/// ```
///
/// where `sign(det[…])` is `plane_pair_dir_sign`, already exact. ★ The derivation is spelled
/// **here and nowhere else**: it used to sit in [`turn_at`]'s doc while the product itself was
/// written out twice, which is the shape this function exists to end.
///
/// ★★ **`0` comes back as `0`, on purpose.** Its two consumers want different things from it:
/// `turn_between` calls it a [`RejectReason::StraightAngle`], and `arrangement`'s `angular_order`
/// buckets it as the `0`/π pole and only rejects a *pair* of them. Deciding here would hand one
/// caller a rule that is not its own. By the same test, [`edge_dir`]'s zero *does* live in the
/// atom — both of its callers answer it identically.
///
/// ★★ **What arcs will widen, and what they will not.** A circle's tangent rides no plane, so an
/// arc arrives as a direction this signature cannot spell: what grows is the **input type**, not
/// the skeleton around it — the identity above becomes a cross product with one quadratic factor,
/// whose sign closes over `QuadVal::sign` (the tangent `m × (p − o)` is *linear* in the branch
/// point, so no biquadratic is needed; two tangents at one vertex would be quadratic, and that is
/// two circles meeting, already refused as `CylinderPairContact`).
///
/// ★★★ **What catches a mistake here, measured.** Negating this product fails **83** tests — but
/// *not* `arrangement`'s `angular_order_…_ccw` nor `a_reflex_node_turns_against_its_ring`, the two
/// that look like its unit goldens. Those build their fixture's ring with the same rule they then
/// read, so a global sign error flips twice and cancels: they pin relative structure, not the
/// convention. The convention is pinned end-to-end (volumes, cavities, nesting) — so a "tidy-up"
/// that drops the `frame_sign` factor will come back red, just not where a reader would look
/// first. (Dropping only `frame_sign` fails 11, all end-to-end — so the corpus does reach
/// `Reversed` faces. What that measures about the unit tests is narrower than it looks: **no unit
/// assertion is sensitive to that factor**, which is not the same as "no unit fixture reaches a
/// reversed face".)
///
/// ★★★ **Three more places read the direction's representation, and they are not this atom.** They
/// used to be inline — this paragraph was the only thing that found them, and it **undercounted**:
/// it listed two, and a sweep of `plane_pair_dir_sign`'s consumers turned up a third in another
/// file. Each is a named function beside this one now, so the next widening is a `match` the
/// compiler points at rather than a list a reader has to trust:
/// - [`antiparallel`] — the π pole in `angular_order` ("same wall, opposite sign"), which for arcs
///   becomes "same circle, opposite tangent";
/// - [`parallel_carriers`] — [`loop_winding`]'s walk-back, which becomes "are the tangents
///   parallel";
/// - [`parallel_carriers`] again — `arrangement::split_at_crossings`' **wall direction families**
///   (`Wall.dir`), the one the old list missed. It sits in a different file and asks the same
///   question with a weaker premise (its two walls share no point), which is why one predicate
///   serves both and the premise stays with the caller.
///
/// ★ **All three now answer for arcs, and the arity question is settled where it belonged.** For a
/// straight edge the tangent is the same at both ends, so "the edge's direction" and "the direction
/// *at* this endpoint" coincide and nothing ever had to tell them apart; an arc's two ends differ.
/// That is [`dir_at`]'s `node`, not a wider argument list here — the direction arrives already
/// bound to the end it was read at, so the sites above compare two directions and nothing else.
pub(crate) fn turn(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<i8, BoolError> {
    let frame = jd.planes[p].frame_sign;
    Ok(match (a, b) {
        (
            EdgeDir::Line {
                carrier: ca,
                sense: sa,
            },
            EdgeDir::Line {
                carrier: cb,
                sense: sb,
            },
        ) => sa * sb * jd.plane_pair_dir_sign(p, *ca, *cb) * frame,
        // ★★★ **A segment against an arc needs no new primitive, and the algebra says why.**
        // A circle bound's plane is ⊥ to the axis (`circle_on_class` answers `None` otherwise), so
        // `n_P ∥ m`. The arc's tangent is `T = ±(m × r)` with `r = x − c`, the segment's direction
        // `d = n_P × n_carrier` is ⊥ to `m`, and BAC-CAB collapses the cross product:
        //
        //   (d × T) · m = (d × (m × r)) · m = (m (d·r) − r (d·m)) · m = (d·r)(m·m)
        //
        // so the whole turn is `sign(d · (x − c))` — **one** `a + b√c` question, and
        // `quad::plane_side` is exactly it: the plane with normal `d` through the centre, measured
        // at the branch point.
        (EdgeDir::Line { carrier, sense }, EdgeDir::Arc(arc)) => {
            sense * arc_side(jd, p, *carrier, arc)?
        }
        (EdgeDir::Arc(arc), EdgeDir::Line { carrier, sense }) => {
            -sense * arc_side(jd, p, *carrier, arc)?
        }
        // Two arcs of one circle at one node are tangent: no turn to read. (Two *different*
        // circles cannot meet at a node — `cylinders_clear` keeps them `r₁+r₂` apart.)
        (EdgeDir::Arc(_), EdgeDir::Arc(_)) => 0,
        // ★★ **A segment against a ruling collapses the same way the arc arm did.** The ruling's
        // direction is `±m` and it lies *in* the class (`m · n_P = 0`), so BAC-CAB leaves
        //
        //   (d × m) · n_out = ((n_P × n_ca) × m) · n_out = (n_ca (m·n_P) − n_P (m·n_ca)) · n_out
        //                   = −(m · n_ca) (n_P · n_out)
        //
        // — one rational dot sign, times the same `frame` factor the line×line arm carries
        // (`n_P` here is the canonical spelling and the turn's reference is the face's outward
        // normal; their sign relation is `frame_sign`, entering **once** because `n_P` appears
        // once). `n_ca` appears once *and* once inside the sense's own definition
        // (`d = sense · (n_P × n_ca)`), so those two flips cancel — the arm needs only that
        // `sense` and the coefficients read the **same** `n_ca` ([`class_coeffs_rat`], the
        // canonical spelling `dir_sign` speaks).
        (EdgeDir::Line { carrier, sense }, EdgeDir::Ruling(r)) => {
            sense * ruling_line_turn(jd, p, *carrier, r)? * frame
        }
        (EdgeDir::Ruling(r), EdgeDir::Line { carrier, sense }) => {
            -sense * ruling_line_turn(jd, p, *carrier, r)? * frame
        }
        // All rulings on one class run along `±m`: parallel, no turn — the `0` bucket, where
        // `antiparallel` separates the π pole (same ruling, opposite travel).
        (EdgeDir::Ruling(_), EdgeDir::Ruling(_)) => 0,
        // A ruling and an arc cannot meet at a node: their classes demand the axis parallel and
        // perpendicular to `P` respectively, so the node would lie on two distinct cylinders —
        // which `cylinders_clear` keeps `r₁+r₂` apart. `0` sends a surprise to the bucket whose
        // walk names it (`UnorderedEdges`) rather than ranking it.
        (EdgeDir::Ruling(_), EdgeDir::Arc(_)) | (EdgeDir::Arc(_), EdgeDir::Ruling(_)) => 0,
    })
}

/// `turn(line, ruling)` with the line's sense factored out — the `−sign(m · n_carrier)` the
/// derivation above collapses to, times the ruling's travel.
///
/// ★ **The overall sign is measured** (cell 4): negating it turns the through-boss volume
/// oracles red (all five production fixtures) — the watcher the 2026-08-24 self-check said was
/// missing. It stood on the BAC-CAB derivation alone while the walk's both-handedness try
/// absorbed a global flip; the panel population reads it for real now.
fn ruling_line_turn(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    carrier: usize,
    r: &RulingDir,
) -> Result<i8, BoolError> {
    debug_assert!(
        class_coeffs_rat(jd, p)
            .map(|n| {
                nacre_scalar::dot_sign_rat(&[n[0], n[1], n[2]], &r.axis)
                    == nacre_scalar::Orient::Zero
            })
            .unwrap_or(true),
        "a ruling direction on a class its axis does not lie in"
    );
    // ★ The **stored** spelling, like every factor the walk's atoms read (`arc_side`'s
    // canonical→stored turn is measured-locked): `sense` is made against the stored frame, so
    // the carrier coefficients here must be too — the canonical name opposes it on half the
    // classes.
    let Some(n) = stored_coeffs_rat(jd, carrier) else {
        return Err(reject(RejectReason::WitnessNotRational));
    };
    let dot = nacre_scalar::dot_sign_rat(&[n[0], n[1], n[2]], &r.axis);
    let up = if r.up { 1 } else { -1 };
    Ok(match dot {
        nacre_scalar::Orient::Positive => -up,
        nacre_scalar::Orient::Negative => up,
        // `m · n_carrier = 0` means the segment's line is parallel to the ruling — no crossing
        // could have put them at one node, so a `0` here is the degenerate concurrency the
        // walk's `UnorderedEdges` bucket names.
        nacre_scalar::Orient::Zero => 0,
    })
}

/// `turn(line, arc)` with the line's sense factored out — the `sign(d·(x−c))` above, times the two
/// frame factors.
///
/// ★★ **It is a `Result`, and that is deliberate.** The first attempt returned `0` both for
/// "collinear" and for "the arithmetic ran out", and the callers read `0` as collinear — so a
/// width failure came back as a *shape* answer and the diagnosis took an extra measurement. The
/// two are different facts and this says so.
///
/// ★★★ **Which of its factors the corpus actually locks, measured by red probe** — because "the
/// suite is green" says nothing about a sign no fixture can see:
///
/// | factor | probe | verdict |
/// |---|---|---|
/// | the whole sign | negate the result | **invisible** — this same atom feeds both the cyclic order and the winding, and the walk tries both handednesses, so a *global* flip is absorbed by trying the other one |
/// | `ccw` | drop it | **locked** — the two arcs at a crossing collapse into one bucket, `UnorderedEdges` |
/// | the canonical→stored turn | use [`class_coeffs_rat`] | **locked** — one class in the corpus disagrees, and that class's walk merges four cells into one 8-half-edge orbit. ★ It used to be *accepted* there: the contour count passes it, and `arrangement::walk_cells`' Euler condition — added because of this probe — is what refuses it |
/// | `axis_up` | drop it | **unexercised** — `true` on every class reached |
/// | `frame_sign` | drop it | **unexercised** — `+1` on every class reached |
///
/// The last two are derived, not guessed (the algebra is above), and a corpus with a face whose
/// stored normal opposes its outward one, or a cylinder pointing the other way, is what would
/// close them.
///
/// ★★★ **The `sense` gap that used to be named here is CLOSED** (2026-08-22). Flipping the sense a
/// split carries onto its sub-segments (`Carrier::Plane::sense`) attaches the arcs to the wrong
/// cells, and nothing in the walk sees it — the cell count and the contour count both come out
/// right. Two guesses at its first reader were wrong in turn (`nest_cells`' root choice, then
/// `label_cells`' keep decision: every order-independent summary of the labels is identical
/// because the two 3-cells *swap* labels). It is `emit_faces`, where the ring comes out the exact
/// reverse, and `bands`' arc fence pins it there through `ClassAudit::outer_rings`.
///
/// ★ The table's first row was re-measured against that new lock and still holds: negating the
/// whole result leaves the whole crate green. A global flip really is absorbed.
fn arc_side(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    carrier: usize,
    arc: &ArcDir,
) -> Result<i8, BoolError> {
    use nacre_scalar::{Orient, Rat};
    let wide = || reject(RejectReason::WitnessNotRational);
    let ArcDir {
        at,
        centre,
        axis_up,
        ccw,
        ..
    } = arc;
    // ★★ **Stored-frame normals, not the canonical ones.** `d` is compared against `frame_sign`
    // below, which speaks the stored frame; `class_coeffs_rat` speaks the canonical one and points
    // the other way on half the classes. See [`stored_coeffs_rat`].
    let (Some(np), Some(nw)) = (stored_coeffs_rat(jd, p), stored_coeffs_rat(jd, carrier)) else {
        return Err(wide());
    };
    let d = cross3_rat(&[np[0], np[1], np[2]], &[nw[0], nw[1], nw[2]]).ok_or_else(wide)?;
    let plane = (|| -> Option<[Rat; 4]> {
        Some([
            d[0],
            d[1],
            d[2],
            Rat::from_int(0).checked_sub(dot3_rat(&d, centre)?)?,
        ])
    })()
    .ok_or_else(wide)?;
    let side = match nacre_scalar::quad::plane_side(&plane, &at.0, &at.1) {
        Orient::Positive => 1i8,
        Orient::Negative => -1,
        // The segment is **tangent** to the circle at this node. The split cuts only at
        // transversal crossings, so this is a shape the producer should not have made — but it is
        // still a *shape* answer, and `0` is what the callers read as "collinear".
        Orient::Zero => 0,
    };
    let up = if *axis_up { 1i8 } else { -1 };
    let way = if *ccw { 1i8 } else { -1 };
    Ok(side * up * way * jd.planes[p].frame_sign)
}

/// Face `f`'s outer-loop vertices as three-plane triples: `f`'s own plane, and the
/// neighbouring planes of the two edges meeting there.
///
/// An original vertex is as implicit a point as a seam node, so a cycle's ring — which mixes
/// them — is one uniform list and [`point_in_ring`] need not know the difference. Two
/// adjacent edges on one neighbour plane would be a straight angle, and it rejects.
#[allow(clippy::too_many_arguments)]
pub(crate) fn face_vertex_triples(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<LoopRing, BoolError> {
    loop_triples(model, &model.faces.get(f).outer, p, inc, jd, plane_ix, cyls)
}

/// One loop in class form with each edge's **carried wall** beside it: `walls[i]` is the plane
/// class the edge `triples[i] → triples[i+1]` rides — the far face's class, read off `inc` where
/// the triples are produced ([`loop_triples`]), never re-derived from the endpoint names. The
/// same trust model as the merge's `Ring { nodes, walls }`: deriving a wall from two names is
/// sound only while every vertex lies on exactly three planes, and the carried value is total
/// even where the *names* degenerate (the fallback-named vertices still know their edges).
///
/// ★★★★★ **A wall is a *carrier*, not a plane index** — the same [`crate::boolean::Wall`] the
/// result side's `Ring` uses, and deliberately not a second vocabulary. An operand can be a
/// previous boolean's result, and then a face's ring runs along a cylinder: a boss on a wall bites
/// an arc out of the plate's caps and splits the wall with its rulings. `usize` had nowhere to
/// write that, which is why the ring was declined there rather than described.
#[derive(Clone, Debug)]
pub(crate) struct NamedRing {
    pub triples: Vec<NodeId>,
    pub walls: Vec<crate::boolean::Wall>,
}

/// One loop of a face, in the vocabulary the tracer speaks (M6-2a): a polygon of three-plane
/// triples, or a **full circle** — one rim edge whose far face is a cylinder, named by that
/// cylinder's class. A circle has no triples, no walls and no endpoints; forcing it through
/// `NamedRing` was `RingNaming`'s job before the vocabulary existed.
#[derive(Clone, Debug)]
pub(crate) enum LoopRing {
    Poly(NamedRing),
    Circle { cyl: usize },
}

impl LoopRing {
    /// The polygon ring, `None` for a circle — the poly-only consumers' filter.
    pub(crate) fn poly(&self) -> Option<&NamedRing> {
        match self {
            LoopRing::Poly(nr) => Some(nr),
            LoopRing::Circle { .. } => None,
        }
    }
}

/// Every loop of one face in class form — **the only thing the tracer needs from `Model`**.
///
/// A loop that could not be named is `None` rather than an error, because the two failures have
/// *different names at the call site* (`OuterRing` vs `HoleRing`) and which one applies is the
/// tracer's to say, not this table's.
#[derive(Clone, Debug, Default)]
pub(crate) struct FaceLoops {
    /// The outer loop, or `None` if [`face_vertex_triples`] declined — and always `None` for a
    /// **lateral** face, whose outer loop is the two rims joined by a self-adjacent seam edge that
    /// no triple names (`crate::planes::CylFaceInfo::span` states that loop instead). Its `holes`
    /// beside it are named like any other face's.
    pub outer: Option<LoopRing>,
    /// One entry per hole ring, or `None` if [`hole_rings`] declined for **any** of them — a hole
    /// that cannot be named is not "no hole".
    pub holes: Option<Vec<LoopRing>>,
}

/// What one boolean's tracer reads instead of the `Model`: every face's loops, plus which slots
/// belong to which operand.
///
/// ★ **Both fields are independent of the plane being traced onto.** They were nevertheless
/// re-derived once per plane class — measured on an 80-fin fold at **2,280,285** calls (one per
/// face per class) for **13,492** distinct answers, 3.3% of the boolean. Hoisting them is what
/// makes the tracer a function of a face table rather than of a topology store, and the 169×
/// reduction comes along for free.
pub(crate) struct TraceInput {
    /// Each operand's faces: the `planes`-table slot and that face's loops, in the order the
    /// shells list them — the walk `trace_one` used to do over `Model`.
    ///
    /// **Compact on purpose.** This was a full-length `Vec<FaceLoops>` beside a list of slots — one
    /// row per table slot whether or not that slot's face was in the input. Pairing the slot with
    /// its loops makes the length the number of faces actually traced, which is what lets a caller
    /// hand the tracer a *subset* without the table's size leaking into the cost.
    pub faces: [Vec<(usize, FaceLoops)>; 2],
    /// **The population gate's own answer, carried — never re-derived** (M6-2 rulings ladder):
    /// the `(plane class, cylinder class)` pairs the gate let through **without** proving the
    /// class's faces clear of the lateral. The tracer's ruling and chord arms fire only on pairs
    /// listed here.
    ///
    /// ★★ **It said "always empty in production" and that went stale** — the gate-opening cell
    /// arrived and the record-and-pass arm fills it for a wall whose plane holds the axis exactly
    /// (`planes.rs`, `crossings.insert`). ☑ Re-measured: listed for the wall/boss pair 2 times in
    /// a plate-and-wall-boss fuse and 16 in the operation after it. A note that says what a
    /// population *is* expires when the population changes; this one is dated by its measurement
    /// instead.
    pub crossings: std::collections::HashSet<(usize, usize)>,
}

/// Derive [`TraceInput`] for one boolean, once.
///
/// Walked exactly as the tracer walked: shell by shell, `surf_ix` naming each face's slot. That is
/// what keeps the table's index space the `planes` one — a face missing from `surf_ix` cannot
/// happen, since `surf_ix` was built from the same two solids.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_input(
    model: &Model,
    operands: [(Handle<Solid>, &EdgeFaces); 2],
    surf_ix: &HashMap<Handle<Face>, usize>,
    n_faces: usize,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    crossings: std::collections::HashSet<(usize, usize)>,
) -> TraceInput {
    let _ = n_faces;
    let mut faces = [Vec::new(), Vec::new()];
    for (side, (solid, inc)) in operands.into_iter().enumerate() {
        for sh in crate::planes::solid_shell_handles(model, solid) {
            for &fh in &model.shells.get(sh).faces {
                let fp = surf_ix[&fh];
                // ★ A **lateral face's outer loop** is not named (M6-2a): it is the two rims joined
                // by the chart's seam, the seam edge is self-adjacent, and no triple describes its
                // corners — `CylFaceInfo::span` states that loop instead, and the tracer's cylinder
                // arm reads the row for it. So `outer` stays `None` here; it is a skip, not a
                // decline.
                //
                // ★★★★★ **Its holes are named like every other face's.** A fuse can burn a hole
                // into a band (a boss straddling a plate's wall), and that loop is an ordinary
                // closed ring whose corners are `plane ∩ plane ∩ cylinder` — the shape
                // [`branch_name_from_def`] restates. The tracer needs it to answer *per angle*
                // instead of claiming the whole circle, and naming it here is what puts the lateral
                // on the same road as everything else: one walk ([`ring_against_plane`]), not a
                // second description of the same loop.
                let loops = if matches!(plane_ix[fp], ClassIx::Cyl(_)) {
                    FaceLoops {
                        outer: None,
                        holes: hole_rings(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                    }
                } else {
                    FaceLoops {
                        outer: face_vertex_triples(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                        holes: hole_rings(model, fh, fp, inc, jd, plane_ix, cyls).ok(),
                    }
                };
                faces[side].push((fp, loops));
            }
        }
    }
    TraceInput { faces, crossings }
}

/// Each hole ring of face `f`, as three-plane triples.
///
/// A rim edge's incidence is `[p, wall]`, so the neighbour plane is the wall on the other
/// side of the rim — the same construction as an outer vertex, and the same rejection of a
/// straight angle. The ring keeps its stored direction: clockwise about `f`'s outward
/// normal, which is what makes it a hole.
#[allow(clippy::too_many_arguments)]
pub(crate) fn hole_rings(
    model: &Model,
    f: Handle<Face>,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Vec<LoopRing>, BoolError> {
    model
        .faces
        .get(f)
        .inner
        .iter()
        .map(|l| loop_triples(model, l, p, inc, jd, plane_ix, cyls))
        .collect()
}

/// A loop's vertices as names — a three-plane triple, or a branch point where a cylinder is one of
/// the three surfaces (the face's own, or a neighbour's).
///
/// The name normally comes from the loop itself — the face's own plane and the two neighbours the
/// meeting edges carry. **That fails when both neighbours lie on one plane**: an earlier boolean can
/// split a plane between two faces with opposite normals (a base's exposed top and the cantilever
/// underside above it), and a vertex where the loop runs straight through their shared line then
/// names one plane twice. Such a triple defines no point — `three_planes` answers `None`, and the
/// exact predicates, whose precondition is `D ≠ 0`, abort on it. Measured 2026-07-22: this is the
/// *only* path by which a degenerate triple reaches them.
///
/// So when the two neighbours are one class, the name is taken from **every plane touching the
/// vertex** ([`vertex_face_indices`]) instead: exactly three classes ⇒ that is the name, and it is
/// the same set whichever face's loop asks, so welding stays consistent (a face whose loop does not
/// degenerate here derives the same three). More than three is a real four-plane concurrency and
/// fewer is a genuine straight angle — both decline. Three *dependent* planes share a line rather
/// than a point, so independence is checked too ([`Judge::plane_pair_dir_sign`], which reads the same
/// un-normalized coefficients the consumer does).
///
/// The common case is untouched, so no existing name moves.
#[allow(clippy::too_many_arguments)]
fn loop_triples(
    model: &Model,
    l: &Loop,
    p: usize,
    inc: &EdgeFaces,
    jd: &Judge<'_, WorkingPlane>,
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
) -> Result<LoopRing, BoolError> {
    let hes = &l.half_edges;
    let edge = |he: &nacre_topo::HalfEdge| -> Result<([Handle<Vertex>; 2], [usize; 2]), BoolError> {
        inc.get(&he.edge)
            .copied()
            .ok_or_else(|| reject(RejectReason::MissingSeam))
    };
    let other = |pair: [usize; 2]| if pair[0] == p { pair[1] } else { pair[0] };
    let n = hes.len();
    // A **full-circle loop** (M6-2a): one rim edge whose far face is a cylinder — no triples,
    // no walls, no endpoints; it is named by the cylinder's class. Detected structurally
    // (`[v, v]` rims are the only single-edge loops a producer makes), so no curve is read.
    if n == 1 {
        let (_, pair) = edge(&hes[0])?;
        if let ClassIx::Cyl(k) = plane_ix[other(pair)] {
            return Ok(LoopRing::Circle { cyl: k });
        }
    }
    let mut out = Vec::with_capacity(n);
    let mut walls = Vec::with_capacity(n);
    for i in 0..n {
        // Vertex `i` starts edge `i` and ends edge `i - 1`.
        let (in_bounds, in_pair) = edge(&hes[(i + n - 1) % n])?;
        let (_, out_pair) = edge(&hes[i])?;
        let (a, b) = (other(in_pair), other(out_pair));
        // ★★★★★ **The filter [`crate::planes::ClassIx::plane`] asks its callers for.** That
        // accessor panics on a cylinder — "a loud panic beats a silently wrong plane" — and it is
        // right to, for the forty-odd callers whose upstream really does filter. This is the one
        // that has to *do* the filtering, because its input is an **operand**, and an operand can
        // be a previous boolean's result: a boss on a wall leaves the plate's caps bitten by an
        // arc and the wall split by two rulings, so those rings run along a cylinder. Both things
        // read below are plane-only — the carried wall and the vertex's three-plane name — so the
        // honest answer is to decline the ring rather than to name it wrongly or to abort.
        //
        // The one curved loop this road *does* speak is the full circle handled above.
        // ★★★★★ **A curved neighbour is described now, not declined.** The two questions are
        // **independent**: edge `i`'s carrier is `b`'s business, and the corner's name is both
        // neighbours' — a ruling can arrive at this corner and a plane leave it. What this road
        // lacked was the vertex's restatement into class space ([`branch_name_from_def`] — the
        // vertex already carries its own name) and somewhere to write a curved carrier
        // ([`NamedRing`]'s `Wall`). The one curved loop answered before this point is the circle.
        //
        // ★★★★★ **And the face itself may be the cylinder now.** A lateral face's *hole* is a loop
        // like any other — a rectangle of two arcs and two rulings in the chart — and its corners
        // are `plane ∩ plane ∩ cylinder`, the very shape [`branch_name_from_def`] restates. The
        // only loop of a lateral face this road still cannot walk is its **outer** one, whose seam
        // edge is self-adjacent (`other` gives back `p`) and whose corners are therefore not
        // three-surface points; `CylFaceInfo::span` states that loop instead.
        // ★★ **The corner is where half-edge `i` starts.** It used to be "the one vertex the two
        // edges share", which is the same vertex wherever that is unique — a loop's edge `i`
        // starts where edge `i − 1` ends (validate's `OpenLoop`) — and no vertex at all for a
        // two-gon (an arc and its chord share both ends). The half-edge already knows; asking the
        // two edges' bounds instead was a second spelling that failed on the one shape where the
        // first does not. `boolean` does not validate its inputs, so the invariant this leans on
        // is restated here, where it is leaned on.
        let corner = crate::he_start(model, hes[i]);
        debug_assert_eq!(
            if hes[(i + n - 1) % n].forward {
                in_bounds[1]
            } else {
                in_bounds[0]
            },
            corner,
            "a loop's edge starts where the previous one ends (OpenLoop)"
        );
        let branch = match plane_ix[p] {
            ClassIx::Plane(near) => match (plane_ix[a], plane_ix[b]) {
                (ClassIx::Cyl(k), ClassIx::Plane(far)) | (ClassIx::Plane(far), ClassIx::Cyl(k)) => {
                    Some(
                        branch_name_from_def(model, jd, corner, k, [near, far])
                            .ok_or_else(|| reject(RejectReason::CurvedOperandBoundary))?,
                    )
                }
                // Two laterals meeting at one corner is M6b's cylinder pair, not this road's.
                (ClassIx::Cyl(_), ClassIx::Cyl(_)) => {
                    return Err(reject(RejectReason::CurvedOperandBoundary));
                }
                (ClassIx::Plane(_), ClassIx::Plane(_)) => None,
            },
            ClassIx::Cyl(k) => match (plane_ix[a], plane_ix[b]) {
                (ClassIx::Plane(x), ClassIx::Plane(y)) => Some(
                    branch_name_from_def(model, jd, corner, k, [x, y])
                        .ok_or_else(|| reject(RejectReason::CurvedOperandBoundary))?,
                ),
                // A lateral's loop running along a second lateral is M6b's cylinder pair too.
                _ => return Err(reject(RejectReason::CurvedOperandBoundary)),
            },
        };
        // Edge `i`'s carried wall: the far face's class, read off `inc` — total even where the
        // vertex *names* below have to fall back or decline (see [`NamedRing`]).
        walls.push(match (plane_ix[b], plane_ix[p]) {
            (ClassIx::Plane(w), _) => crate::boolean::Wall::Plane(w),
            (ClassIx::Cyl(k), ClassIx::Plane(near)) => {
                let end = branch.expect("a curved edge's corner is a branch point");
                curved_wall(model, jd, cyls, &hes[i], k, near, end)?
            }
            (ClassIx::Cyl(_), ClassIx::Cyl(_)) => {
                unreachable!("a lateral face beside a lateral neighbour was rejected above")
            }
        });
        if let Some(n) = branch {
            out.push(n);
            continue;
        }
        let (ClassIx::Plane(near), ClassIx::Plane(wall), ClassIx::Plane(far)) =
            (plane_ix[p], plane_ix[b], plane_ix[a])
        else {
            unreachable!("the curved arms are handled above")
        };
        // `inc` names faces, so `other` matches by face — but the triple names *planes*, and a
        // consumer's `==` on it must mean "same plane". Canonize here, once, at the source.
        let mut t = [near, far, wall];
        t.sort_unstable();
        if t[0] != t[1] && t[1] != t[2] {
            out.push(NodeId::three_planes(t));
            continue;
        }
        // Both neighbours are one plane: the corner's own faces say what else names it.
        let mut classes: Vec<usize> = vertex_face_indices(corner, inc)
            .into_iter()
            .map(|k| plane_ix[k].plane())
            .collect();
        classes.sort_unstable();
        classes.dedup();
        if classes.len() != 3 {
            // >3: a genuine four-plane concurrency, which this substrate cannot name.
            // <3: the vertex lies on fewer than three planes, so no triple names it — a genuine
            // straight angle, where dropping the vertex would preserve the polygon but this brick
            // does not drop vertices.
            //
            // A vertex can also lack a triple by sitting on a **coplanar seam** (an edge between
            // two coplanar, same-facing faces of one solid): its only edges are the seam's, so it
            // touches one plane class. That is *not* a straight angle — the vertex is a corner and
            // dropping it would change the shape — but no producer makes such an edge since
            // `ImprintSketch` was retired, so it cannot reach here. See the 2026-07-22 dev-log
            // cells if one ever does: the answer is that a whole-seam ring is not a boundary.
            return Err(reject(if classes.len() > 3 {
                RejectReason::FourPlane
            } else {
                // Fewer than three classes: there is no triple to name this vertex by.
                RejectReason::RingNaming
            }));
        }
        // `plane_ix[p].plane()`, not `p`. An earlier revision kept `p` raw because consumers still matched the
        // face's own plane by raw index; they now compare classes (2026-07-22), and a triple that
        // mixed one face index with two class indices was exactly the ambiguity this brick exists
        // to remove.
        let mut t = [plane_ix[p].plane(), 0, 0];
        let mut k = 1;
        for &c in &classes {
            if c != plane_ix[p].plane() {
                t[k] = c;
                k += 1;
            }
        }
        if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
            return Err(reject(RejectReason::ThreePlanes)); // three planes through one line, not one point
        }
        out.push(NodeId::three_planes(t));
    }
    Ok(LoopRing::Poly(NamedRing {
        triples: out,
        walls,
    }))
}

/// The exact side of plane `q` that the implicit point `t` lies on: `0` means *on* it.
///
/// `+1` is the side the witness triangle's right-hand normal points to — that is `n_out(q)`, the face's
/// **outward** side, since `outer_tri` winds the triangle outward.
///
/// ★ **That is not the frame the arrangement's labels are stated in.** A plane class's
/// `[*_above, *_below]` labels are about the class root's **stored surface normal** — the convention
/// `SegKind::Seated{body_above}` and `emit_faces`' `flip` are written against — and the two frames
/// differ by [`crate::planes::FaceInfo::orient_sign`], which is `-1` exactly when the root face
/// is `Reversed`. No
/// `add_cuboid` face ever is, but a face an earlier boolean re-emitted flipped is (a pocket wall),
/// so **a producer that turns raw `side_of` into an above/below *label* silently flips its bit on
/// such a class**; multiply by `orient_sign(q)` if that is what you are computing. Reading a sign
/// *difference* (does this edge cross `W`?) is frame-free and needs no correction.
/// ★★ **`None` where the node is not three planes.** A [`NodeId::Branch`] *is* a point, but its
/// coordinates are quadratic-irrational and `orient3d` is the plane-triple judge — so this says
/// "not mine to answer" rather than guessing. Callers turn that into their own vocabulary (the
/// tracer a [`crate::DeclineKind`], the ray caster a reject), which is why it is not a reject here.
///
/// ★★★★★ **That arm used to be unexercised, and the rung that took `plane_ring`'s checks away
/// fired it — with a wrong sign.** The scan road now walks rings whose corners a cylinder made,
/// and the first thing that came back was a branch corner on the *opposite* side of its own face's
/// plane from its four plane-named neighbours (measured: `sides = [1, -1, -1, -1, -1, 1]`, the two
/// `1`s being the branch nodes). The scan read those as crossings that are not there and named one
/// with a plane triple that never met, and `orient3d` answered `D = 0`. ★ I reported that panic as
/// a hole upstream in the naming; every ring name measured correct and the fault was here.
pub(crate) fn side_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
    q: usize,
) -> Option<i8> {
    match n {
        NodeId::ThreePlane(t) => Some(jd.orient3d(t[0], t[1], t[2], q)),
        // ★★ **A branch point's side of a plane is one `a + b√c` sign.** The point is
        // `line.base() + s·line.dir()`, the plane's coefficients are rational, and
        // `quad::plane_side` is that sign — the same predicate `ruling_side` reads. `None` is a
        // missing description (a class with no world name, a cylinder with no world statement),
        // never a shape this cannot answer.
        NodeId::Branch { cyl, .. } => {
            let (line, sv) = branch_meet(jd, cyl, &cyls.get(cyl)?.def, n)?;
            let co = class_coeffs_rat(jd, q)?;
            let raw = jd.planes[q].plane.coefficients();
            // ★★★★★ **`world_rat` is the plane's *name*, not an oriented normal.** `orient3d`
            // answers against the class's **outward** normal, and `world_rat` may be any nonzero
            // multiple of the stored one — including a negative. Both describe one plane, so they
            // are proportional; the sign of that constant is read off the first component
            // `world_rat` makes nonzero, and `frame_sign` carries stored → outward.
            // ★ Both must be nonzero, not just the rational one: they are proportional so their
            // zero sets agree *exactly*, but `raw` is `f64` and a component it rounds to zero would
            // make `raw[i] > 0.0` false and hand back a sign with nothing behind it. Requiring both
            // turns that into a refusal.
            let zero = nacre_scalar::Rat::from_int(0);
            let i = (0..4).find(|&i| co[i] != zero && raw[i] != 0.0)?;
            let k = if (co[i] > zero) == (raw[i] > 0.0) {
                1
            } else {
                -1
            };
            let fix = k * jd.planes[q].frame_sign;
            Some(
                fix * match nacre_scalar::quad::plane_side(&co, &line, &sv) {
                    nacre_scalar::Orient::Positive => 1,
                    nacre_scalar::Orient::Negative => -1,
                    nacre_scalar::Orient::Zero => 0,
                },
            )
        }
    }
}

/// **What [`ring_against_plane`] found** — three outcomes, and they are three because collapsing
/// any two would put one name on unlike facts.
///
/// ★ `AllOn` used to be the walk's `None`, and a node it cannot read would have had to share it.
/// One says *the ring lies in the plane* (a shape), the other *this walk has no vocabulary for a
/// node* (a road) — and the three callers do different things with each.
pub(crate) enum RingWalk {
    /// Where the ring meets the line, in ring order from the first off-`q` node.
    Met(Vec<Feature>),
    /// Every node **and every edge** lies on `q`: a ring in the plane has no flanks to be decided
    /// by.
    ///
    /// ★★ **The second half is new and it is not pedantry.** This used to fire on "every node lies
    /// on `q`", which a ring with a curved edge satisfies while *leaving* the plane — and one
    /// consumer ([`every_ray`]) skips such a ring entirely, dropping its crossings from a parity
    /// count. A ring whose nodes are all on `q` while an edge departs is
    /// [`RingWalk::CurvedDeparture`] instead.
    AllOn,
    /// **The ring leaves the meet between two on-meet nodes, and which side it leaves to is not a
    /// question this walk can answer.**
    ///
    /// ★★★★★ Two nodes on `q` do not put the edge between them on `q`: two points fix a straight
    /// line, so a straight edge is on it and a **curved** one departs and returns. Where that
    /// happens the run is cut into on-meet stretches, which is enough for **extent** — but not for
    /// **parity's position**. The count of crossings over such a run is `flanks_differ` whichever
    /// side the arc bulges to (a circle meets a plane twice, so the departure is entirely on one
    /// side, and `(s_a ≠ σ) + (σ ≠ s_b)` has the parity of `s_a ≠ s_b` either way) — but *which* of
    /// the two ends carries it does depend on σ. So a departing run that must flip, a stretch with
    /// an arc on both sides, and a ring whose every node is on `q` while an edge departs are all
    /// answered by name rather than guessed.
    ///
    /// ★ The value it wants exists in pieces (`arc_side` beside `ccw` and the axis sense); it is
    /// not built because nothing asks yet, and an unmeasured sign is worse than an honest refusal.
    CurvedDeparture,
    /// A node whose side this walk cannot answer.
    ///
    /// ★★ **It used to mean "a [`NodeId::Branch`]", and it does not any more.** [`side_of`]'s
    /// branch arm answers, so a cylinder's corner is read like any other; what is left here is
    /// that arm's own `None` — a class with no narrow rational description, a cylinder missing
    /// from the table. ☑ Still never produced: measured **0** across the workspace suite.
    Unnameable,
}

/// Where a ring meets the line that `q` cuts its plane along — see [`ring_against_plane`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Feature {
    /// Edge `edge` (node `edge` → node `edge + 1`) crosses the line strictly inside: both its
    /// endpoints are off `q` and on opposite sides of it.
    ///
    /// `from` is the side the edge **leaves**, in [`side_of`]'s frame — never `0`, since a crossing
    /// has both ends off the line. See [`Self::Run`]'s `flank` for why a side travels with a
    /// feature at all.
    Crossing { edge: usize, from: i8 },
    /// `len` consecutive nodes from `first` lie *on* `q`, **and so do the edges between them** —
    /// an on-line interval rather than a point.
    ///
    /// ★★★★★ **That second half is a fact the walk now establishes, not one it used to assume.**
    /// It read "`len >= 2` means the edges between them lie on the line too", which is a claim about
    /// *nodes* being enough: two points fix a straight line, so a straight edge between two on-`q`
    /// nodes is on it — and a **curved** one departs and comes back, touching `q` only at its ends.
    /// A face's boundary that runs `… → node → arc → node → …` was read as one interval and stated
    /// a graze over ground it does not bound. The run is cut at each departure now
    /// ([`RingWalk::CurvedDeparture`] for what cutting cannot settle), so the sentence above is
    /// true again.
    ///
    /// `flanks_differ` is the whole decision: the two off-line neighbours bracketing the run sit on
    /// **opposite** sides, so the ring genuinely crosses the line here; equal sides mean it touched
    /// and turned back, and nothing crossed.
    ///
    /// ★★★ **`flank` is the side itself, and it travels here because asking again is a second
    /// walk.** A consumer that needs *which* side (not just whether the two agree) would otherwise
    /// call [`side_of`] over a ring member of its own — which is exactly the shape this walk exists
    /// to prevent, and a source-level lock says so
    /// (`rotation_sweep`'s `no_production_code_walks_a_ring_past_the_shared_walk`). It is the side
    /// of the off-line neighbour **before** the run; with `flanks_differ` false the one after is the
    /// same, and with it true the other is its negation, so one number carries both. Never `0`.
    ///
    /// ★★★★★ **`None` where there is no such neighbour** — a stretch that begins where the ring
    /// *returned* to the meet is preceded by the departure itself, and that side is σ
    /// ([`RingWalk::CurvedDeparture`]). The type says so rather than carrying a plausible number:
    /// the one reader today refuses on it, and the next one must face the same choice.
    Run {
        first: usize,
        len: usize,
        flanks_differ: bool,
        flank: Option<i8>,
    },
}

/// Read a ring against one plane: where it meets the line, and whether it crosses or only touches.
///
/// ★ **One walk, two consumers.** `arrangement::trace_transversal_face` clips a face's ring against a cut
/// plane and [`every_ray`] casts a parity ray along `P ∩ Q_a`; both must answer the same question
/// first — *does the boundary cross this line here?* — and a node sitting **on** the line is the
/// only hard part of it. The tracer had the rule (look at the node's two off-line neighbours:
/// opposite sides is a crossing, equal sides a touch) inlined in its scan, entangled with naming,
/// alias recording and decline kinds; the ray caster had no rule at all and threw such a candidate
/// away. Resilience that lives in one consumer is resilience the other does not have — the same
/// shape [`ring_in_ring`]'s retry was in.
///
/// What each consumer does *with* a feature stays its own: the tracer turns it into a **named
/// point** (four-plane aliases, `DeclineKind`, occupancy), the ray caster into one bit ("ahead of
/// `v`?"). That is where the sharing stops.
///
/// `AllOn` when every node **and edge** lies on `q` — a ring in the plane has no flanks to be
/// decided by.
///
/// **Features come out in ring order from the first off-`q` node.** That is the order the tracer's
/// scan produced them in, and its naming step records aliases into a union-find as it goes, so the
/// order is contract, not incident.
///
/// ★★ **`on_meet(i)` answers "does ring edge `i` lie on what `q` cuts here?"** — edge `i` runs from
/// `nodes[i]` to `nodes[i + 1]`. Two on-`q` nodes do **not** settle it: two points fix a straight
/// line, so a straight edge between them is on it, and a **curved** one leaves and comes back. What
/// counts as "the meet" is the caller's, because it differs by road — a plane class cuts a planar
/// face in a **line** (so the test is "not an arc") and a lateral face in a **circle** (so the test
/// is "carried by that very class", since an arc *can* lie on a circle).
///
/// ★ Today every caller answers `true` for every edge, which is exactly the claim the walk has been
/// making silently; the rung that reads it makes the claim honest.
pub(crate) fn ring_against_plane(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    nodes: &[NodeId],
    q: usize,
    on_meet: impl Fn(usize) -> bool,
) -> RingWalk {
    let n = nodes.len();
    let Some(side) = (0..n)
        .map(|i| side_of(jd, cyls, nodes[i], q))
        .collect::<Option<Vec<i8>>>()
    else {
        return RingWalk::Unnameable;
    };
    let Some(start) = side.iter().position(|&s| s != 0) else {
        // ★ Every *node* on `q` is not "the ring lies in `q`": an edge may still depart. Only when
        // every edge stays does `AllOn`'s sentence hold.
        return if (0..n).all(&on_meet) {
            RingWalk::AllOn
        } else {
            RingWalk::CurvedDeparture
        };
    };
    let mut out = Vec::new();
    let mut j = 0;
    while j < n {
        let i = (start + j) % n;
        if side[i] != 0 {
            let ni = (i + 1) % n;
            if side[ni] != 0 && side[ni] != side[i] {
                out.push(Feature::Crossing {
                    edge: i,
                    from: side[i],
                });
            }
            j += 1;
        } else {
            // A maximal run of on-line vertices. Two is the common case, but a vertex whose name
            // had to be taken from its touching planes (`loop_triples`) stays in the ring even
            // when the loop runs straight through it, so a run can be longer.
            let first = i;
            let mut len = 0usize;
            while j < n && side[(start + j) % n] == 0 {
                len += 1;
                j += 1;
            }
            let before = side[(first + n - 1) % n];
            let after = side[(start + j) % n];
            let flanks_differ = before != after;
            // ★★★★★ **The run is cut where the ring leaves the meet.** `len` on-`q` nodes give
            // `len - 1` edges; each that departs breaks the on-meet interval in two. A run with no
            // break is the whole thing, exactly as before.
            let mut breaks: Vec<usize> = (0..len.saturating_sub(1))
                .filter(|&k| !on_meet((first + k) % n))
                .collect();
            if breaks.is_empty() {
                out.push(Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    flank: Some(before),
                });
                continue;
            }
            // A break makes parity's *position* σ's question (see `CurvedDeparture`), and a stretch
            // with a departure on both sides has no flank of its own for the same reason.
            if flanks_differ {
                return RingWalk::CurvedDeparture;
            }
            breaks.push(len - 1);
            let mut from = 0usize;
            for cut in breaks {
                let stretch = cut + 1 - from;
                if stretch < 2 {
                    return RingWalk::CurvedDeparture;
                }
                out.push(Feature::Run {
                    first: (first + from) % n,
                    len: stretch,
                    flanks_differ: false,
                    flank: (from == 0).then_some(before),
                });
                from = cut + 1;
            }
        }
    }
    RingWalk::Met(out)
}

/// Is the implicit point `v` inside the simple ring `ring`, both on face plane `p`?
///
/// **A ray, cast along a line we already have.** Every ring edge lies on `P ∩ R`, and `v`
/// lies on `P ∩ Q_a` for either of its own two planes. Those two lines meet at
/// `X = {P, Q_a, R}`, which is *itself* a three-plane point — so "is `X` inside the edge"
/// and "is `X` ahead of `v`" are both [`order_along`], the comparator cell 3d already built
/// for two three-plane points on one line. **No coordinate is read and no point is built.**
///
/// **The flanks delete the special case, not the choice of ray.** A ring node *on* the ray's line
/// leaves no room for "is `X` inside the edge" to decide anything, and this used to abandon the
/// candidate; with both of `v`'s candidates abandoned the question came back `no_clear_ray`, and a
/// band of rotation angles died of it. But the node is not ambiguous at all — its two off-line
/// **neighbours** settle it: opposite sides and the boundary crossed here, equal sides and it
/// touched and turned back. That is the rule `trace_transversal_face` has always read a ring with,
/// and [`ring_against_plane`] is now where both get it.
///
/// So a `Feature::Run` — one node, or a whole edge of the ring lying on the line — contributes one
/// crossing iff its flanks differ, and a `Feature::Crossing` contributes one where it always did.
/// Nothing is counted twice: a crossing's endpoints are both off the line by construction.
///
/// Candidates are each node's two non-`P` planes, in ring order, `+d` before `-d`; the first
/// usable one wins, which keeps the answer deterministic. `no_clear_ray` survives for the two
/// cases nothing can name: a ring lying wholly on `Q_a` (no flanks), and an on-line node whose own
/// two walls are both parallel to the line (nothing pins it there). The answer must not depend on
/// which candidate was chosen, and a golden says so.
///
/// `v` must not lie *on* `ring` — a hole ring never touches the outer ring it sits in, and a seam
/// loop never touches `∂f` — and this is where it is finally checked: an intersection at
/// `X == v` strictly inside an edge is the `POINT_ON_RING` reject.
pub(crate) fn point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    v: [usize; 3],
    ring: &[RingEdge],
) -> Result<bool, BoolError> {
    every_ray(jd, p, v, ring)?
        .first()
        .copied()
        .ok_or_else(|| reject(RejectReason::NoClearRay))
}

/// Is the ring whose nodes are `probes` inside `outer` — **asked of the ring, not of one point**.
///
/// [`point_in_ring`] answers for a single vertex and can honestly fail there: every ray it can
/// cast from that vertex may have a ring node on its line. That is a fact about the *probe*, not
/// about the two rings, so the question is retried from the next node and only an exhausted ring
/// is a real degeneracy.
///
/// ★ **The retry belongs here, not in the callers.** It used to live in two of them, spelled two
/// different ways, and the third — the coplanar-merge cleaning pass — never got it: it probed
/// `nodes[0]` alone and rejected the whole boolean when that one vertex happened to be spoiled.
/// Measured, that lost a *band* of rotation angles at a stroke, because a small change of angle
/// leaves the topology (and therefore the unlucky first node) exactly where it was. Resilience
/// that lives in a caller is resilience the next caller does not have.
///
/// The probes are tried in order, so the answer is deterministic. Callers decide what *adjacency*
/// means for them (see `arrangement`'s `nest_cells` and `innermost_host`, which disagree) and ask
/// this only about rings they have already established are disjoint.
pub(crate) fn ring_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    probes: &[[usize; 3]],
    outer: &[RingEdge],
) -> Result<bool, BoolError> {
    for &v in probes {
        if let Ok(hit) = point_in_ring(jd, p, v, outer) {
            return Ok(hit);
        }
    }
    Err(reject(RejectReason::NoClearRay))
}

/// The parity every clear ray reports. The ring is simple, so they must all agree; a golden
/// says so, which is a second machine for free.
pub(crate) fn every_ray(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    v: [usize; 3],
    ring: &[RingEdge],
) -> Result<Vec<bool>, BoolError> {
    // The vertex name is a plane triple, so it obeys the same rule as a ring's: class roots only.
    // A caller holding face indices (a hand-built table, a test) is normalized here rather than
    // silently comparing a face against a class.
    let mut v = [v[0], v[1], v[2]];
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // ★★ **The unnameable check is the walk's now, so it is asked per `qa` rather than once up
    // front.** That is a real shift and it is spelled rather than glossed: a ring carrying a node
    // the walk cannot read used to be refused here even when no ray was cast, and is now refused
    // by the first ray that actually asks. The two differ only when *every* `qa` answers `AllOn`
    // — a ring lying in both cut planes — and that pairs with a branch node nothing produces here
    // yet, so the difference is unreachable twice over. Recorded because it will stop being.
    //
    // ★ A **ring**, not a probe list: `ring_against_plane` reads it as a cyclic sign sequence, so
    // a dropped member would be a different polygon answered about confidently — which is why the
    // walk is handed the ring **whole** and answers `Unnameable` for the ring rather than letting
    // a caller drop a node (see [`three_plane_probes`], where dropping *is* honest).
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let mut out = Vec::new();
    for &qa in v.iter().filter(|&&x| x != p) {
        // Where the ring meets the line — the walk `trace_transversal_face` reads too.
        //
        // ★ **A `Crossing` is exactly what this loop used to derive per edge.** It said "is
        // `X = {P, Q_a, R}` strictly inside the edge" with two `order_along`s: `a` is the sign of
        // `X − From` along `P ∩ R` and `b` that of `X − To`, both normalized to the same direction
        // on the same line, so `a·b < 0` iff `From` and `To` lie on opposite sides of `Q_a` — which
        // is what the walk already knows from their sides. The parallel guard goes with it: an edge
        // whose line is parallel to `P ∩ Q_a` has both endpoints on one side and is not a crossing.
        // ★ No cylinder table on this road: a branch node in a **result** cell's ring
        // declines here exactly as it did before, and threading one is the arc road's business.
        // ★ The meet here is a **line** (`p ∩ Q_a`), and two points fix a line — so a straight
        // edge between two on-line nodes is on it and a curved one is not. Same reading as
        // `arrangement::trace_transversal_face`'s, which walks the same kind of ring.
        let on_meet = |i: usize| !matches!(ring[i].carrier, Carrier::Arc(_));
        let features = match ring_against_plane(jd, &[], &nodes, qa, on_meet) {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::BranchVertexUnnamed)),
            RingWalk::AllOn => continue, // the whole ring lies on `Q_a`
            RingWalk::CurvedDeparture => return Err(reject(RejectReason::CurvedDeparture)),
        };
        let qb = *v
            .iter()
            .find(|&&x| x != p && x != qa)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        // A node on the line is a crossing point in its own right, so it must be nameable as one:
        // some plane of its own, off the line, pins it there. (`third_on_l` picks a handle the
        // same way, and for the same reason — one parallel to the line names no point on it.)
        // ★ The projection is safe *here* and nowhere earlier: the walk answered `Met`, which it
        // only does when every node is three planes.
        let namer = |i: usize| pin_on_line(jd, p, qa, three_plane_name(nodes[i])?);
        // Which side of `v` each run sits on, and whether the ring crossed the line there at all.
        // A run is one interval of the line and the ring is simple, so it cannot double back
        // inside itself: its two ends bracket it, and ends that disagree mean `v` is *between*
        // them — on the ring.
        let mut run_hits: Vec<i8> = Vec::new();
        let mut unnameable = false;
        for f in &features {
            let Feature::Run {
                first,
                len,
                flanks_differ,
                ..
            } = *f
            else {
                continue;
            };
            let ends = [first, (first + len - 1) % nodes.len()];
            let (Some(lo), Some(hi)) = (namer(ends[0]), namer(ends[1])) else {
                unnameable = true;
                break;
            };
            let o = [
                order_along(jd, p, qa, lo, qb),
                order_along(jd, p, qa, hi, qb),
            ];
            if o[0] != o[1] || o[0] == 0 {
                return Err(reject(RejectReason::PointOnRing)); // `v` inside the run, or one of it
            }
            run_hits.push(if flanks_differ { o[0] } else { 0 });
        }
        if unnameable {
            continue;
        }
        for dir in [1i8, -1] {
            let mut crossings = run_hits.iter().filter(|&&o| o == dir).count();
            for f in &features {
                let Feature::Crossing { edge, .. } = *f else {
                    continue; // runs are counted above
                };
                // Strictly ahead of `v` along `dir · (n_P × n_Qa)`?
                let carrier = ring[edge]
                    .carrier
                    .wall()
                    .ok_or_else(|| reject(RejectReason::RingNaming))?;
                match order_along(jd, p, qa, carrier, qb) {
                    0 => return Err(reject(RejectReason::PointOnRing)), // `X == v`, inside an edge
                    o if o == dir => crossings += 1,
                    _ => {}
                }
            }
            out.push(crossings % 2 == 1);
        }
    }
    Ok(out)
}

/// **Does the open segment between two named points on one line meet this face's material?**
///
/// [`every_ray`]'s sibling, and not its copy. That one asks about a **point** and may bail with
/// `PointOnRing` when the point lands on the boundary — here the two endpoints are *expected* to,
/// because the defect this answers is an edge whose ends sit on a face's ring while its middle
/// crosses the interior. An interval query has no candidate to fall back to, so every case that one
/// declines has to become a value.
///
/// ★ **One algorithm, no special cases.** The rings meet the line `P ∩ w` at a set of places; the
/// two unbounded ends of the line are outside the face and each genuine crossing flips that, so the
/// line reads **outside / inside / outside / …**. The answer is whether any *inside* stretch
/// overlaps the open `(u, v)`. Written as branches — "is there a crossing between them", "is it in a
/// hole" — it was twice wrong, because each branch re-derived a piece of that structure and lost
/// another. Read as one alternation it also subsumes the endpoint test: an endpoint strictly inside
/// puts `u` in an inside stretch.
///
/// **Holes come along for free.** All rings go into one bag: the rings of a face are disjoint, so
/// even-odd over the union *is* the material region (a point inside a hole has crossed twice) — the
/// rule `design.md` states one dimension down for 2D sketches.
///
/// ★★ **A `Run` is a stretch, not a place.** `Run { len >= 2 }` means the boundary *lies along* the
/// line, so it occupies an interval where the segment would be **on** the face rather than inside
/// it — two faces sharing an edge, which is ordinary adjacency. Those stretches are boundary and
/// are not counted as inside; `flanks_differ` still says whether crossing the run flips the side,
/// which is the rule the tracer and the ray caster already read a ring with.
pub(crate) fn segment_meets_face(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    w: usize,
    u: [usize; 3],
    v: [usize; 3],
    rings: &[Vec<RingEdge>],
) -> Result<bool, BoolError> {
    // A point on `P ∩ w` is named there by a third plane of its own that **cuts** the line; one
    // parallel to it names nothing (the same duty `every_ray`'s `namer` states).
    let handle = |t: [usize; 3]| -> Option<usize> {
        t.iter()
            .copied()
            .find(|&x| x != p && x != w && jd.plane_pair_dir_sign(p, w, x) != 0)
    };
    let (Some(hu), Some(hv)) = (handle(u), handle(v)) else {
        return Err(reject(RejectReason::RingNaming));
    };
    // Every place a ring meets the line: `[lo, hi]` handles (equal for a crossing at a point) and
    // whether passing it flips inside/outside.
    let mut events: Vec<([usize; 2], bool)> = Vec::new();
    for ring in rings {
        // ★ A ring, whole: the walk reads it as a cyclic sign sequence — see [`three_plane_probes`]
        // for where dropping a node *is* honest.
        let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
        let on_meet = |i: usize| !matches!(ring[i].carrier, Carrier::Arc(_));
        let features = match ring_against_plane(jd, &[], &nodes, w, on_meet) {
            RingWalk::Met(f) => f,
            RingWalk::Unnameable => return Err(reject(RejectReason::BranchVertexUnnamed)),
            // The whole ring lies on `w`: this face's boundary is the line itself, and the
            // alternation has no crossings to read. Refusing to guess.
            RingWalk::AllOn => return Err(reject(RejectReason::PointOnRing)),
            RingWalk::CurvedDeparture => return Err(reject(RejectReason::CurvedDeparture)),
        };
        for f in &features {
            match *f {
                Feature::Crossing { edge, .. } => {
                    let h = ring[edge]
                        .carrier
                        .wall()
                        .ok_or_else(|| reject(RejectReason::RingNaming))?;
                    events.push(([h, h], true));
                }
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    ..
                } => {
                    let ends = [first, (first + len - 1) % nodes.len()];
                    // ★ The projection is safe *here*: the walk answered `Met`, which it only does
                    // when every node is three planes.
                    let name = |i: usize| -> Option<usize> { handle(three_plane_name(nodes[i])?) };
                    let (Some(a), Some(b)) = (name(ends[0]), name(ends[1])) else {
                        return Err(reject(RejectReason::RingNaming));
                    };
                    let lo_first = order_along(jd, p, w, a, b) <= 0;
                    events.push((if lo_first { [a, b] } else { [b, a] }, flanks_differ));
                }
            }
        }
    }
    events.sort_by(|x, y| match order_along(jd, p, w, x.0[0], y.0[0]) {
        -1 => std::cmp::Ordering::Less,
        1 => std::cmp::Ordering::Greater,
        _ => std::cmp::Ordering::Equal,
    });
    // Walk the line: outside before the first event, flipping as each genuine crossing is passed.
    // The stretch between two events is a cell; an inside cell that overlaps the open `(u, v)` is
    // the surface meeting itself.
    let (lo, hi) = if order_along(jd, p, w, hu, hv) <= 0 {
        (hu, hv)
    } else {
        (hv, hu)
    };
    let mut inside = false;
    for i in 0..events.len() {
        if events[i].1 {
            inside = !inside;
        }
        if !inside {
            continue;
        }
        // The cell runs from this event's far end to the next event's near end.
        let cell_start = events[i].0[1];
        let Some(next) = events.get(i + 1) else {
            // ★ Reaching the unbounded tail while *inside* means the rings crossed the line an odd
            // number of times, which a closed curve cannot do. Reading it as "outside" would let a
            // real contact past in silence, so it is named instead — the same rule the rest of this
            // engine follows for an invariant it cannot verify.
            return Err(reject(RejectReason::RingParity));
        };
        let cell_end = next.0[0];
        // Overlap with the **open** interval: strictly, so touching at `u` or `v` is not inside.
        if order_along(jd, p, w, cell_start, hi) < 0 && order_along(jd, p, w, cell_end, lo) > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Every **face** incident to `vh`, as `planes`-table slots. (Was `vertex_plane_indices`; it
/// returns `inc`'s pairs verbatim, and those are faces. Its one caller maps them through
/// `plane_ix`.)
pub(crate) fn vertex_face_indices(vh: Handle<Vertex>, inc: &EdgeFaces) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::new();
    for (bounds, pair) in inc.values() {
        if bounds.contains(&vh) {
            for &p in pair {
                if !out.contains(&p) {
                    out.push(p);
                }
            }
        }
    }
    out.sort_unstable();
    out
}

/// Whether the point named by plane triple `v` lies on any edge of `ring` (a ring on plane `p`).
/// Used by [`point_in_component`] to abandon a non-generic ray rather than guess on a boundary.
fn point_on_ring(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    mut v: [usize; 3],
    ring: &[RingEdge],
) -> Result<bool, BoolError> {
    v.sort_unstable();
    if ring.len() < 3 {
        return Err(reject(RejectReason::DegenerateRing));
    }
    for e in ring {
        let (Some(si), Some(sj)) = (e.from_h.class(), e.to_h.class()) else {
            return Err(reject(RejectReason::RingNaming));
        };
        let r = e
            .carrier
            .wall()
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        // ☑ A literal triple, so `side_of`'s branch arm is unreachable here and the empty
        // cylinder table is never consulted — this asks about a *point*, not a ring.
        if side_of(jd, &[], NodeId::three_planes(v), r) != Some(0) {
            continue; // `v` is not even on the edge's line
        }
        // Name `v` as a point of that line: `{p, r, s}` for one of its own planes `s` off the line.
        let s = *v
            .iter()
            .find(|&&x| x != p && jd.plane_pair_dir_sign(p, r, x) != 0)
            .ok_or_else(|| reject(RejectReason::RingNaming))?;
        let (a, b) = (order_along(jd, p, r, s, si), order_along(jd, p, r, s, sj));
        if a * b <= 0 {
            return Ok(true); // between the endpoints (or on one)
        }
    }
    Ok(false)
}

/// Whether the point named by plane triple `query` lies inside a connected component — the exact
/// 3D lift of [`point_in_ring`]. A winding-parity ray whose line `L = a ∩ b` is built from two of
/// the query's own planes (never an arbitrary direction): each crossing with a face on plane `q`
/// is the exact three-plane point `{a,b,q}`, so the whole test is on the index-plane substrate —
/// no coordinate read, no f64. Non-convex is native (parity, not a convex test).
///
/// `faces` is the component as `(plane, rings)` per face — `rings[0]` = outer, `rings[1..]` = holes.
/// Only this component's faces are summed, so testing a void's vertex against a material component
/// reads `true` iff that material's outer shell nests the void.
///
/// A ray grazing a face boundary, or an undecidable in-face containment, is abandoned for the
/// query's next plane pair. **`Ok(None)` means every pair of this query's planes was blocked** —
/// this *node* cannot decide the question, which is a fact about the node, not an error: the
/// callers hold other nodes to try, and "every node abstained" is *their* proposition to reject
/// (`NoClearRay`, raised where the retries actually run out). `Err` is reserved for the
/// judgement itself failing (`JudgeExhausted`, a ring that cannot be named, …) — those must
/// propagate, never be traded for the next node: an abstention has other nodes as its remedy,
/// a failed judgement does not, and retrying it would let a real cause masquerade as
/// "no clear ray" once every node hit it.
/// One boundary of a component's face, with a polygon's edges already derived.
///
/// ★★★ **It mirrors [`crate::boolean::Bound`] on purpose.** The probe used to read a *flattened*
/// projection — `Vec<Vec<RingEdge>>`, which can only spell a polygon — so every circular and
/// banded boundary was **dropped on the way in** (`LocalFace::poly_rings`) and the component the
/// ray counted was not the component. Reading the engine's own boundary vocabulary is what lets
/// the ray count what is actually there.
#[derive(Clone, Debug)]
pub(crate) enum BoundEdges {
    /// A polygon, its edges carrying their walls (the only shape the ray counts today).
    Ring(Vec<RingEdge>),
    /// A whole circle, by the cylinder whose surface it rides — a disk face's outer bound, or a
    /// bored face's hole. Boxed for the same reason [`CompSurf::Cylinder`] is.
    Circle(Box<nacre_topo::CylinderDef>),
    /// A lateral band's boundary: the two plane classes its rims sit on.
    Band { lo: usize, hi: usize },
    /// A lateral face's boundary the ray does not count — a panel or hole ring on the cylinder,
    /// or a band with a chain rim. Every consumer abstains on it (`Ok(None)`), by name rather
    /// than by a projection the face does not have.
    Lateral,
}

/// A component face's surface, as the ray needs it.
///
/// ★ A cylinder carries its **truth**, not its class index: the crossings are solved against
/// `origin/dir/radius`, and the producer (`boolean`) is the one holding the class table. The
/// consumer should not have to look anything up — the same rule `RingEdge` states about walls.
#[derive(Clone, Debug)]
pub(crate) enum CompSurf {
    Plane(usize),
    /// Boxed because a `CylinderDef` is four rationals wide and every *planar* face would
    /// otherwise carry that much dead space — and planar faces are nearly all of them.
    Cylinder(Box<nacre_topo::CylinderDef>),
}

/// One face of a component — its surface, and its boundaries.
#[derive(Clone, Debug)]
pub(crate) struct CompFace {
    pub(crate) surf: CompSurf,
    pub(crate) outer: BoundEdges,
    pub(crate) inner: Vec<BoundEdges>,
}

/// How a ray met one cylindrical face.
enum CurvedHit {
    /// This many crossings lie on the **counted half** of the ray. ★ Which half that is belongs to
    /// the caller, not here: it hands in the plane, and this counts that plane's negative side.
    /// The two roads pick opposite halves — the named one counts behind its query, the coordinate
    /// one ahead of its origin — and the parity is the same either way, so a name that said
    /// "behind" would be false for one of them.
    Counted(usize),
    /// The ray touched a boundary — tangent, on a ruling, or exactly through the query. Abandon
    /// this ray and let the caller try the next plane pair, the same policy `point_on_ring` sets
    /// for a polygon.
    Graze,
}

/// **Is a face's material there**, given a way to ask one of its bounds — inside the outer bound
/// and outside every hole.
///
/// ★★★ **The combination is one rule; only the *asking* is two.** A named point answers a polygon
/// by a ring walk and a circle by its radial side; a rational one projects through a chart. Those
/// genuinely differ — the point arrives differently. But "outer and not any hole, and abandon the
/// moment a bound says *on me*" is the same sentence for both roads, and a sentence written twice
/// is one that drifts ([[rule-lives-inline-next-door]] is this repository's most-repeated defect).
/// It was written twice for one commit; this is that commit's correction.
///
/// ★★ **The hole clause fires and no test would notice if it stopped.** Measured: 7 entries with
/// a hole across `nacre-ops`, **6 of which subtract** — and stubbing the loop away leaves every
/// target green. That is a property of today's population, not of the rule: a hole here is a
/// **through** bore, so a ray that passes through it pierces the *pair* of annular caps and a
/// wrongly-counted crossing is wrongly counted **twice**, leaving the parity alone. A blind bore
/// on a classified component breaks the pairing and there is no such fixture yet. Extraction is
/// what guards this, not a lock.
fn material_of(
    f: &CompFace,
    mut ask: impl FnMut(&BoundEdges) -> Result<Option<bool>, BoolError>,
) -> Result<Option<bool>, BoolError> {
    let Some(mut here) = ask(&f.outer)? else {
        return Ok(None);
    };
    if here {
        for hole in &f.inner {
            match ask(hole)? {
                Some(true) => {
                    here = false;
                    break;
                }
                Some(false) => {}
                None => return Ok(None),
            }
        }
    }
    Ok(Some(here))
}

/// **A component's probe** — a point to ask "how deep is this component nested" from.
///
/// ★ The two variants are two *descriptions of a point*, not two policies: a three-plane name is
/// exact without coordinates at all (so it survives a rotated class, where no rational coordinate
/// exists), and rational coordinates are what a component describes when **none of its faces
/// carries a vertex** — a lone cylinder's boundary is two disks and a band.
///
/// ★★ **Both are points *on* the component's boundary**, which is what makes their depths
/// comparable: today's named probe is a face vertex, and [`coord_probes`] takes a cap disk's
/// **centre**, which lies on that face. An *interior* witness — the axis midpoint, say — would
/// also answer containment, but it is a different kind of point and could be separated from the
/// boundary by another component's wall.
pub(crate) enum Probe {
    Named([usize; 3]),
    Coord {
        p: [nacre_scalar::Rat; 3],
        dir: [nacre_scalar::Rat; 3],
    },
}

/// **Is this probe inside the component?** — the one door in front of the two roads.
pub(crate) fn probe_in_component(
    jd: &Judge<'_, WorkingPlane>,
    probe: &Probe,
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    match probe {
        Probe::Named(x) => point_in_component(jd, *x, faces),
        Probe::Coord { p, dir } => point_in_faces_rat(jd, p, dir, faces),
    }
}

/// **Coordinate probes of a component whose faces carry no vertex.**
///
/// The witness is a **cap disk's centre**: the face is planar with a circular outer bound, so the
/// centre is the axis point at that plane's own axis parameter — the rule `bands.rs` already reads
/// a class's position along an axis with ([`axis_param_of_plane`](crate::planes::axis_param_of_plane)).
/// It is strictly inside the circle for any positive radius, so nothing needs to test that.
///
/// ★ **A holed cap is passed over rather than guessed at.** The centre of an annulus is in its
/// hole, not on the face, and a witness that is not on the boundary is a confidently wrong depth
/// rather than an abstention. Another face — or, if there is none, the caller's reject — is the
/// remedy.
///
/// Several directions per point, because one ray can graze and the remedy is another direction;
/// the order is not load-bearing.
pub(crate) fn coord_probes(jd: &Judge<'_, WorkingPlane>, faces: &[CompFace]) -> Vec<Probe> {
    use nacre_scalar::Rat;
    let zero = Rat::from_int(0);
    let nonzero = |v: &[Rat; 3]| v.iter().any(|c| *c != zero);
    let neg = |v: &[Rat; 3]| -> Option<[Rat; 3]> {
        Some([
            zero.checked_sub(v[0])?,
            zero.checked_sub(v[1])?,
            zero.checked_sub(v[2])?,
        ])
    };
    let mut out = Vec::new();
    for f in faces {
        let (CompSurf::Plane(q), BoundEdges::Circle(def)) = (&f.surf, &f.outer) else {
            continue;
        };
        if !f.inner.is_empty() {
            continue;
        }
        let Some(p) = class_coeffs_rat(jd, *q)
            .and_then(|coeffs| crate::planes::axis_param_of_plane(&coeffs, def))
            .and_then(|t| {
                let (o, m) = (def.origin(), def.dir());
                let mut p = [zero; 3];
                for k in 0..3 {
                    p[k] = o[k].checked_add(t.checked_mul(m[k])?)?;
                }
                Some(p)
            })
        else {
            continue;
        };
        let m = def.dir();
        let mut dirs = vec![m];
        dirs.extend(neg(&m));
        for k in 0..3 {
            let mut e = [zero; 3];
            e[k] = Rat::from_int(1);
            match cross3_rat(&e, &m) {
                Some(w) if nonzero(&w) => {
                    dirs.push(w);
                    dirs.extend(neg(&w));
                    break;
                }
                _ => {}
            }
        }
        out.extend(dirs.into_iter().map(|dir| Probe::Coord { p, dir }));
    }
    out
}

/// **Two rational planes whose meet is the line through `p` along `dir`.**
///
/// ★ The normals are **basis crosses** (`ê_k × dir`), whose components are a shuffle of `dir`'s —
/// no products, nothing to overflow. That is the same rule `SketchPlane::normal_def` states, and
/// for the same reason: the "obvious" second direction `n × u` squares the inputs' widths.
fn planes_through_line(
    p: &[nacre_scalar::Rat; 3],
    dir: &[nacre_scalar::Rat; 3],
) -> Option<[[nacre_scalar::Rat; 4]; 2]> {
    use nacre_scalar::Rat;
    let zero = Rat::from_int(0);
    let basis = |k: usize| -> [Rat; 3] {
        let mut e = [zero; 3];
        e[k] = Rat::from_int(1);
        e
    };
    let mut ns: Vec<[Rat; 3]> = Vec::new();
    for k in 0..3 {
        if let Some(n) = cross3_rat(&basis(k), dir) {
            if n.iter().any(|c| *c != zero)
                && (ns.is_empty()
                    || cross3_rat(&ns[0], &n).is_some_and(|c| c.iter().any(|v| *v != zero)))
            {
                ns.push(n);
            }
        }
        if ns.len() == 2 {
            break;
        }
    }
    let [n0, n1] = <[[Rat; 3]; 2]>::try_from(ns).ok()?;
    let plane = |n: [Rat; 3]| -> Option<[Rat; 4]> {
        Some([n[0], n[1], n[2], zero.checked_sub(dot3_rat(&n, p)?)?])
    };
    Some([plane(n0)?, plane(n1)?])
}

/// **Is the point `p` inside this component**, asked along the ray `p + t·dir` — the coordinate
/// twin of [`point_in_component`].
///
/// ★★★ **Two roads, one crossing rule.** The two differ only in how the *point* arrives: a
/// three-plane name (exact without coordinates, so it survives a rotated class) or rational
/// coordinates (the only thing a curved component can offer, since none of its faces carries a
/// vertex). Everything they ask about a **face** — a circle's radial side, a band's roots and
/// axial span, the abandon-on-boundary policy — is the same rule and is called, not restated.
/// Restating it is how the two would come to disagree about a graze.
///
/// This road was `98e949a`'s casualty: K1 retired it when the band pass stopped needing it, and
/// its own note said it answered *"a 3D containment question — the shape the component probe
/// asks"*. It is back for that shape, with the curved arms it never had.
///
/// ★★ **What is measured, and what is owed.** `an_enclosed_cylindrical_void_is_a_cavity` is the
/// fixture that makes this road say `true`, and it goes red if the road is stubbed to `false` —
/// two *disjoint* bodies would not, since their answer is "outside" whatever the road does.
///
/// ★ The **curved** arm, though, is barely loaded from here, and the counts say so: this road is
/// entered **9 times** in the suite and **7 of those look at no cylinder face at all** (the other
/// component is all planes — the void fixture's is a box). The two that do are
/// `two_cylinders_with_coplanar_caps_fuse_apart`, and both are a **miss**. Nor can a fixture with
/// parallel axes do better: `cylinders_clear` refuses any boolean whose two cylinders are within
/// `r₁+r₂` of each other, so a ray from one cap's centre crosses the other band **0 or 2 times**
/// and the parity is the same with the arm stubbed to zero (measured). Its `k ≠ 0` behaviour is
/// load-bearing through [`point_in_component`], which is the reason it is *called* here rather
/// than restated: one arm, measured once. A fixture that loads it from **this** caller wants
/// **skew** axes and is owed.
///
/// `Ok(None)` = this ray grazed; the caller has other directions to try.
pub(crate) fn point_in_faces_rat(
    jd: &Judge<'_, WorkingPlane>,
    p: &[nacre_scalar::Rat; 3],
    dir: &[nacre_scalar::Rat; 3],
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
    use nacre_scalar::Rat;
    let not_rational = || reject(RejectReason::WitnessNotRational);
    let zero = Rat::from_int(0);
    if dir.iter().all(|c| *c == zero) {
        return Err(reject(RejectReason::DegenerateWitness));
    }
    // The half this road counts is **ahead** of the origin, so its cut plane's negative side is
    // `dir·(x − p) > 0` — the mirror of the named road, which counts behind. Either half has the
    // same parity; what matters is that one road picks one.
    let ahead = (|| {
        let n = [
            zero.checked_sub(dir[0])?,
            zero.checked_sub(dir[1])?,
            zero.checked_sub(dir[2])?,
        ];
        Some([n[0], n[1], n[2], zero.checked_sub(dot3_rat(&n, p)?)?])
    })()
    .ok_or_else(not_rational)?;
    let line = planes_through_line(p, dir).ok_or_else(not_rational)?;

    let mut count = 0usize;
    for f in faces {
        let q = match &f.surf {
            CompSurf::Cylinder(def) => {
                let BoundEdges::Band { lo, hi } = f.outer else {
                    return Ok(None);
                };
                match cylinder_face_crossings(jd, [&line[0], &line[1]], def, [lo, hi], &ahead) {
                    Some(CurvedHit::Counted(k)) => {
                        count += k;
                        continue;
                    }
                    Some(CurvedHit::Graze) | None => return Ok(None),
                }
            }
            CompSurf::Plane(q) => *q,
        };
        let coeffs = class_coeffs_rat(jd, q).ok_or_else(not_rational)?;
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let nd = dot3_rat(&n, dir).ok_or_else(not_rational)?;
        let residual = dot3_rat(&n, p)
            .and_then(|v| v.checked_add(coeffs[3]))
            .ok_or_else(not_rational)?;
        // ★ **Parallel is the origin's question, not the ray's.** `nd == 0` and a nonzero
        // residual is a plain miss. `nd == 0` with a zero residual is the ray lying **in** this
        // plane, where it never passes from one side of the face to the other — zero crossings —
        // *unless* the origin is on the face itself, and then "is the origin inside the component"
        // has no answer. That is exactly the `t == 0` policy below, so the two are one arm: a
        // coplanar ray is the origin's own crossing.
        let t = if nd == zero {
            if residual != zero {
                continue;
            }
            zero
        } else {
            (|| -> Option<Rat> {
                zero.checked_sub(residual)?
                    .checked_mul(Rat::new(nd.denom(), nd.numer())?)
            })()
            .ok_or_else(not_rational)?
        };
        if t < zero {
            continue; // behind the origin
        }
        let x = (|| -> Option<[Rat; 3]> {
            let mut x = *p;
            for k in 0..3 {
                x[k] = x[k].checked_add(t.checked_mul(dir[k])?)?;
            }
            Some(x)
        })()
        .ok_or_else(not_rational)?;
        let chart = Chart2dRat::of_normal(&n).ok_or_else(not_rational)?;
        let x2 = chart.project(&x).ok_or_else(not_rational)?;
        // Each bound answers by its own kind — the same split the named road makes.
        let inside = |b: &BoundEdges| -> Result<Option<bool>, BoolError> {
            match b {
                BoundEdges::Ring(r) => {
                    // ★ Under three nodes `point_in_ring_2d_rat` answers `Outside` by contract,
                    // which would make a degenerate ring *invisible* to the parity instead of
                    // loud. A face of a valid solid has no such ring, so saying so is free.
                    if r.len() < 3 {
                        return Err(reject(RejectReason::DegenerateRing));
                    }
                    // The branch check is spelled before the chart rather than left to
                    // `node_coords_rat`'s `None`, so an unnamed vertex is reported as itself and
                    // not as arithmetic that ran out of room — the same split `every_ray` makes.
                    let nodes: Vec<NodeId> = r
                        .iter()
                        .map(|e| three_plane_name(e.node).map(NodeId::ThreePlane))
                        .collect::<Option<_>>()
                        .ok_or_else(|| reject(RejectReason::BranchVertexUnnamed))?;
                    let ring2 = chart.ring(jd, &nodes).ok_or_else(not_rational)?;
                    Ok(match point_in_ring_2d_rat(x2, &ring2) {
                        RingSide::Inside => Some(true),
                        RingSide::Outside => Some(false),
                        RingSide::OnBoundary => None,
                    })
                }
                BoundEdges::Circle(def) => Ok(point_in_disk(&x, def)),
                BoundEdges::Band { .. } | BoundEdges::Lateral => Ok(None),
            }
        };
        let Some(material) = material_of(f, inside)? else {
            return Ok(None);
        };
        if t == zero {
            // The crossing is the ray's own origin — either the ray meets this plane there, or it
            // lies in it. On the face's material the query sits on the component's boundary, where
            // "inside" has no answer; off it, the plane is touched (or run along) at points
            // outside the face and counts for nothing.
            if material {
                return Ok(None);
            }
            continue;
        }
        if material {
            count += 1;
        }
    }
    Ok(Some(count % 2 == 1))
}

/// **Is `p` strictly inside this circle?** — `None` when it lies *on* the circle (non-generic,
/// abandon) or when the arithmetic could not answer.
///
/// ★ It takes **coordinates**, so both roads ask it: the named probe realizes its three-plane
/// crossing first, the coordinate road already has one. The rule must not be written twice.
///
/// ★ The rule is `cylinder_radial_side`'s, the one `arrangement::node_in_circle` already reads —
/// a circle bound is `cylinder ∩ plane`, so "inside the disk" is "inside the cylinder's radius".
fn point_in_disk(p: &[nacre_scalar::Rat; 3], def: &nacre_topo::CylinderDef) -> Option<bool> {
    match nacre_scalar::cylinder_radial_side(p, &def.origin(), &def.dir(), def.radius()) {
        nacre_scalar::Orient::Negative => Some(true),
        nacre_scalar::Orient::Positive => Some(false),
        nacre_scalar::Orient::Zero => None,
    }
}

/// [`cylinder_face_crossings`] for a probe ray named by three plane classes: it states the ray's
/// two planes and the "behind" cut plane from the class table and hands them over.
///
/// `Ok(None)` — as `Some(None)` here — is the abandon-this-ray answer; `Err` never happens because
/// every failure this can meet is arithmetic, and arithmetic that cannot answer is an abstention.
fn curved_count(
    jd: &Judge<'_, WorkingPlane>,
    a: usize,
    b: usize,
    c: usize,
    def: &nacre_topo::CylinderDef,
    span: [usize; 2],
) -> Result<Option<usize>, BoolError> {
    // ★ The cylinder gate refuses any class that is rotated or has no narrow rational name
    // (`planes.rs`), so in a boolean that has a cylinder these are always `Some`. A `None` here
    // would mean that gate let something through — assert it, then abstain rather than guess.
    let Some(((ca, cb), cc)) = class_coeffs_rat(jd, a)
        .zip(class_coeffs_rat(jd, b))
        .zip(class_coeffs_rat(jd, c))
    else {
        debug_assert!(
            false,
            "the cylinder gate is supposed to make these rational"
        );
        return Ok(None);
    };
    // The cut plane through the ray's origin, normal `n_a × n_b` — the very direction
    // `plane_plane_cylinder` gives its meet line, so "negative side" is "behind the query".
    let Some(half) = (|| {
        let d = cross3_rat(&[ca[0], ca[1], ca[2]], &[cb[0], cb[1], cb[2]])?;
        let at = nacre_scalar::three_planes_rat([ca, cb, cc])?;
        let d0 = nacre_scalar::Rat::from_int(0).checked_sub(dot3_rat(&d, &at)?)?;
        Some([d[0], d[1], d[2], d0])
    })() else {
        return Ok(None);
    };
    Ok(
        match cylinder_face_crossings(jd, [&ca, &cb], def, span, &half) {
            Some(CurvedHit::Counted(k)) => Some(k),
            Some(CurvedHit::Graze) | None => None,
        },
    )
}

/// **Where the line `line[0] ∩ line[1]` crosses one lateral band, and how many of those lie on
/// the half the caller is counting.**
///
/// ★★★ **One copy, two roads.** The named-point probe below and the coordinate-point road both
/// state their ray as *two rational planes*, so both ask this. Writing the arms twice is how the
/// two would come to disagree about a graze.
///
/// `half` is a rational plane through the ray's origin, and a crossing counts exactly when it is
/// on that plane's **negative** side — so each road states its own half and nothing here has a
/// front or a back:
///
/// - the named road's normal is the line's own direction `n_a × n_b`, which counts **behind** the
///   query. ★ That is not a convention to be guessed: `plane_plane_cylinder` builds the meet
///   line's `dir` as `cross3(n_a, n_b)`, the very `d` that `order_along` — and therefore that
///   probe's `fwd == 1` — measures against.
/// - the coordinate road's normal is `−dir`, which counts **ahead** of the origin.
///
/// `span` are the two rim plane classes. The test for "on this band" does **not** read their
/// stored normals: it takes each rim's axis parameter ([`crate::planes::axis_param_of_plane`]) and
/// rebuilds a plane there with the **axis** as normal, so "between the rims" is simply "opposite
/// sides of two planes that point the same way". Reading the stored normals is the mistake
/// `plus_t_is_above`'s doc records as having turned 36 tests red.
///
/// `None` is checked-`Rat` arithmetic that could not answer — an honest decline, never a guess.
fn cylinder_face_crossings(
    jd: &Judge<'_, WorkingPlane>,
    line: [&[nacre_scalar::Rat; 4]; 2],
    def: &nacre_topo::CylinderDef,
    span: [usize; 2],
    half: &[nacre_scalar::Rat; 4],
) -> Option<CurvedHit> {
    use nacre_scalar::Orient;
    use nacre_scalar::quad::{CylinderMeet, QuadVal};
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let (meet, roots): (_, [QuadVal; 2]) =
        match nacre_scalar::quad::plane_plane_cylinder(line[0], line[1], &o, &m, r)? {
            CylinderMeet::Pair { line, s } => (line, s),
            // A tangency touches the surface without crossing it — parity has no answer there.
            CylinderMeet::Tangent { .. } => return Some(CurvedHit::Graze),
            // ★ The ray lies **on** the cylinder: every point of it is on the face. Not a crossing
            // count at all.
            CylinderMeet::OnRuling(_) => return Some(CurvedHit::Graze),
            // Parallel to the axis and off the surface, or missing the quadric outright: no crossing.
            CylinderMeet::AxisParallelMiss(_) | CylinderMeet::Miss(_) => {
                return Some(CurvedHit::Counted(0));
            }
            // ★ Two classes can name one plane (aliasing), and then `a ∩ b` is not a line — there is
            // no ray to count along. `circle_crossings` proves these away for a class ⊥ to the axis;
            // a probe ray has no such proof, so they are answered rather than `unreachable!`d.
            CylinderMeet::CoincidentPlanes | CylinderMeet::ParallelPlanes => {
                return Some(CurvedHit::Graze);
            }
        };
    // The two rims, restated with the **axis** as their normal so "between" is a sign difference.
    let mut rim = Vec::with_capacity(2);
    for c in span {
        let coeffs = class_coeffs_rat(jd, c)?;
        let t = crate::planes::axis_param_of_plane(&coeffs, def)?;
        // n = m, d0 = −m·(o + t·m)
        let mut d0 = nacre_scalar::Rat::from_int(0);
        for k in 0..3 {
            let at = o[k].checked_add(t.checked_mul(m[k])?)?;
            d0 = d0.checked_sub(m[k].checked_mul(at)?)?;
        }
        rim.push([m[0], m[1], m[2], d0]);
    }
    let mut count = 0usize;
    for s in &roots {
        let sides = [
            nacre_scalar::quad::plane_side(&rim[0], &meet, s),
            nacre_scalar::quad::plane_side(&rim[1], &meet, s),
        ];
        // On a rim: the crossing is on the face's own boundary — non-generic, abandon.
        if sides.contains(&Orient::Zero) {
            return Some(CurvedHit::Graze);
        }
        if sides[0] == sides[1] {
            continue; // outside the band's axial span — this face is not crossed there
        }
        match nacre_scalar::quad::plane_side(half, &meet, s) {
            // The crossing **is** the ray's origin, and it is on this face: the query sits on the
            // other component's surface, which the parity cannot answer. Same as the planar
            // `fwd == 0 && in_g` arm.
            Orient::Zero => return Some(CurvedHit::Graze),
            Orient::Negative => count += 1,
            Orient::Positive => {}
        }
    }
    Some(CurvedHit::Counted(count))
}

/// A component as its faces, each boundary already carrying its edges' walls.
pub(crate) type ComponentFaces = Vec<CompFace>;

pub(crate) fn point_in_component(
    jd: &Judge<'_, WorkingPlane>,
    query: [usize; 3],
    faces: &[CompFace],
) -> Result<Option<bool>, BoolError> {
    let mut vplanes = query.to_vec();
    vplanes.sort_unstable();
    vplanes.dedup();

    // One ray attempt along `L = a ∩ b`, located by the query point `{a,b,c}`. `Ok(None)` = the
    // line grazed and the caller should try the next plane pair.
    let attempt = |a: usize, b: usize, c: usize| -> Result<Option<bool>, BoolError> {
        let mut count = 0usize;
        for f in faces {
            // ★ A cylindrical face is counted by its own arm — the crossings are roots of a
            // quadratic, not three-plane points, and "inside the face" is an axial span rather
            // than a ring walk.
            if let CompSurf::Cylinder(def) = &f.surf {
                let BoundEdges::Band { lo, hi } = f.outer else {
                    // A cylinder face bounded by anything but a band has no producer; refusing to
                    // guess costs the caller another node, never a wrong answer.
                    return Ok(None);
                };
                match curved_count(jd, a, b, c, def, [lo, hi])? {
                    Some(k) => {
                        count += k;
                        continue;
                    }
                    None => return Ok(None),
                }
            }
            let CompSurf::Plane(q) = &f.surf else {
                return Ok(None);
            };
            let q = *q;
            // Inside `q`'s material at the point `x`: inside the outer bound, outside every hole.
            // ★ Each bound answers by its own kind — a polygon by a ring walk, a circle by its
            // radial side — and either can say "the point is *on* me", which is the same
            // non-generic abandon `point_on_ring` has always raised.
            let material = |x: [usize; 3]| -> Result<Option<bool>, BoolError> {
                let inside = |b: &BoundEdges| -> Result<Option<bool>, BoolError> {
                    match b {
                        BoundEdges::Ring(r) => {
                            if point_on_ring(jd, q, x, r)? {
                                return Ok(None);
                            }
                            Ok(every_ray(jd, q, x, r)?.first().copied())
                        }
                        BoundEdges::Circle(def) => Ok(node_coords_rat(jd, NodeId::three_planes(x))
                            .and_then(|p| point_in_disk(&p, def))),
                        // A band bounds a cylinder, never a plane — a producer error, not an input.
                        BoundEdges::Band { .. } | BoundEdges::Lateral => Ok(None),
                    }
                };
                material_of(f, inside)
            };
            if jd.plane_pair_dir_sign(a, b, q) == 0 {
                // `L` is parallel to `q` — and possibly **in** it, which is not the same thing.
                // A coplanar ray never passes from one side of this face to the other, so zero
                // crossings is the right count and always was. What the old `continue` also
                // swallowed is the *query*: if it lies on this face, "is the query inside the
                // component" has no answer at all, and skipping the face answers it anyway. That
                // is the same proposition the `fwd == 0` arm below abandons for a transversal
                // plane — one rule that had only one of its two spellings.
                //
                // ★ **Unfired, and measured to be.** 26 rays in the suite lie in a face's plane
                // and **none** of them is on that face's material. (An earlier count said two;
                // it read the *outer bound* rather than the face, and both were in a hole.) It is
                // here because the sentence is true, not because a fixture is red.
                let mut vq = [query[0], query[1], query[2]];
                vq.sort_unstable();
                if side_of(jd, &[], NodeId::three_planes(vq), q) == Some(0)
                    && material(vq)? != Some(false)
                {
                    return Ok(None);
                }
                continue;
            }
            let mut x = [a, b, q];
            x.sort_unstable();
            let Some(in_g) = material(x)? else {
                return Ok(None);
            };
            let fwd = order_along(jd, a, b, c, q);
            if fwd == 0 {
                // The crossing is the ray origin itself (query on plane `q`). Inside `q`'s
                // material ⇒ query on the component's surface ⇒ undecidable → abandon.
                if in_g {
                    return Ok(None);
                }
                continue;
            }
            if fwd == 1 && in_g {
                count += 1; // one side of the line — parity is the same on either half
            }
        }
        Ok(Some(count % 2 == 1))
    };

    for i in 0..vplanes.len() {
        for j in (i + 1)..vplanes.len() {
            let (a, b) = (vplanes[i], vplanes[j]);
            // Locator `c`: a plane of the query off the line `a ∩ b`, so `{a,b,c}` is the query.
            let Some(&c) = vplanes
                .iter()
                .find(|&&x| x != a && x != b && jd.plane_pair_dir_sign(a, b, x) != 0)
            else {
                continue;
            };
            if let Some(inside) = attempt(a, b, c)? {
                return Ok(Some(inside));
            }
        }
    }
    // Every plane pair of this query was blocked: the node abstains. The inner `attempt`
    // already speaks this language per pair (`Ok(None)`); the boundary now keeps it instead of
    // dressing the abstention up as an error for the caller to catch and swallow.
    Ok(None)
}
/// **A ring node's coordinate, in whichever world names it** — the key
/// [`loop_winding`]'s lexicographic scan orders by.
///
/// ★★★ **The two arms are not two widths of one thing, they are two *kinds*.** A three-plane node
/// is the rational meet of three planes; a branch node's coordinate is `a + b√c` and no rational
/// vessel holds it. That is why this is an enum and not a `[Rat; 3]` with a decline: the second
/// arm is not a precision failure to be lifted, it is a different number.
enum CoordKey {
    /// The three classes, handed to `Judge::cmp_coord` — which keeps its toleranced ladder and its
    /// escalation, so the existing population's answers are bit-identical to before.
    Three([usize; 3]),
    /// The point a plane pair cuts out of a cylinder: `base + s·dir` with `s = a + b√c`.
    /// Boxed: this arm is an order of magnitude wider than a name, and a ring of names is the
    /// common case.
    Branch(Box<(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal)>),
}

/// Build one ring node's key.
///
/// ★ **The cylinder comes from the ring's own arcs, not from a class-table lookup.** A crossing
/// node has *both* an arc and a segment on it, so "the edge at index `i`" is the wrong place to
/// ask — the segment there names no cylinder. The name says which cylinder (`NodeId::Branch`
/// carries the class), and the arc that rides it is somewhere in this ring by construction: a
/// branch node is an arc endpoint.
fn coord_key(
    jd: &Judge<'_, WorkingPlane>,
    ring: &[RingEdge],
    i: usize,
) -> Result<CoordKey, BoolError> {
    let node = ring[i].node;
    let NodeId::Branch { cyl, .. } = node else {
        // A `match` and not a fallback: a third variant must light this up rather than fall in
        // here (`let`-`else` is what hid a new variant once already).
        return match node {
            NodeId::ThreePlane(t) => Ok(CoordKey::Three(t)),
            NodeId::Branch { .. } => unreachable!("the let-else above took every branch node"),
        };
    };
    let def = ring
        .iter()
        .find_map(|e| match &e.carrier {
            Carrier::Arc(a) if a.cyl == cyl => Some(&a.def),
            Carrier::Ruling(r) if r.cyl == cyl => Some(&r.def),
            Carrier::Plane { .. } | Carrier::Arc(_) | Carrier::Ruling(_) => None,
        })
        .ok_or_else(|| reject(RejectReason::BranchVertexUnnamed))?;
    let (line, s) =
        branch_meet(jd, cyl, def, node).ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
    Ok(CoordKey::Branch(Box::new((line, s))))
}

/// Order two ring nodes along one world axis — `+1` when `a`'s coordinate is the larger.
///
/// The mixed pair is `nacre_scalar::quad::cmp_coord_meet_branch`, which is exact and **total**:
/// both coordinates lift to one first-storey sign. It wants the three-plane side as a `MeetPoint`,
/// and this crate builds one directly — `MeetPoint::Narrow` is a public variant, so no door has to
/// be opened in `nacre-scalar` for it.
fn cmp_key(
    jd: &Judge<'_, WorkingPlane>,
    a: &CoordKey,
    b: &CoordKey,
    axis: usize,
) -> Result<i8, BoolError> {
    use nacre_scalar::{Orient, quad};
    let sign = |o: Orient| match o {
        Orient::Positive => 1i8,
        Orient::Negative => -1,
        Orient::Zero => 0,
    };
    // ★ Through [`node_coords_rat`], not a second copy of the three-plane solve — the rational
    // meet is stated once. The round-trip through the name is free: the solve is symmetric in its
    // three planes, so canonical order changes nothing.
    let meet = |t: [usize; 3]| {
        node_coords_rat(jd, NodeId::three_planes(t))
            .map(nacre_scalar::MeetPoint::Narrow)
            .ok_or_else(|| reject(RejectReason::WitnessNotRational))
    };
    Ok(match (a, b) {
        (CoordKey::Three(x), CoordKey::Three(y)) => jd.cmp_coord(*x, *y, axis),
        (CoordKey::Three(x), CoordKey::Branch(b)) => {
            sign(quad::cmp_coord_meet_branch(&meet(*x)?, &b.0, &b.1, axis))
        }
        (CoordKey::Branch(b), CoordKey::Three(y)) => {
            -sign(quad::cmp_coord_meet_branch(&meet(*y)?, &b.0, &b.1, axis))
        }
        (CoordKey::Branch(a), CoordKey::Branch(b)) => {
            sign(quad::cmp_coord_branch((&a.0, &a.1), (&b.0, &b.1), axis))
        }
    })
}

/// **The branch nodes lying strictly between two ring nodes, in ring order** — the exact half of
/// the split-twin subdivision (`boolean::name_result_vertices`' opening pass). The caller has
/// already matched the candidates' plane pair to the edge's `{own, wall}`, so by name every
/// candidate lies on the edge's own carrier line and single-axis order *is* order along it.
///
/// ★ **The axis is "wherever the endpoints differ", not the line's direction.** Two distinct
/// points of one line differ on some axis, the line is strictly monotone on that axis, and
/// betweenness is direction-blind — so no direction vector is read at all, and the choice is
/// deterministic (first differing axis). This is [`cmp_key`]'s vocabulary end to end; nothing new
/// is exact here.
///
/// `None` when an order cannot be formed (a coordinate outside the rational vessel, an
/// escalation, a class with no coefficients). The caller leaves such an edge **unsplit**, which
/// is today's behaviour exactly — the far-plane road starves there and the walls-fallback net
/// answers — so the conservative arm degrades to the state this pass was built to improve, never
/// to something new.
pub(crate) fn branch_between(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    a: NodeId,
    b: NodeId,
    candidates: &[NodeId],
) -> Option<Vec<NodeId>> {
    let key = |n: NodeId| -> Option<CoordKey> {
        match n {
            NodeId::ThreePlane(t) => Some(CoordKey::Three(t)),
            NodeId::Branch { cyl, .. } => {
                let (line, s) = branch_meet(jd, cyl, &cyls.get(cyl)?.def, n)?;
                Some(CoordKey::Branch(Box::new((line, s))))
            }
        }
    };
    let (ka, kb) = (key(a)?, key(b)?);
    let axis = (0..3).find(|&ax| matches!(cmp_key(jd, &ka, &kb, ax), Ok(s) if s != 0))?;
    // `+1` = "the first is larger" (`cmp_key`), so `dir` is the a→b slope's sign on this axis and
    // "strictly between" is both hops running the same way.
    let dir = cmp_key(jd, &ka, &kb, axis).ok()?;
    let mut mid: Vec<(NodeId, CoordKey)> = Vec::new();
    for &c in candidates {
        if c == a || c == b {
            continue;
        }
        let kc = key(c)?;
        if cmp_key(jd, &ka, &kc, axis).ok()? == dir && cmp_key(jd, &kc, &kb, axis).ok()? == dir {
            mid.push((c, kc));
        }
    }
    // Ring order: nearer to `a` first — along `dir`, the larger-toward-`a` side leads. A pair
    // this cannot strictly order (an escalation, or two distinct nodes at one coordinate — which
    // on a shared line would be two names for one point) makes the whole edge unsplittable.
    let mut sortable = true;
    mid.sort_by(|(_, x), (_, y)| match cmp_key(jd, x, y, axis) {
        Ok(s) if s == dir => std::cmp::Ordering::Less,
        Ok(s) if s == -dir => std::cmp::Ordering::Greater,
        _ => {
            sortable = false;
            std::cmp::Ordering::Equal
        }
    });
    if !sortable {
        return None;
    }
    Some(mid.into_iter().map(|(n, _)| n).collect())
}

/// An ordered ring's winding about the face's outward normal: `-1` clockwise — the material
/// is *outside* the ring, so it bounds a hole — and `+1` counter-clockwise, an island.
///
/// The turn at a convex-hull vertex is the winding, and the lexicographically smallest node
/// is one: it is an extreme point of the node set, which is planar, so it is a vertex of the
/// ring's hull. Finding it is the **only** thing here that needs two implicit points in one
/// decision, and [`three_plane_cmp_coord`](nacre_geom::intersect::three_plane_cmp_coord) is that predicate.
///
/// A shortcut dies here, and is recorded so it is not walked twice: a *supporting edge* —
/// one whose plane `Q_j` has every other node on one side — would give a hull vertex from
/// the one-implicit `three_plane_orient3d` alone. But a simple polygon need not have an edge
/// on its hull (fold each side of a pentagon slightly inward), so no such edge is guaranteed.
/// A hull *vertex* always exists.
///
/// ★★ **Where arcs touched this** (M6-2b, done): the turn is read at **one** node, so a curved edge
/// needs no angle sum — only its tangent's direction at that node. The other two sites that read
/// the direction's *representation* each got their own answer: the walk-back below asks
/// [`continuation`], which has a curved arm ("are the tangents parallel" is "same circle, same
/// travel"), and the lexicographic minimum above runs on [`CoordKey`], which holds a branch node's
/// `a + b√c` coordinate beside a name and compares across the two through the quad tower.
///
/// A ring may be *non-simple* — visiting one node twice — and still be a legitimate face: the
/// unbounded contour of two cells that meet at a single point pinches through that point, tracing
/// a figure-8. The winding is read from the turn at the lexicographically smallest node, and a
/// coincidence elsewhere in the ring does not affect that turn, so a repeated node is not by itself
/// an error. (The old check rejected on the first coincidence with the running minimum, which made
/// the verdict depend on the ring's arbitrary start index — one operand order rejected a pinch the
/// other accepted.) Only a pinch *at* the extreme node itself leaves the turn ambiguous; that stays
/// a `LOOP_ORIENT_MISMATCH`, decided by exact equality rather than by a tolerance.
pub(crate) fn loop_winding(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
) -> Result<i8, BoolError> {
    // ★ **The floor reads the carriers.** Two straight edges between two points are one edge traced
    // twice; two arcs are a lens and an arc with its chord is a circular segment. See the same rule
    // at the walk's orbit length (`arrangement::walk_cells`).
    if ring.len()
        < if ring.iter().any(|e| matches!(e.carrier, Carrier::Arc(_))) {
            2
        } else {
            3
        }
    {
        return Err(reject(RejectReason::DegenerateRing));
    }
    // ★ **The comparator, in one place, reading the identity directly.** This is the second of the
    // two sites that deliberately do not go through [`three_plane_name`]: `Judge::cmp_coord` speaks
    // three plane indices (it lives in `nacre-cip`, below this crate, so the name cannot travel
    // there), and a branch point's coordinate is `a + b√c` with its own total comparators
    // (`nacre_scalar::quad::cmp_coord_meet_branch` and `cmp_coord_branch`). **This dispatch is
    // where that decision belongs** — the door's single answer is the wrong one here, and
    // `Judge::cmp_coord`'s four-rung ladder speaks neither `MeetLine` nor `QuadVal`.
    //
    // ★★ **The keys are materialized, and that is a change of shape, not just of type.** The old
    // spelling read `[usize; 3]` out of the name on every comparison — free. A branch key is a
    // *solve* ([`branch_meet`] re-derives the point from the name), and the scan below asks for
    // each node's key on the order of six times, so re-solving per read would multiply the exact
    // work by that. One pass, one `Vec`.
    let keys: Vec<CoordKey> = (0..ring.len())
        .map(|i| coord_key(jd, ring, i))
        .collect::<Result<_, BoolError>>()?;
    let key = |i: usize| &keys[i];
    // Lexicographically smallest node — a hull vertex, hence a valid turn site. A coincidence with
    // the running minimum just means "not strictly smaller", so keep it; do not reject.
    let mut lo = 0usize;
    for i in 1..ring.len() {
        let mut order = 0i8;
        for axis in 0..3 {
            order = cmp_key(jd, key(i), key(lo), axis)?;
            if order != 0 {
                break;
            }
        }
        if order == -1 {
            lo = i;
        }
    }
    // ★★★ **The scan above is only a minimum if the relation is an order, and that is an
    // assumption about the predicates, not about this loop.** It composes per-axis comparisons
    // lexicographically; a per-axis answer that is not a fact about the geometry — one plane
    // described two ways, say — makes the composition intransitive, and then a forward scan can
    // stop at a node with something smaller behind it. That node is not extreme, the turn read
    // there is not the winding, and the wrong sign comes back **confident**: the engine noticed
    // only two layers later, as "no outer contour", and named the symptom.
    //
    // ★ This is the postcondition the algorithm actually needs — cheaper than asking whether the
    // relation is transitive (`O(n)` against `O(n³)`, and rings here reach 95 nodes) and closer to
    // the point. **It holds however the predicates behave; it is the net under them.**
    // ★ A declining comparison makes no claim, so it cannot witness a violation either: `0` is
    // "says nothing" here, not "equal".
    let lex = |i: usize, j: usize| -> i8 {
        (0..3)
            .map(|axis| cmp_key(jd, key(i), key(j), axis).unwrap_or(0))
            .find(|&c| c != 0)
            .unwrap_or(0)
    };
    debug_assert!(
        !(0..ring.len()).any(|i| i != lo && lex(i, lo) == -1),
        "the lexicographic scan did not find a minimum — the comparison is not an order here"
    );
    // The turn is read at `lo`; if that exact point recurs the corner is a pinch and its turn is
    // ambiguous — honest-reject rather than guess.
    // ★ Stops at the first pinch, as the `any` it replaced did. Running on would ask comparisons
    // the old spelling never made, and one of those could *decline* — turning a `CoincidentNodes`
    // that was already decided into a width reject.
    for i in 0..ring.len() {
        let mut same = i != lo;
        for axis in 0..3 {
            same = same && cmp_key(jd, key(i), key(lo), axis)? == 0;
        }
        if same {
            return Err(reject(RejectReason::CoincidentNodes));
        }
    }
    // **A ring node need not be a corner.** The arrangement names a point wherever another feature
    // crosses an edge, and `loop_triples` keeps such a vertex even when the loop runs straight
    // through it — so one edge of the polygon can arrive as several collinear ring edges. Reading
    // the turn at `lo` against its immediate predecessor then asks about two halves of one
    // straight edge, which has no turn to give.
    //
    // The turn to read is the one between the directions the loop **actually** arrives and leaves
    // on: walk back past the edges the loop runs straight through. `lo` stays a hull vertex — the
    // stretch lies on one **line** through it, so the region is still on one side of that line.
    // ★ That argument is the straight one, and the guard below is where it stops: two arcs of one
    // circle are also "straight through", and a *circle* through `lo` does not put the region on
    // one side of anything.
    //
    // **Only while the stretch keeps going the same way.** A collinear edge traversed the *other*
    // way means the ring doubles back along the line it came in on — an antenna, whose tip has no
    // turn and whose neighbours' turn belongs to a different vertex. Skipping past that would
    // read a turn from somewhere else and call it this vertex's: a wrong winding, silently.
    //
    // ★★★ **The comparison is between neighbours, at the node they share** — it used to be between
    // the candidate and `ring[lo]`, which is the same answer for straight edges (parallel and
    // same-sense are both transitive along a chain of shared points) and **meaningless** the moment
    // an edge is curved: a far arc's tangent is not `lo`'s tangent, so comparing them asks about
    // two different places. An earlier note here worried that neighbour-only would *weaken* "the
    // whole stretch runs one way"; transitivity is why it does not.
    //
    // ★★★★ **And the value carried out is the one already read at a shared node** — never an
    // edge's direction at its own far start. For a straight edge the two are the same value, which
    // is why the older spelling stood; on a **diameter** chord they are exactly opposite, and that
    // is what a boss straddling a plate edge measured (both half-disks are `+1`; the half whose arc
    // *arrives* came back `−1`). Where the stretch is curved the equality that licenses stepping at
    // all fails, and the walk refuses by name rather than reading a winding from the wrong place.
    let n = ring.len();
    let leaving = dir_at(jd, cyls, p, &ring[lo], ring[lo].node)?;
    let mut back = (lo + n - 1) % n;
    // ★ The direction the loop arrives on, carried out of the walk — the edge that ends it is the
    // one the turn is read against, and its direction is already in hand.
    let arriving = loop {
        let ahead = (back + 1) % n;
        let shared = ring[ahead].node;
        let earlier = dir_at(jd, cyls, p, &ring[back], shared)?;
        let later = dir_at(jd, cyls, p, &ring[ahead], shared)?;
        match continuation(jd, p, &earlier, &later) {
            // ★★★★ **`earlier`, and not `ring[back]`'s direction at its own start.** The two are
            // the same value for a straight edge — a line's tangent does not change along it — and
            // that equality is what let the older spelling stand. On an arc they differ by the
            // whole turn of the arc, and on a *diameter* chord they are exactly opposite: measured,
            // the two half-disks of a boss straddling a plate edge came back `+1` and `−1` where
            // both are `+1`, because the half whose arc **arrives** read its tangent at the far
            // end. This one is read at the node the loop actually passes through.
            Continuation::Turns => break earlier,
            Continuation::DoublesBack => return Err(reject(RejectReason::StraightAngle)),
            // ★★★ **The step is licensed by the stretch being a *line*.** What the walk carries out
            // is `ring[back]`'s direction at its own start, and that equals `lo`'s arriving
            // direction only because a line's tangent is the same everywhere on it. Two arcs of one
            // circle are tangent-continuous, so `continuation` answers `Straight` for them too —
            // and stepping there would read the winding from a different point of the ring.
            // ★★★★★ **At `lo` itself a smooth boundary still states a winding — by curvature.**
            // The step past a straight run is licensed by the stretch being a *line*, and two arcs
            // of one circle are tangent-continuous without being one; stepping there would read the
            // winding from a different point of the ring. But at `lo` there is nothing to step
            // past: the loop **is** smooth at the extreme node, and a smooth extremum's winding is
            // the arc's own rotation. See [`smooth_extremum_winding`].
            // Not curved here: an ordinary straight run through `lo`, walked back as before.
            Continuation::Straight
                if ahead == lo
                    && let (EdgeDir::Arc(e), EdgeDir::Arc(l)) = (&earlier, &later) =>
            {
                return smooth_extremum_winding(jd, p, e, l);
            }
            Continuation::Straight
                if matches!(earlier, EdgeDir::Arc(_)) || matches!(later, EdgeDir::Arc(_)) =>
            {
                return Err(reject(RejectReason::CurvedStraightRun));
            }
            Continuation::Straight => {}
        }
        back = (back + n - 1) % n;
        if back == lo {
            // Every edge of the ring lies on one line: it bounds nothing. ★ This is the *loop's*
            // termination, not one of `continuation`'s answers — it is about having walked the
            // whole ring, not about what any one edge does.
            return Err(reject(RejectReason::DegenerateRing));
        }
    };
    turn_between(jd, p, &arriving, &leaving)
}

/// **The winding of a ring that runs *smooth* through its extreme node** — curvature, not a turn.
///
/// ★★★★★ **This is not [`turn`]'s question, and that is why it is not [`turn`]'s arm.** `turn`
/// answers "how much does the direction rotate at this node", and for two arcs of one circle the
/// honest answer is `0`: they are tangent-continuous, nothing rotates *at* the node. What
/// [`loop_winding`] needs there is a different fact — which way the boundary **curves** — and it is
/// available only because the node is the ring's lexicographic minimum.
///
/// **Why the minimum makes it answerable.** `lo` is a hull vertex: the whole ring lies on one side
/// of a supporting line through it. The boundary there is an arc, so the arc curves off that line
/// into the side the ring is on — the region is locally convex at `lo`, and a locally convex point's
/// turn carries the ring's orientation. For an arc that "turn" is spread along the arc rather than
/// concentrated at a vertex, but its **sign** is the arc's own rotation, which is exactly the
/// winding.
///
/// **The sign, in three factors.** Travel rotates about `s·m`, `s = +1` when `ccw`. The turn's
/// reference is the face's **outward** normal, `n_out = frame_sign · n_P`, and `n_P · m > 0` is
/// `axis_up`. So
///
/// ```text
///   (s·m) · n_out = s · frame_sign · (n_P · m)
///     ⇒  winding = ccw · axis_up · frame_sign
/// ```
///
/// — no coordinate, no predicate, three signs the directions already carry.
///
/// The two arcs come from [`Continuation::Straight`], which for arcs means *one cylinder and the
/// same travel sense*, so both agree on every factor; they are re-checked here rather than assumed,
/// because this function's answer is a sign and a wrong one is silent.
///
/// ☑ **Which factors the population locks.** Negating the product, dropping `ccw`, and dropping
/// `axis_up` each move the chained fixtures' wall (to `NonManifoldResultEdge`,
/// [`RejectReason::RingOrientation`] and [`RejectReason::RulingBoundNotYet`] respectively), so the
/// lock names them. Dropping `frame_sign` changes nothing: it is `+1` on every class that reaches
/// this rule today, which is the same shape [`turn`]'s own note records for its factor — the
/// difference being that `turn`'s corpus does reach `Reversed` faces and this rule's does not yet.
fn smooth_extremum_winding(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    earlier: &ArcDir,
    later: &ArcDir,
) -> Result<i8, BoolError> {
    if earlier.cyl != later.cyl || earlier.ccw != later.ccw || earlier.axis_up != later.axis_up {
        return Err(reject(RejectReason::CurvedStraightRun));
    }
    let sign = |b: bool| if b { 1i8 } else { -1 };
    Ok(sign(later.ccw) * sign(later.axis_up) * jd.planes[p].frame_sign)
}

/// `sign((n_P × n_Q) · N_R)`, where `N_R` is the right-hand normal of `R.tri`.
///
/// [`plane_pair_dir_sign`](nacre_geom::intersect::plane_pair_dir_sign) gives the sign against `R`'s *stored* normal, exactly.
/// That normal is parallel to `N_R` but may oppose it on a `Reversed` face, so we
/// correct with their dot — two parallel unit vectors, `|·| ≈ 1`, nowhere near the
/// sign boundary.
///
/// The correction *is* the face's stated flag — since the stored-orientation
/// cutover, `frame_sign` is `Forward`/`Reversed` as a sign, and
/// "`Reversed` ⇔ `n_out = −plane.normal()`" holds by construction rather than by
/// hope. What keeps it honest is the winding: `collect_planes` debug_asserts the
/// witness triangle against `n_out`, and `validate` pins the loop itself as
/// `FaceMisoriented`.
pub(crate) fn dir_sign(jd: &Judge<'_, WorkingPlane>, p: usize, q: usize, r: usize) -> i8 {
    let planes = jd.planes;
    jd.plane_pair_dir_sign(p, q, r) * planes[r].frame_sign
}

// ---------------------------------------------------------------------------
// The **rational road** (M6-2a C4a): the same containment questions as above, asked about a
// point that has coordinates instead of a name.
//
// Everything else in this module names a point by three planes, which is what makes it exact.
// A cylinder's band has no such name to offer — its witness is a rational point on the axis —
// so the uniform-slab theorem needs a road that starts from coordinates and stays exact anyway.
// It does: a rational point, a rational direction, plane classes with narrow rational
// descriptions, and `point_in_ring_2d_rat` for the in-face parity. No `f64` decides anything
// here either.
// ---------------------------------------------------------------------------

/// A plane class's exact description, or `None` where it has none to give (a rotated class'
/// realized coefficients are not its truth, so they are refused rather than read).
/// **The carrier of a ring edge that rides a cylinder** — an arc or a ruling, filled or refused.
///
/// ★★ **The direction bits do not come from `edge_at`.** That convention ("even half-edge = CCW /
/// up") is the *arrangement's* edge indexing, and an operand's face loop has no such index. What is
/// true here is the model's own:
///
/// * an **arc** has a stated convention — `derive_edge_curve`'s (Plane, Cylinder) arm: "on a circle
///   carrier the vertex *order* says which arc; `[A, B]` is A to B **counter-clockwise about the
///   axis**". So walking the edge `forward` is walking it CCW.
/// * a **ruling** has none — the same arm says a plane parallel to the axis meets the lateral along
///   rulings and "the endpoints decide". `MergedRuling::end` ascending the axis is the
///   *arrangement's* convention, so `up` is **derived** here from the two endpoints' axial
///   coordinate rather than read off a rule.
///
/// `side` is [`crate::arrangement::ruling_side`]'s one spelling, and it needs a point on the ruling
/// *exactly* — which is why the branch name comes in: [`branch_meet`] realizes it as the `(line, s)`
/// that function takes. A ruling whose end is not a branch point (a seam end) has no such point and
/// is refused rather than guessed.
fn curved_wall(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    he: &nacre_topo::HalfEdge,
    cyl: usize,
    near: usize,
    end: NodeId,
) -> Result<crate::boolean::Wall, BoolError> {
    let curved = || reject(RejectReason::CurvedOperandBoundary);
    match model.edge_curve(he.edge) {
        nacre_geom::Curve::Circle(_) => Ok(crate::boolean::Wall::Arc {
            cyl,
            ccw: he.forward,
        }),
        nacre_geom::Curve::Line(_) => {
            let def = cyls.get(cyl).ok_or_else(curved)?.def.clone();
            let at = branch_meet(jd, cyl, &def, end).ok_or_else(curved)?;
            let w = class_coeffs_rat(jd, near).ok_or_else(curved)?;
            let side =
                crate::arrangement::ruling_side(&w, &def, (&at.0, &at.1)).ok_or_else(curved)?;
            // Which way travel runs along the axis: the stored edge ascends when its second
            // endpoint does, and `forward` says whether this half-edge walks it that way.
            //
            // ★★★ **There is nothing to read, so it is derived — from the two ends' definitions,
            // never from their realizations.** A ruling edge's stored pair carries no order:
            // `edge_for` keys it `unordered(va, vb)` ("a ruling edge is straight, so the unordered
            // pair orders it"), unlike an arc, whose `[A, B]` *is* the CCW convention. So the two
            // ends have to be compared — and each end is a branch point of `near`, the cylinder,
            // and **one other plane**. `near` *holds* the ruling, so it is that other plane that
            // **cuts** it, and where it crosses the axis is a rational question
            // ([`crate::planes::axis_param_of_plane`]).
            //
            // ☑ **The two parameters cannot tie**: the gate admits a plane that meets this
            // cylinder only parallel to the axis or perpendicular to it (anything between is
            // `ObliqueCylinderCut`), and a parallel one cannot cut a ruling — so both cutting
            // planes are caps, distinct caps cross the axis at distinct parameters. The strict
            // `>` therefore restates the comparison it replaces exactly, rather than growing a
            // decline for a case that has none.
            let other_param = |v: Handle<Vertex>| -> Option<nacre_scalar::Rat> {
                let nacre_topo::VertexDef::Branch { planes, .. } = model.vertices.get(v).def else {
                    return None;
                };
                let mut cut = None;
                for &h in &planes {
                    let c = *model.world_plane_name(h)?.narrow()?;
                    // `near`'s own class, in whichever of the two spellings this vertex carries.
                    if plane_sense(&c, &w).is_some() {
                        continue;
                    }
                    if cut.replace(c).is_some() {
                        return None;
                    }
                }
                crate::planes::axis_param_of_plane(&cut?, &def)
            };
            let [v0, v1] = model.edges.get(he.edge).vertices;
            let (t0, t1) = (
                other_param(v0).ok_or_else(curved)?,
                other_param(v1).ok_or_else(curved)?,
            );
            let ascends = t1 > t0;
            Ok(crate::boolean::Wall::Ruling {
                cyl,
                side,
                up: ascends == he.forward,
            })
        }
    }
}

/// **Do these two exact 4-vectors describe one plane, and does `b`'s normal point the same way?**
///
/// `Some(+1)` same plane same sense, `Some(-1)` same plane opposite sense, `None` different planes
/// (or overflow). The comparison is cross-multiplication against a nonzero component — the idiom
/// [`nacre_scalar::quad::plane_plane_cylinder`]'s parallel arm already uses ("coincident iff the
/// full 4-vectors are proportional"), spelled once here because a *second* caller now needs it.
fn plane_sense(a: &[nacre_scalar::Rat; 4], b: &[nacre_scalar::Rat; 4]) -> Option<i8> {
    let zero = nacre_scalar::Rat::from_int(0);
    let i = (0..3).find(|&k| a[k] != zero)?;
    if b[i] == zero {
        return None;
    }
    for j in 0..4 {
        if b[j].checked_mul(a[i])? != a[j].checked_mul(b[i])? {
            return None;
        }
    }
    Some(if (a[i] > zero) == (b[i] > zero) {
        1
    } else {
        -1
    })
}

/// **An operand vertex's own branch name, restated in this arrangement's class space.**
///
/// ★★★★★ **The second half of a correspondence the forward direction gets for free.** When a
/// boolean *mints* a [`nacre_topo::VertexDef::Branch`] it writes the two planes as **its own
/// classes' representative surfaces**, so restating class order as handle order is the only
/// correction it needs ([`nacre_topo::QuadRoot::canonical`], which `assemble` calls). Coming back
/// the other way the handles are **given**, and they may be a surface that merged into a class
/// under a different representative — and, because a class holds faces whose normals oppose, under
/// the **opposite sign**. `VertexDef::Branch`'s own doc says this correspondence "has to be
/// established a second time"; this is that time.
///
/// ★★ **Both corrections are the same rule.** `Lo`/`Hi` are the order along `ℓ = n₁ × n₂`, and
/// `plane_plane_cylinder` fixes the base by `{n₁·x = −d₁, n₂·x = −d₂, ℓ·x = 0}` — a condition
/// `−ℓ` satisfies identically. So negating **either** normal leaves the two planes, and the base,
/// exactly where they were and only reverses `ℓ`: the two roots trade places, which is what
/// `flipped` says. Swapping the pair reverses `ℓ` too (the derivation `canonical` already carries).
/// ⇒ **flip once per reversal, and an even number of reversals is no flip at all.** The swap is
/// [`NodeId::branch`]'s to count; the two signs are this function's.
/// ☑ The cylinder needs no correction: negating its axis direction does not move the surface, so
/// the two roots are the same two points in the same order.
pub(crate) fn branch_name_from_def(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    v: Handle<Vertex>,
    cyl: usize,
    candidates: [usize; 2],
) -> Option<NodeId> {
    let nacre_topo::VertexDef::Branch { planes, root, .. } = model.vertices.get(v).def else {
        return None;
    };
    // Which candidate class each stored handle *is*, and with which sense. The match decides the
    // correspondence and the sign in one comparison — asking them separately would be two chances
    // to disagree.
    let mut seen: [Option<(usize, i8)>; 2] = [None, None];
    for (i, &h) in planes.iter().enumerate() {
        let name = model.world_plane_name(h)?;
        let c = name.narrow()?;
        for &k in &candidates {
            let Some(sense) = plane_sense(c, &class_coeffs_rat(jd, k)?) else {
                continue;
            };
            if seen[i].is_some() {
                return None; // one handle answering to both classes is not a correspondence
            }
            seen[i] = Some((k, sense));
        }
    }
    let ((k0, s0), (k1, s1)) = (seen[0]?, seen[1]?);
    if k0 == k1 {
        return None;
    }
    let root = if s0 == s1 { root } else { root.flipped() };
    Some(NodeId::branch(k0, k1, cyl, root))
}

pub(crate) fn class_coeffs_rat(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
) -> Option<[nacre_scalar::Rat; 4]> {
    jd.planes[c].world_rat
}

/// A node's exact coordinates: the rational meet of its three classes' descriptions. `None` when
/// any class lacks a narrow rational description — the caller declines rather than guessing.
///
/// ★ **This reads the identity directly rather than going through [`three_plane_name`]**, because
/// it is one of the two places whose answer for a second variant is *its own*: a branch point's
/// coordinate is `a + b√c`, not a rational meet, so this `match` is where that decision belongs —
/// not behind a door whose one answer is "no three-plane name, nothing to give".
///
/// ★★ **A branch node's `None` is a type fact, not a width decline** — the vessel is rational and
/// the coordinate is not, so no amount of precision reaches it. Resist answering `Some` for a
/// tangency because *its* coordinate happens to be rational: a function right for one root and not
/// the other is the "sometimes right" trap.
///
/// ★★★ **And the live caller names that `None` wrongly for this cause — reachably, now.**
/// `arrangement::split_circles` asks it of both ends of every segment and turns `None` into
/// `WitnessNotRational`, whose sentence is "a wider rational would lift this" — **false** for a
/// branch point, which has no rational coordinate at any width. This is where every chained-cylinder
/// operand stops today, so the wrong sentence is the one a user meets. The fix is structural rather
/// than a rename: that caller is being taught to ask for an **order** instead of a coordinate, and
/// then it will not ask this at all. (The predecessor this paragraph used to name,
/// `circles_meet_no_segment`, no longer exists.)
pub(crate) fn node_coords_rat(
    jd: &Judge<'_, WorkingPlane>,
    n: NodeId,
) -> Option<[nacre_scalar::Rat; 3]> {
    match n {
        NodeId::Branch { .. } => None,
        NodeId::ThreePlane(t) => nacre_scalar::three_planes_rat([
            class_coeffs_rat(jd, t[0])?,
            class_coeffs_rat(jd, t[1])?,
            class_coeffs_rat(jd, t[2])?,
        ]),
    }
}

/// **A node's realized coordinate, whichever kind of name it is** — the `f64` sibling of
/// [`node_coords_rat`] and [`branch_point`], which each answer for one variant only.
///
/// ★★ **It lives here because the `match` does.** Reaching into a [`NodeId`] variant outside this
/// file (and `reuse`) is what [`three_plane_name`]'s gate forbids, and the first draft of this
/// function sat in `arrangement` and broke it — with the whole suite green, exactly as that gate's
/// doc predicts. The two roads it dispatches between are already both here, so this is where the
/// third question about the same name belongs.
///
/// ★ `cfg(test)` only while the audit is its one consumer. The seam table's `VertexDef` minting
/// asks for the same pair of roads and will want it in production.
#[cfg(test)]
pub(crate) fn node_point_f64(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
) -> Option<[f64; 3]> {
    match n {
        NodeId::ThreePlane(t) => nacre_geom::intersect::three_planes(
            &jd.planes[t[0]].plane,
            &jd.planes[t[1]].plane,
            &jd.planes[t[2]].plane,
        )
        .map(|p| p.as_array()),
        NodeId::Branch { cyl, .. } => branch_point(jd, cyl, &cyls[cyl].def, n),
    }
}

/// **A branch node's realized coordinate, derived from the name.**
///
/// The sibling of [`node_coords_rat`] for the other variant, and the two return types *are* the
/// distinction: a branch coordinate is `a + b√c`, so the rational vessel next door cannot hold it
/// and only the cache can.
///
/// ★★ **It re-solves from the name rather than taking the producer's `(line, s)`** — the truth is
/// the definition and the coordinate is its cache. That is also what makes canonicalization
/// load-bearing: a name whose root failed to follow its pair through the sort designates the
/// *other* crossing, and the point moves where a test can see it.
///
/// ★ **Strict about a tangency.** `Double` is answered only by `Tangent` and `Lo`/`Hi` only by
/// `Pair`, and the mismatches are `None` rather than a nearest guess — so a second name for a
/// tangency's one point is unrepresentable here, not merely discouraged.
///
/// `def` must be the cylinder class `n` names; the caller holds the class table's row and this has
/// no way to look one up, so `cyl` comes with it and the two are checked against the name rather
/// than promised — a mismatched pair would otherwise realize a real point of the *wrong* cylinder.
/// `None` is a three-plane node, a class with no rational description, a root the meet does not
/// have, or checked-`Rat` overflow.
pub(crate) fn branch_point(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<[f64; 3]> {
    let (line, s) = branch_meet(jd, cyl, def, n)?;
    Some(nacre_scalar::quad::branch_point_f64(&line, &s))
}

/// **The exact half of [`branch_point`]** — the `(line, s)` the name designates, before it is
/// realized.
///
/// ★ «진실은 정의, 좌표는 캐시»: the pair *is* the point and the `[f64; 3]` beside it is its
/// realization, so the two are one function split in the middle rather than two solves. Every
/// exact question about a branch point — its order along the line, its side of a plane, its θ about
/// the seam — takes this and never the realization.
pub(crate) fn branch_meet(
    jd: &Judge<'_, WorkingPlane>,
    cyl: usize,
    def: &nacre_topo::CylinderDef,
    n: NodeId,
) -> Option<(nacre_scalar::quad::MeetLine, nacre_scalar::quad::QuadVal)> {
    use nacre_scalar::quad::{CylinderMeet, QuadVal};
    use nacre_topo::QuadRoot;
    let (planes, root) = match n {
        NodeId::Branch {
            planes,
            cyl: named,
            root,
        } => {
            debug_assert_eq!(
                named, cyl,
                "the def handed in is not the cylinder the name says"
            );
            (planes, root)
        }
        NodeId::ThreePlane(_) => return None,
    };
    let (p1, p2) = (
        class_coeffs_rat(jd, planes[0])?,
        class_coeffs_rat(jd, planes[1])?,
    );
    let (o, m, r) = (def.origin(), def.dir(), def.radius());
    let (line, s) = match (
        nacre_scalar::quad::plane_plane_cylinder(&p1, &p2, &o, &m, r)?,
        root,
    ) {
        (CylinderMeet::Pair { line, s }, QuadRoot::Lo) => (line, s[0]),
        (CylinderMeet::Pair { line, s }, QuadRoot::Hi) => (line, s[1]),
        (CylinderMeet::Tangent { line, s }, QuadRoot::Double) => (line, QuadVal::from_rat(s)),
        _ => return None,
    };
    Some((line, s))
}

pub(crate) fn dot3_rat(
    x: &[nacre_scalar::Rat; 3],
    y: &[nacre_scalar::Rat; 3],
) -> Option<nacre_scalar::Rat> {
    x[0].checked_mul(y[0])?
        .checked_add(x[1].checked_mul(y[1])?)?
        .checked_add(x[2].checked_mul(y[2])?)
}

pub(crate) fn cross3_rat(
    x: &[nacre_scalar::Rat; 3],
    y: &[nacre_scalar::Rat; 3],
) -> Option<[nacre_scalar::Rat; 3]> {
    Some([
        x[1].checked_mul(y[2])?
            .checked_sub(x[2].checked_mul(y[1])?)?,
        x[2].checked_mul(y[0])?
            .checked_sub(x[0].checked_mul(y[2])?)?,
        x[0].checked_mul(y[1])?
            .checked_sub(x[1].checked_mul(y[0])?)?,
    ])
}

/// **A rational 2D chart of a plane**, for running a parity test in it.
///
/// `e₁ = ê_k × n` for the first basis axis giving a nonzero cross, `e₂ = n × e₁`. The chart is
/// deliberately **not** orthonormal: crossing parity is invariant under any affine isomorphism of
/// the plane, and demanding unit vectors would need square roots that leave the rationals. One
/// copy of this rule, because a second spelling of it is how two consumers start disagreeing
/// about which side of a ring a point is on.
pub(crate) struct Chart2dRat {
    e1: [nacre_scalar::Rat; 3],
    e2: [nacre_scalar::Rat; 3],
}

impl Chart2dRat {
    /// The chart of the plane with normal `n`. `None` on a zero normal or on checked-`Rat`
    /// overflow.
    pub(crate) fn of_normal(n: &[nacre_scalar::Rat; 3]) -> Option<Self> {
        let zero = nacre_scalar::Rat::from_int(0);
        let basis = |k: usize| -> [nacre_scalar::Rat; 3] {
            let mut e = [zero; 3];
            e[k] = nacre_scalar::Rat::from_int(1);
            e
        };
        let e1 = (0..3)
            .filter_map(|k| cross3_rat(&basis(k), n))
            .find(|e| e.iter().any(|c| *c != zero))?;
        let e2 = cross3_rat(n, &e1)?;
        Some(Chart2dRat { e1, e2 })
    }

    /// A point's chart coordinates.
    pub(crate) fn project(&self, p: &[nacre_scalar::Rat; 3]) -> Option<[nacre_scalar::Rat; 2]> {
        Some([dot3_rat(p, &self.e1)?, dot3_rat(p, &self.e2)?])
    }

    /// The chart's two in-plane axes — the mixed-ring parity walks its ray along `e1` (one
    /// decision rule with this chart, not a second spelling of a basis).
    pub(crate) fn axes(&self) -> (&[nacre_scalar::Rat; 3], &[nacre_scalar::Rat; 3]) {
        (&self.e1, &self.e2)
    }

    /// A ring of nodes in chart coordinates.
    pub(crate) fn ring(
        &self,
        jd: &Judge<'_, WorkingPlane>,
        ring: &[NodeId],
    ) -> Option<Vec<[nacre_scalar::Rat; 2]>> {
        ring.iter()
            .map(|&n| self.project(&node_coords_rat(jd, n)?))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The canonical form is the sorted triple, and the constructor is the only way to get one —
    /// so the same three planes named in any order are **one** name.
    ///
    /// The expected value is written out rather than derived, because a test that re-derives it
    /// through the constructor would agree with the constructor however the constructor behaved.
    #[test]
    fn three_planes_names_one_vertex_however_it_is_spelled() {
        let canonical = NodeId::ThreePlane([2, 5, 9]);
        for spelling in [
            [2, 5, 9],
            [2, 9, 5],
            [5, 2, 9],
            [5, 9, 2],
            [9, 2, 5],
            [9, 5, 2],
        ] {
            assert_eq!(
                NodeId::three_planes(spelling),
                canonical,
                "{spelling:?} names the same vertex as [2, 5, 9]"
            );
        }
        assert_eq!(
            three_plane_name(canonical),
            Some([2, 5, 9]),
            "and the door agrees"
        );
    }

    /// ★ **`Ord` is the bare triple's lexicographic order** — the proposition the whole migration
    /// to this type stands on. Four rules read it and would answer differently if it moved:
    /// `Aliases::union_point`'s "smallest name wins", `merge_component`'s sorted ring starts,
    /// `reuse::canonical`'s ring rotation, and the reject witness's `Break::key`.
    ///
    /// Each expected sign is written by hand from the pair, not computed from either operand, so
    /// the test cannot agree with a wrong implementation by sharing its derivation.
    #[test]
    fn a_name_orders_like_the_triple_it_is() {
        use std::cmp::Ordering::{Equal, Greater, Less};
        for (a, b, want) in [
            ([0, 1, 2], [0, 1, 3], Less),    // last component decides
            ([0, 1, 3], [0, 1, 2], Greater), // …and antisymmetrically
            ([0, 2, 9], [0, 3, 4], Less),    // middle decides before last
            ([1, 0, 0], [0, 9, 9], Greater), // first decides before middle
            ([4, 4, 4], [4, 4, 4], Equal),
        ] {
            assert_eq!(
                NodeId::ThreePlane(a).cmp(&NodeId::ThreePlane(b)),
                want,
                "{a:?} vs {b:?}"
            );
        }
    }

    /// **A branch point named either way round is one name** — the [`NodeId::three_planes`]
    /// proposition for the second variant, where the canonicalization is bigger than a sort.
    ///
    /// The expected values are written from the rule's sentence, not by calling the constructor a
    /// second time. The rule itself is [`nacre_topo::QuadRoot::canonical`]'s, and its own locks
    /// (in `nacre-topo`) are what say the *geometry* agrees; this only says the name reads it.
    #[test]
    fn a_branch_is_one_name_however_the_pair_is_handed_in() {
        use nacre_topo::QuadRoot::{Double, Hi, Lo};
        for (root, other) in [(Lo, Hi), (Hi, Lo), (Double, Double)] {
            assert_eq!(
                NodeId::branch(9, 2, 7, root),
                NodeId::branch(2, 9, 7, other),
                "{root:?} handed in descending is {other:?} in stored order"
            );
        }
        // ★ And the two roots of one pair stay **two** names — a canonicalization that collapsed
        // them would make this pass by making everything equal.
        assert_ne!(NodeId::branch(2, 9, 7, Lo), NodeId::branch(2, 9, 7, Hi));
        // A different cylinder through the same two planes is a different point.
        assert_ne!(NodeId::branch(2, 9, 7, Lo), NodeId::branch(2, 9, 8, Lo));
    }

    /// ★ **The door tells the two variants apart, and `Ord` places them deterministically.**
    ///
    /// The order itself is a *choice* — derived `Ord` puts every `ThreePlane` before every
    /// `Branch`, so "the smallest name wins" gains a systematic lean toward three-plane names at
    /// the six rules that read it. Nothing can observe it yet (no ring holds a branch node), which
    /// is exactly why it is written down here rather than left to be discovered.
    #[test]
    fn the_two_variants_are_told_apart_and_ordered() {
        let three = NodeId::three_planes([9, 2, 5]);
        let branch = NodeId::branch(2, 9, 7, nacre_topo::QuadRoot::Lo);
        assert_eq!(three_plane_name(three), Some([2, 5, 9]));
        assert_eq!(
            three_plane_name(branch),
            None,
            "a branch point has no triple"
        );
        assert!(
            three < branch,
            "declaration order: ThreePlane before Branch"
        );
        // The probe helper drops exactly the branch node, and keeps the ring's order otherwise.
        assert_eq!(
            three_plane_probes([branch, three, branch]),
            vec![[2, 5, 9]],
            "a probe list may lose a member; that is its licence"
        );
    }
}
