use super::*;
/// A solid's edges, each with its endpoints and **the indices of its two faces**.
///
/// Named for what it holds: the pair is `planes`-table slots, i.e. *faces*, not plane classes —
/// `edge_faces` builds it from `surf_ix: HashMap<Handle<Face>, usize>`. It was `EdgePlanes`, and
/// that name is how a face index gets read as a plane one. The face→plane step is `plane_ix`, and
/// it happens in `loop_triples`, nowhere else.
///
/// ★ It also holds **every vertex's incident faces**, read off the same walk. That table
/// is what lets a ring vertex be named from *all* the planes through it rather than from the one
/// face loop asking — and that set takes no `O(V·P)` query: the
/// operand's own topology knows it, and the edge walk visits every (vertex, face) pair anyway.
/// One edge's incidence: its two bounding vertices and the two faces (table slots) using it.
pub(crate) type EdgeIncidence = ([Handle<Vertex>; 2], [usize; 2]);

pub(crate) struct EdgeFaces {
    edges: HashMap<Handle<Edge>, EdgeIncidence>,
    vertex_faces: HashMap<Handle<Vertex>, Vec<usize>>,
}

impl EdgeFaces {
    /// An edge's two bounding vertices and its two faces (table slots).
    pub(crate) fn get(&self, e: &Handle<Edge>) -> Option<&EdgeIncidence> {
        self.edges.get(e)
    }
    /// Every edge's `(bounds, faces)`, in no particular order — the audits' walk.
    #[cfg(test)]
    pub(crate) fn edges(&self) -> impl Iterator<Item = &EdgeIncidence> {
        self.edges.values()
    }
    /// The faces (table slots) with an edge at `vh`, sorted — empty for a vertex no edge of this
    /// solid touches.
    pub(crate) fn faces_at(&self, vh: Handle<Vertex>) -> &[usize] {
        self.vertex_faces.get(&vh).map(Vec::as_slice).unwrap_or(&[])
    }
}

/// Index a solid's edges by handle — [`edge_incidence`] keyed for lookup.
///
/// Hole-ring edges are in here too: `edge_incidence` walks every loop of every
/// face, and rejects an edge whose incidence is not a pair.
pub(crate) fn edge_faces(
    model: &Model,
    solid: Handle<Solid>,
    surf_ix: &HashMap<Handle<Face>, usize>,
) -> Result<EdgeFaces, BoolError> {
    let mut edges = HashMap::new();
    let mut vertex_faces: HashMap<Handle<Vertex>, Vec<usize>> = HashMap::new();
    for (eh, bounds, inc) in edge_incidence(model, solid, surf_ix)? {
        for &vh in &bounds {
            let at = vertex_faces.entry(vh).or_default();
            for &f in &inc {
                if !at.contains(&f) {
                    at.push(f);
                }
            }
        }
        edges.insert(eh, (bounds, inc));
    }
    for at in vertex_faces.values_mut() {
        at.sort_unstable();
    }
    Ok(EdgeFaces {
        edges,
        vertex_faces,
    })
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
    /// Must be a [`NodeKind::Pierce`]; the quadric road reads its cylinder and root.
    Pierce(NodeId),
}

impl PointOn {
    /// ★★★★ **A cylinder pin arrives with a pierce name — a producer's invariant, so this stops
    /// loudly rather than naming a refusal.** Every site that writes `EndPin::Cylinder` writes the
    /// `NodeId::Pierce` beside it, so the two disagreeing is a defect in this kernel rather than a
    /// property of the model, and [`RejectReason::PierceVertexUnnamed`] would be a false sentence
    /// for it: that one says the point *is* exactly named and only the path lacks a way to carry
    /// it. The rule is `crate::planes::ClassIx::plane`'s — *"a loud panic beats a silently wrong
    /// plane."*
    /// ☑ Measured unexercised over the whole suite and the ignored sweep before it was made loud.
    #[inline]
    fn of(pt: (NodeId, EndPin)) -> PointOn {
        match pt.1 {
            EndPin::Class(r) => PointOn::Class(r),
            EndPin::Cylinder => match pt.0.kind() {
                NodeKind::Pierce { .. } => PointOn::Pierce(pt.0),
                NodeKind::ThreePlane(_) => {
                    unreachable!(
                        "a cylinder pin was written beside a three-plane name: {:?}",
                        pt.0
                    )
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
/// — the shape of a real bug (passing the wrong `q`). Built here, `ds`
/// cannot disagree with the `(p, q)` it was built for. And the symmetry pays again: a segment
/// endpoint's `Located` is built once per segment, so a pierce end's [`pierce_meet`] is never
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
    Pierce {
        name: NodeId,
        meet: Box<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)>,
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
    Pierce {
        name: NodeId,
        meet: Box<(nacre_exact::quad::MeetLine, nacre_exact::quad::QuadVal)>,
    },
}

impl Located<'_> {
    /// The name this point is known by, given the line it sits on — derived for a plane pin.
    pub(crate) fn name(&self, p: usize, q: usize) -> NodeId {
        match self {
            Located::Class { pin, .. } => NodeId::three_planes(Canon3::three([p, q, *pin])),
            Located::Pierce { name, .. } => *name,
        }
    }
}

impl OnLine {
    /// The name this point is known by, given the line it sits on — derived for a plane pin.
    pub(crate) fn name(&self, p: usize, q: usize) -> NodeId {
        match self {
            OnLine::Class { pin, .. } => NodeId::three_planes(Canon3::three([p, q, *pin])),
            OnLine::Pierce { name, .. } => *name,
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
        (Located::Pierce { name: x, .. }, OnLine::Pierce { name: y, .. }) => x == y,
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
        PointOn::Pierce(name) => {
            // ★ A `match` and not a fallback: a third variant must light this up rather than
            // fall in with the three-plane one (`winding`'s `coord_key` states the rule).
            let cyl = match name.kind() {
                NodeKind::Pierce { cyl, .. } => cyl,
                NodeKind::ThreePlane(_) => return None,
            };
            let meet = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, name)?;
            Some(OnLine::Pierce {
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
        PointOn::Pierce(name) => {
            // ★ A `match` and not a fallback: a third variant must light this up rather than
            // fall in with the three-plane one (`winding`'s `coord_key` states the rule).
            let cyl = match name.kind() {
                NodeKind::Pierce { cyl, .. } => cyl,
                NodeKind::ThreePlane(_) => return None,
            };
            let meet = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, name)?;
            Some(Located::Pierce {
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
    use nacre_exact::Orient;
    use nacre_exact::quad::{cmp_coord_branch, cmp_coord_meet_branch};
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
    let rat = |n: NodeId| -> Option<nacre_exact::MeetPoint> {
        node_coords_rat(jd, n).map(nacre_exact::MeetPoint::Narrow)
    };
    let cmp = match (a, b) {
        // ☑ Unreachable — two plane-pinned points returned through the integer road above. Spelled
        // as a decline rather than a panic or a second rational comparison nobody would exercise.
        (Located::Class { .. }, OnLine::Class { .. }) => return None,
        (Located::Class { .. }, OnLine::Pierce { meet, .. }) => sign(cmp_coord_meet_branch(
            &rat(a.name(p, q))?,
            &meet.0,
            &meet.1,
            k,
        )),
        (Located::Pierce { meet, .. }, OnLine::Class { .. }) => -sign(cmp_coord_meet_branch(
            &rat(b.name(p, q))?,
            &meet.0,
            &meet.1,
            k,
        )),
        (Located::Pierce { meet: m1, .. }, OnLine::Pierce { meet: m2, .. }) => {
            sign(cmp_coord_branch((&m1.0, &m1.1), (&m2.0, &m2.1), k))
        }
    };
    Some(cmp * dsign)
}
