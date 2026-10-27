use super::*;
/// A solid's edges, each with its endpoints and **the indices of its two faces**.
///
/// Named for what it holds: the pair is `planes`-table slots, i.e. *faces*, not plane classes —
/// `edge_faces` builds it from `surf_ix: HashMap<Handle<Face>, usize>`. A name like `EdgePlanes`
/// is how a face index gets read as a plane one. The face→plane step is `plane_ix`, and
/// it happens in `loop_triples`, nowhere else.
///
/// It also holds **every vertex's incident faces**, read off the same walk. That table
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
/// `V_j ∈ R_j`, [`Judge::orient3d`] gives `sign((V_i − V_j)·N_j)` for `N_j` the outward normal of
/// `R_j`'s witness triangle; multiplying by `sign(d·N_j)` ([`dir_sign`]) recovers the order.
/// Both factors are judged on the truth, so the comparator is a true total order.
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

/// **The order of two points along `L = P ∩ Q`, whatever pins them** — [`order_on`] for an end
/// carried as a (name, pin) pair.
///
/// **Which direction is `order_along`'s, derived rather than matched.** [`Judge::orient3d`]
/// is *not* a four-plane determinant — it is the point `V = {p, q, i}` against the **triangle** `j`,
/// so it is symmetric in its first three. With `Vᵢ − Vⱼ = t·(n_p × n_q)` (both lie on `L`),
///
/// ```text
/// order_along = sign((Vᵢ−Vⱼ)·N_out(j)) · sign((n_p×n_q)·N_out(j)) = sign(t·(…)²) = sign(t)
/// ```
///
/// — so it orders along **`n_p × n_q` over the stored normals**, and `j`'s own orientation cancels
/// as a square. That is why the second road below reads its axis and sign off
/// [`stored_line_dir`] — that same product, exactly, over the classes' stored-oriented rational
/// coefficients. A cross product of the classes' `world_rat` normals is **not** it: those are
/// canonical and may oppose the stored ones per plane, which flips the direction (measured on a
/// `world_rat` cross: 80 of 160 comparisons against `order_along` came out reversed).
///
/// **One rule, two numbers.** A pair with a plane pin is ordered by that formula — the plane-pinned
/// point's side through the integer predicates, a pierce point's through one `a + b√c` sign
/// ([`PiercePoint::side`]); only two pierce points, neither with a pin, go through the `a + b√c` tower on
/// one coordinate axis. `None` is a missing description (a pin with no exact description, a
/// cylinder with no world statement), never a shape this cannot order.
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

/// **The order of two [`PointOn`]s along `L = P ∩ Q`** — the rule [`order_pinned`] states: two plane
/// pins by [`order_along`], any pair with a pierce point through the located forms
/// ([`order_located`]).
///
/// **Two plane pins take `order_along` here, before anything is located.** The located road gives
/// the same sign (its first arm is `order_along` with the Cramer parts cached), but building the
/// forms for a single comparison is not free: a variant where the ray caster and the chord test
/// asked their plane pairs that way — the ray even locating `v` once per ray — ran 2–4% slower on
/// the axis-aligned 80-fin fold, with and without one codegen unit.
#[inline]
pub(crate) fn order_on(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    q: usize,
    a: PointOn,
    b: PointOn,
) -> Option<i8> {
    if let (PointOn::Class(i), PointOn::Class(j)) = (a, b) {
        return Some(order_along(jd, p, q, i, j));
    }
    let (a, b) = (locate(jd, cyls, p, q, a)?, on_line(jd, cyls, p, q, b)?);
    order_located(jd, p, q, &a, &b)
}

/// **What a point on `L = P ∩ Q` is given as** — the input [`locate`] and [`on_line`] take.
///
/// **A plane pin brings no name, and that is the point.** `{p, q, pin}` *is* the name, so
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
    /// **A cylinder pin arrives with a pierce name — a producer's invariant, so this stops
    /// loudly rather than naming a refusal.** Every site that writes `EndPin::Cylinder` writes the
    /// `NodeId::Pierce` beside it, so the two disagreeing is a defect in this kernel rather than a
    /// property of the model, and [`RejectReason::PierceVertexUnnamed`] would be a false sentence
    /// for it: that one says the point *is* exactly named and only the path lacks a way to carry
    /// it. The rule is `crate::planes::ClassIx::plane`'s — *"a loud panic beats a silently wrong
    /// plane."*
    /// Measured unexercised over the whole suite and the ignored sweep.
    #[inline]
    pub(crate) fn of(pt: (NodeId, EndPin)) -> PointOn {
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
/// **This exists so the rule can have one implementation *and* keep its hoists.** The
/// arrangement's interval overlay is the boolean's hottest named phase (12–24% of it, 1.4–3.2M
/// containment questions per 60-fin fold), and its speed comes from building a point's Cramer
/// parts **once per wall pair** and its `dir_sign` **once per segment endpoint**, then asking many
/// questions of them. A rule that takes bare `(name, pin)` throws both away on every call. The
/// alternative — a fast inlined copy beside the general one — is the defect shape this codebase
/// keeps finding in itself, so the hoist becomes an *argument* instead of a second spelling.
///
/// **Both sides are `Located`, and the symmetry is load-bearing.** A draft took one side bare
/// plus a separate `b_ds: Option<i8>`, which is a value that has to match arguments it cannot see
/// — the shape of a real bug (passing the wrong `q`). Built here, `ds`
/// cannot disagree with the `(p, q)` it was built for. And the symmetry pays again: a segment
/// endpoint's `Located` is built once per segment, so a pierce end's [`pierce_meet`] is never
/// re-solved per (split point × segment).
pub(crate) enum OnLine {
    /// A point three planes name. `ds` is `dir_sign(p, q, pin)`.
    ///
    /// `ds` is read when this point is the **second** argument — the formula is
    /// `at_a.orient3d(b.pin) × ds_b`, which is [`order_along`] with `a`'s half hoisted.
    /// **No `name`.** A plane-pinned point on `L = P ∩ Q` *is* `{p, q, pin}`, so the name is
    /// derivable, and storing it would put a 32-byte `NodeId` in a `Vec` the containment loop reads
    /// millions of times. What a name would be for — identity against the other side — needs none:
    /// a plane-pinned point can never *equal* a cylinder-pinned one, so cross-variant identity is
    /// `false` by type, and a pierce point is ordered against this one by its side of `pin`.
    Class { pin: usize, ds: i8 },
    /// A point a cylinder pins: the meet line of its two planes and the root along it.
    ///
    /// **Boxed, and the profile is why** — the same reason [`EdgeDir::Arc`] is: a `MeetLine`
    /// and a `QuadVal` are an order of magnitude wider than a plane pin's two words, and the
    /// arrangement stores one of these **per segment endpoint** while every question on the hot
    /// path is a plane pin. Unboxed it moved `split_at_crossings` from 23.2% of a rotated fold to
    /// 28.2%, and 12.1% of an axis-aligned one to 19.0% — with the *same trip counts*, so the cost
    /// was the width and nothing else. One allocation on the cold road buys that back.
    Pierce { name: NodeId, at: Box<PierceAt> },
}

/// **What a pierce point carries to be ordered** — the cold payload of [`OnLine::Pierce`] and
/// [`Located::Pierce`].
pub(crate) struct PierceAt {
    /// The meet line of the name's two planes and the root along it.
    pub(crate) meet: PiercePoint,
    /// The classes the gate proved the point's whole cylinder beside, with the side
    /// ([`crate::planes::WorkingCyl::beside`]) — the side of such a class this point lies on, which
    /// its coefficients cannot say: the class has no world statement. Empty, and unallocated,
    /// unless the boolean has one.
    pub(crate) beside: Vec<(usize, i8)>,
}

impl PierceAt {
    fn of(
        jd: &Judge<'_, WorkingPlane>,
        cy: &crate::planes::WorkingCyl,
        cyl: usize,
        name: NodeId,
    ) -> Option<Box<PierceAt>> {
        Some(Box::new(PierceAt {
            meet: pierce_meet(jd, cyl, &cy.def, name)?,
            beside: cy.beside_all().collect(),
        }))
    }

    /// The side of plane `q` this point lies on, in [`side_of`]'s convention — the gate's record
    /// where it holds one, else `q`'s statement ([`PiercePoint::side`]); one reading, as `side_of`'s
    /// pierce arm makes it.
    fn side(&self, jd: &Judge<'_, WorkingPlane>, q: usize) -> Option<i8> {
        match self.beside.iter().find(|&&(c, _)| c == q) {
            Some(&(_, side)) => Some(side),
            None => self.meet.side(jd, q),
        }
    }
}

/// A point in the form a comparison's **first** argument needs: an [`OnLine`] plus the Cramer
/// handle the plane road asks `orient3d` of.
///
/// **The asymmetry is measured, not assumed.** Making both sides `Located` for the symmetry does
/// not pay: an endpoint's handle is **never read** — the plane road is
/// `at_a.orient3d(b.pin) × ds_b` — while `ImplicitPoint` is wide enough (a `OnceCell` of Cramer
/// parts) that storing two per segment moves `split_at_crossings` from 23.2% of a rotated fold to
/// 27.2%, and 12.1% of an axis-aligned one to 17.8%. First arguments are built a few times per
/// wall; second arguments are stored per segment and asked millions of times. What the symmetry
/// would protect holds anyway: `ds` is built beside the `pin` it belongs to, by one constructor, so
/// it cannot disagree with the `(p, q)` it was made for.
pub(crate) enum Located<'a> {
    /// **No `ds` here, and that is the point.** The plane road is
    /// `at_a.orient3d(b.pin) × ds_b` — a *first* argument's `dir_sign` is never read. Computing it
    /// anyway cost a `plane_pair_dir_sign` per wall pair and per sub-interval: measured, the cover
    /// loop went 48.6→90.8ms on an axis-aligned fold with the trip counts unchanged.
    Class {
        pin: usize,
        at: crate::tolerant::ImplicitPoint<'a, WorkingPlane>,
    },
    Pierce {
        name: NodeId,
        at: Box<PierceAt>,
    },
}

/// Whether two located points are **the same point by identity**, without asking a predicate.
///
/// Two plane-pinned points compare by **pin**, which is the integer test the arrangement's
/// overlay has always used: on one line `{p, q, pin}` is injective in `pin`, so equal pins are one
/// point. It is a *sufficient* condition and not the whole test — where four planes meet, one point
/// wears two pins — which is why callers still ask the order afterwards.
/// Anything a cylinder pins has no pin payload to compare, so those compare by **name**.
#[inline]
pub(crate) fn same_point(a: &Located<'_>, b: &OnLine) -> bool {
    match (a, b) {
        (Located::Class { pin: x, .. }, OnLine::Class { pin: y, .. }) => x == y,
        (Located::Pierce { name: x, .. }, OnLine::Pierce { name: y, .. }) => x == y,
        // A point three planes name and a point a cylinder pins are never the same *name* — the
        // two are different variants of `NodeId`. Whether they are the same **point** is a
        // different question, and the order predicate below is what asks it.
        _ => false,
    }
}

/// **Is `at` within the closed extent an edge's two ends mark out on `L = P ∩ Q`?** — endpoints
/// included.
///
/// On an endpoint it is contained, checked by **identity**, because `order_along(x, x)` is not
/// defined to return 0. Otherwise it
/// is contained iff it lies between the two — opposite order signs — or a `0` says it *is* one of
/// them by a second name.
///
/// **Asked before the predicates, not after.** Identity is a *sufficient* condition for
/// containment, so answering it first skips both orientations. It is not the whole test: where four
/// planes meet, one point wears two pins, and `at` may be the group's representative while the edge
/// still remembers the other — which is what the `== 0` arms catch. Subsumed, not dropped; the
/// order between them is free.
///
/// **`ends` must have been located for the very `(p, q)` handed in here.** Each [`OnLine`]
/// carries a `dir_sign` — or a `(line, s)` — made for one plane pair, and asking about a point
/// located for a *different* pair silently compares two different lines. The type does not say so:
/// `OnLine` records no pair, and giving it one widens the vector the overlay's cover loop reads
/// millions of times. **Callers hold the pair beside the ends.**
///
/// **The point is a parameter, and that is what lets a caller hoist it.** The point asked about
/// is constant for every edge on one wall, so a caller locates it once per **wall pair** and every
/// edge there shares its Cramer parts. Taking a bare class instead would cap the sharing at one
/// edge's two endpoints.
///
/// `None` is a description that could not be formed exactly, never a shape this cannot answer;
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
/// **Infallible for a plane pin** — `Judge::point` and [`dir_sign`] are total — so a plane pin
/// never makes a comparison fail.
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
            // A `match` and not a fallback: a third variant must light this up rather than
            // fall in with the three-plane one (`winding`'s `coord_key` states the rule).
            let cyl = match name.kind() {
                NodeKind::Pierce { cyl, .. } => cyl,
                NodeKind::ThreePlane(_) => return None,
            };
            debug_assert!(names_the_line(name, p, q), "{name:?} is not on {p} ∩ {q}");
            Some(OnLine::Pierce {
                name,
                at: PierceAt::of(jd, cyls.get(cyl)?, cyl, name)?,
            })
        }
    }
}

/// Whether a pierce point's name puts it on the line `p ∩ q` — the pair the two locators are handed.
///
/// **A producer's invariant, asserted rather than trusted.** A pierce point is located as one
/// only on the line its own pair names: a site that mints the point on a line names it by that
/// line's pair, and a site that pins an end arriving under a name asks [`pin_for`], which answers
/// the quadric for the name's own pair alone and the point's other plane for any other line. The
/// quadratic road leans on it: the point's meet line is solved from the name's planes, and the axis two pierce points are
/// compared on comes from the line's (`stored_line_dir`). A pair that disagreed would compare a
/// point against a line it does not lie on. `false` for a three-plane name, which no pierce arm hands in.
fn names_the_line(name: NodeId, p: usize, q: usize) -> bool {
    pierce_name(name).is_some_and(|(mut planes, ..)| {
        planes.sort_unstable();
        planes == [p.min(q), p.max(q)]
    })
}

/// Put a point into the **first**-argument form — an [`on_line`] with the Cramer handle beside it.
///
/// **Infallible for a plane pin**, like [`on_line`].
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
            // A `match` and not a fallback: a third variant must light this up rather than
            // fall in with the three-plane one (`winding`'s `coord_key` states the rule).
            let cyl = match name.kind() {
                NodeKind::Pierce { cyl, .. } => cyl,
                NodeKind::ThreePlane(_) => return None,
            };
            debug_assert!(names_the_line(name, p, q), "{name:?} is not on {p} ∩ {q}");
            Some(Located::Pierce {
                name,
                at: PierceAt::of(jd, cyls.get(cyl)?, cyl, name)?,
            })
        }
    }
}

/// **The order of two located points along `L = P ∩ Q`** — [`order_on`]'s road for a pair with a
/// pierce point, and the whole rule for a caller that hands its hoists in rather than rebuilding them.
///
/// The plane/plane test is the **first statement** on purpose: everything after it computes
/// the line's rational direction and a radical comparison, and the overlay asks this millions of
/// times on pairs that never reach there.
/// **Hot and cold are separate functions, and the profile is what says so.** With the
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
    // `order_along(jd, p, q, i, j)` exactly: `ImplicitPoint::orient3d(j)` is
    // `Judge::orient3d(p, q, pin, j)` on the same four planes with the Cramer parts cached, and
    // `ds` is that call's `dir_sign` factor.
    if let (Located::Class { at, .. }, OnLine::Class { pin: j, ds, .. }) = (a, b) {
        return Some(at.orient3d(*j) * ds);
    }
    order_located_quad(jd, p, q, a, b)
}

/// [`order_located`]'s other road: anything a cylinder pins. Kept out of line so the plane road
/// above stays small enough to be.
///
/// **A pierce point against a plane pin is the plane road's own formula.** The plane road's
/// `at.orient3d(j) · ds` is «the side of the first point against `j`, times `dir_sign(p, q, j)`»
/// ([`order_pinned`]'s derivation — `j`'s orientation cancels as a square), and a pierce point has a
/// side of a plane too ([`PiercePoint::side`]). So the mixed pair reads no coordinate of the plane-pinned
/// point and no direction of the line: `(X, j) ↦ side_X(j) · ds_j` and `(i, X) ↦ −side_X(i) · ds_i`.
/// `None` where `side_X` has no description to read. Only two pierce points compare on a
/// coordinate axis through the `a + b√c` tower, below.
#[inline(never)]
fn order_located_quad(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    a: &Located<'_>,
    b: &OnLine,
) -> Option<i8> {
    use nacre_exact::Orient;
    use nacre_exact::quad::cmp_coord_branch;
    let (m1, m2, meets) = match (a, b) {
        // Unreachable — two plane-pinned points returned through the integer road above. Spelled
        // as a decline rather than a panic or a second rational comparison nobody would exercise.
        (Located::Class { .. }, OnLine::Class { .. }) => return None,
        (Located::Pierce { at, .. }, OnLine::Class { pin, ds }) => {
            return Some(at.side(jd, *pin)? * ds);
        }
        (Located::Class { pin, .. }, OnLine::Pierce { at, .. }) => {
            return Some(-at.side(jd, *pin)? * dir_sign(jd, p, q, *pin));
        }
        (Located::Pierce { at: a1, .. }, OnLine::Pierce { at: a2, .. }) => {
            match (a1.meet.exact(), a2.meet.exact()) {
                (Some(m1), Some(m2)) => (m1, m2, (&a1.meet, &a2.meet)),
                // A point with no exact meet is ordered by its interval box along the line's
                // stored direction — the same `sign((a − b)·d)` the exact arm reads per axis.
                _ => return order_by_boxes(jd, p, q, &a1.meet, &a2.meet),
            }
        }
    };
    let sign = |o: Orient| -> i8 {
        match o {
            Orient::Positive => 1,
            Orient::Negative => -1,
            Orient::Zero => 0,
        }
    };
    // The axis is the first component the line's direction has, and the order along the line is
    // the order on that axis times that component's sign. The direction is the rational one:
    // on an axis where the true component is `0` every point of the line has the same coordinate,
    // so an axis picked off the rounded plane caches there reads two distinct points as one.
    let d = stored_line_dir(jd, p, q)?;
    let zero = nacre_exact::Rat::from_int(0);
    let (k, dsign) = (0..3).find_map(|k| {
        let s = if d[k] > zero {
            1i8
        } else if d[k] < zero {
            -1
        } else {
            0
        };
        (s != 0).then_some((k, s))
    })?;
    let order = sign(cmp_coord_branch(m1, m2, k)) * dsign;
    // The box road's convention (`sign((a − b)·d)` along the stored direction), checked against
    // the exact one wherever both answer — the road a pair with a nameless point takes, so a
    // difference between the two would order such a point against its named neighbours wrongly.
    debug_assert!(
        order_by_boxes(jd, p, q, meets.0, meets.1).is_none_or(|b| b == order),
        "the interval order along a line disagrees with the exact one"
    );
    Some(order)
}

/// `sign((a − b) · d)` for two points on the line `p ∩ q`, `d` its stored direction
/// (`n_p × n_q`, [`crate::planes::WorkingPlane::stored_hp`]) — read off interval boxes. `None`
/// where they cannot separate the two: distinct names, so never a tie.
fn order_by_boxes(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
    a: &PiercePoint,
    b: &PiercePoint,
) -> Option<i8> {
    let prec = jd.standard.prec;
    let (x, y) = (a.point_hp(jd)?, b.point_hp(jd)?);
    let (np, nq) = (jd.stored_hp(p, prec), jd.stored_hp(q, prec));
    let d = [
        np[1].mul(&nq[2], prec).sub(&np[2].mul(&nq[1], prec), prec),
        np[2].mul(&nq[0], prec).sub(&np[0].mul(&nq[2], prec), prec),
        np[0].mul(&nq[1], prec).sub(&np[1].mul(&nq[0], prec), prec),
    ];
    (0..3)
        .map(|k| x[k].sub(&y[k], prec).mul(&d[k], prec))
        .reduce(|s, t| s.add(&t, prec))?
        .sign()
        .map(|pos| if pos { 1 } else { -1 })
}
