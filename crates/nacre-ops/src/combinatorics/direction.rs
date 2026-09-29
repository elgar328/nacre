use super::*;
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
/// bare `i8`, so no caller pairs them by hand. (`EdgeDir::new`
/// still states a pair directly — that is for fixtures, and it is the one place a wrong pairing is
/// still spellable.)
///
/// ★★ **Made once.** `order_along` is `sign((V_i − V_j)·d)`, so the direction of travel is either
/// `−order_along(from, to)` (invert the result) or `order_along(to, from)` (swap the arguments) —
/// the same value by the antisymmetry of a difference. Two spellings would agree only because two
/// independent inversions cancel, and a change to `order_along`'s convention would move one and
/// not the other.
///
/// ★ **The zero policy lives here, and that is not a matter of taste**: both callers answered a
/// coincidence the same way (`CoincidentNodes`). Where two consumers want *different* answers —
/// the turn's zero, which is a rejection to one and a bucket to the other — the policy stays with
/// them and only the sign is shared (see [`turn`]).
///
/// ★★ **Both ends arrive as a name *and* a pin, because a cylinder-pinned one needs both.** The
/// pin says which kind of thing holds the point on `L`; the name says which point. A plane pin
/// carries its own name in its payload, so for it one argument would do — but
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
    // ★★ **`None` is `WitnessNotRational`, not `RingNaming`.** A cylinder-pinned end is answered,
    // so what reaches here is [`order_pinned`]'s own `None`: a class with no narrow rational
    // description, a coordinate past `Rat`. The names chain perfectly; the *value* could not be
    // formed. That is the sentence `arc_at` next door uses for the same cause.
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
/// comparison (measured on `split_at_crossings`, which asks the same kind of question: 1.5M
/// `dir_sign` calls where 113k are distinct).
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
        // *cannot* answer a cut end — that was true while a pierce point had no order, and
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
        // order is asked (its ends are pierce points, which have no third plane to order by).
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
/// Everything here is a *cache of the node's name*: [`pierce_meet`] re-solves the point from the
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
    // ★★ **Two causes, two names.** `pierce_meet` folds four `None`s into one, and they are not the
    // same fact: a three-plane node on an arc is a *naming* failure (the split shipped the
    // two-names case it exists to refuse), while a class with no rational description, a meet that
    // does not solve, and a root the meet does not have are all "the exact route declined". Asking
    // the kind first is what keeps `RingNaming`'s sentence true where it is raised.
    //
    // ★ Only the first is reachable from here: the other three would have stopped the split that
    // built this arc, since it re-solves the same pair of classes for the same cylinder.
    // ★ A `match` and not a fallback: a third variant must light this up rather than be
    // rejected under a name that is not its own (`winding`'s `coord_key` states the rule).
    match node.kind() {
        NodeKind::Pierce { .. } => {}
        NodeKind::ThreePlane(_) => return Err(reject(RejectReason::RingNaming)),
    }
    let at = pierce_meet(jd, a.cyl, &a.def, node).ok_or_else(undecided)?;
    let coeffs = class_coeffs_rat(jd, p).ok_or_else(undecided)?;
    // The circle's centre: where the axis pierces this plane. ★ The canonical sign the coefficients
    // carry cancels in the parameter (numerator and denominator both flip), so this one does not
    // need the stored-frame turn that `arc_side` does.
    let t = crate::planes::axis_param_of_plane(&coeffs, &a.def).ok_or_else(undecided)?;
    let (o, m) = (a.def.origin(), a.def.dir());
    let centre: [nacre_exact::Rat; 3] = (|| {
        let mut c = [nacre_exact::Rat::from_int(0); 3];
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
        axis: m,
        axis_up: crate::planes::plus_t_is_above(
            jd.planes[p].world.as_ref().ok_or_else(undecided)?,
            &a.def,
        ),
        ccw: a.ccw,
    })))
}

/// A class's exact description **turned to face the way the plane does** — the frame every
/// direction sign in this file is written in ([`turn`]'s `frame_sign` bridge, the cell labels'
/// «above»).
///
/// ★★★ [`class_coeffs_rat`] hands back the class's *canonical* name — first nonzero component
/// positive — which points the **other way** from the plane on half the classes, so a rule spelled
/// against it reads backwards on exactly those classes. The turn is the truth's
/// ([`WorkingPlane::stored_world_rat`] — the world name's sense, read off the plane's `sense` and
/// points), never a comparison with the plane cache: that cache is a rounded image, and the name
/// needs no second description to say which way it faces.
pub(crate) fn stored_coeffs_rat(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
) -> Option<[nacre_exact::Rat; 4]> {
    jd.planes[c].stored_world_rat()
}

/// **The direction of the line `L = P ∩ Q`, exactly** — `n_p × n_q` over the two classes' stored
/// normals ([`stored_coeffs_rat`]), which is the direction [`order_along`] orders the line by. One
/// spelling for every reader that needs the line's direction as a vector rather than as a sign
/// against a third plane: the quadratic order's axis, an arc's tangency, an arc's side.
///
/// ★★ **Rational, never read off the plane caches.** A component the truth makes exactly `0`
/// stays `0` here; a cross product of the rounded caches can leave it at a few ulps, and a reader
/// that picks "the first nonzero component" as its axis then compares points that all share that
/// coordinate — two distinct points read as one. `None` when a class has no world description or
/// the product overflows.
pub(crate) fn stored_line_dir(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
) -> Option<[nacre_exact::Rat; 3]> {
    let (np, nq) = (stored_coeffs_rat(jd, p)?, stored_coeffs_rat(jd, q)?);
    nacre_exact::cross3_rat(&[np[0], np[1], np[2]], &[nq[0], nq[1], nq[2]])
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
// Production reads the turn through `turn_between`, which lets a caller skip a straight
// stretch, so this per-node spelling has no product caller at all -- `loop_winding`'s tests
// are the only ones, and they ask it at a hull vertex where the turn and the winding agree.
#[cfg(test)]
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
/// ★★ **Measured: nothing in the suite reaches `DoublesBack`** — `straight_angle` is
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
) -> Result<Continuation, BoolError> {
    Ok(match (earlier, later) {
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
        // ★ A line and an arc **tangent** at the shared node (a fillet's smooth
        // corner): the turn is `0`, and whether the ring runs on or doubles back is the travel
        // directions' agreement ([`tangent_travel_agrees`]).
        (EdgeDir::Line { carrier, sense }, EdgeDir::Arc(arc))
        | (EdgeDir::Arc(arc), EdgeDir::Line { carrier, sense })
            if turn(jd, p, earlier, later)? == 0 =>
        {
            match tangent_travel_agrees(jd, p, *carrier, *sense, arc) {
                Some(true) => Continuation::Straight,
                Some(false) => Continuation::DoublesBack,
                None => return Err(reject(RejectReason::WitnessNotRational)),
            }
        }
        // ★ A line and a ruling that are **one** line (a prism's edge running on past the
        // cylinder's cap along a ruling): the straight reading once more, by the travel
        // directions ([`line_along_ruling`]).
        (EdgeDir::Line { .. }, EdgeDir::Ruling(_)) | (EdgeDir::Ruling(_), EdgeDir::Line { .. })
            if turn(jd, p, earlier, later)? == 0 =>
        {
            match line_along_ruling(jd, p, earlier, later)? {
                Some(true) => Continuation::Straight,
                Some(false) => Continuation::DoublesBack,
                None => return Err(reject(RejectReason::WitnessNotRational)),
            }
        }
        // A line and an arc at a **transversal** crossing turn — that is what transversal means.
        _ => Continuation::Turns,
    })
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
pub(super) fn turn_between(
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
/// **here and nowhere else**, and this function is the one place the product computes it.
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
/// whose sign closes over `QuadVal::sign` (the tangent `m × (p − o)` is *linear* in the pierce
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
/// ★★★ **Three more places read the direction's representation, and they are not this atom.** Each
/// is a named function beside this one, so the next widening is a `match` the compiler points at
/// rather than a list a reader has to trust:
/// - [`antiparallel`] — the π pole in `angular_order` ("same wall, opposite sign"), which for arcs
///   becomes "same circle, opposite tangent";
/// - [`parallel_carriers`] — [`loop_winding`]'s walk-back, which becomes "are the tangents
///   parallel";
/// - [`parallel_carriers`] again — `arrangement::split_at_crossings`' **wall direction families**
///   (`Wall.dir`). It sits in a different file and asks the same
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
        // A circle bound's plane is ⊥ to the axis, so `n_P ∥ m` — see
        // [`class_carries_circle`] for what actually holds that (the gate, not this module). The
        // arc's tangent is
        // `T = ±(m × r)` with `r = x − c`, the segment's direction
        // `d = n_P × n_carrier` is ⊥ to `m`, and BAC-CAB collapses the cross product:
        //
        //   (d × T) · m = (d × (m × r)) · m = (m (d·r) − r (d·m)) · m = (d·r)(m·m)
        //
        // so the whole turn is `sign(d · (x − c))` — **one** `a + b√c` question, and
        // `quad::plane_side` is exactly it: the plane with normal `d` through the centre, measured
        // at the pierce point.
        (EdgeDir::Line { carrier, sense }, EdgeDir::Arc(arc)) => {
            sense * arc_side(jd, p, *carrier, arc)?
        }
        (EdgeDir::Arc(arc), EdgeDir::Line { carrier, sense }) => {
            -sense * arc_side(jd, p, *carrier, arc)?
        }
        // Two arcs of one circle at one node are tangent: no turn to read. (Two *different*
        // circles cannot meet at a node: faces of different operands are proved apart by the
        // gate, per face pair, and one valid operand's own faces do not cross —
        // the sketch refuses an arc–arc join, `ArcsMeetAtVertex`.)
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
        // faces the gate proved apart (different operands) or faces of one valid
        // operand, which do not cross. `0` sends a surprise to the bucket whose
        // walk names it (`UnorderedEdges`) rather than ranking it.
        (EdgeDir::Ruling(_), EdgeDir::Arc(_)) | (EdgeDir::Arc(_), EdgeDir::Ruling(_)) => 0,
    })
}

/// `turn(line, ruling)` with the line's sense factored out — the `−sign(m · n_carrier)` the
/// derivation above collapses to, times the ruling's travel.
///
/// ★ **The overall sign is measured**: negating it turns the through-boss volume
/// oracles red (all five production fixtures). The BAC-CAB derivation alone would not lock it —
/// the walk's both-handedness try
/// absorbs a global flip; the panel population reads it for real.
fn ruling_line_turn(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    carrier: usize,
    r: &RulingDir,
) -> Result<i8, BoolError> {
    debug_assert!(
        class_coeffs_rat(jd, p)
            .map(|n| {
                nacre_exact::dot_sign_rat(&[n[0], n[1], n[2]], &r.axis) == nacre_exact::Orient::Zero
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
    let dot = nacre_exact::dot_sign_rat(&[n[0], n[1], n[2]], &r.axis);
    let up = if r.up { 1 } else { -1 };
    Ok(match dot {
        nacre_exact::Orient::Positive => -up,
        nacre_exact::Orient::Negative => up,
        // `m · n_carrier = 0` means the segment's line is parallel to the ruling — no crossing
        // could have put them at one node, so a `0` here is the degenerate concurrency the
        // walk's `UnorderedEdges` bucket names.
        nacre_exact::Orient::Zero => 0,
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
/// | the canonical→stored turn | use [`class_coeffs_rat`] | **locked** — one class in the corpus
/// disagrees, and that class's walk merges four cells into one 8-half-edge orbit. ★ The contour
/// count passes it; `arrangement::walk_cells`' Euler condition is what refuses it |
/// | `axis_up` | drop it | **locked** — the unmoved corpus already: the disk-side watcher (`disk_side_probe::record`, *"the disk-side rule and the cell's own corners disagree"*) dies on the first boolean |
/// | `frame_sign` | drop it | **locked** — the commuting oracle's always-on subset turns red on 250 of 396 cells (`t(−4,−4,−2)` puts the plate's caps on a seed plane, `frame_sign = −1`), at the same watcher |
///
/// The last two are derived, not guessed (the algebra is above). ★ Both rows read
/// «unexercised — `+1`/`true` on every class reached» without a fixture made for them: a face
/// lies on a seed plane
/// with its outward along +axis under a translation as ordinary as `t(−4,−4,−2)`, but the
/// corpus does not contain the population a rule needs until a fixture adds it.
///
/// ★★★ **The `sense` a split carries is locked too.** Flipping the sense a
/// split carries onto its sub-segments (`Carrier::Plane::sense`) attaches the arcs to the wrong
/// cells, and nothing in the walk sees it — the cell count and the contour count both come out
/// right. Its first reader is neither `nest_cells`' root choice nor
/// `label_cells`' keep decision (every order-independent summary of the labels is identical
/// because the two 3-cells *swap* labels). It is `emit_faces`, where the ring comes out the exact
/// reverse, and `bands`' arc fence pins it there through `ClassAudit::outer_rings`.
///
/// ★ The table's first row holds against that lock too: negating the
/// whole result leaves the whole crate green. A global flip really is absorbed.
/// **At a node where a line and an arc are tangent, do their travel directions agree?** — the
/// sign of `d · t`: `d` the line's travel direction (`sense` along [`stored_line_dir`] — the
/// direction [`order_pinned`] orders by and [`arc_side`] reads), and
/// `t = way · (m × (N − c))` the arc's travel tangent at its node `N` (`way` is `+1` for
/// counter-clockwise travel about the axis `m`), formed in one radical from `N = base + s·dir`.
///
/// Asked only where [`turn`] is `0` for the pair — the line is tangent to the circle at `N`, so
/// the two directions are parallel and the dot decides. `Some(false)`: **opposite** — a smooth
/// join, the ring runs straight through, two departures are a half turn apart. `Some(true)`:
/// the **same** way — the ring doubles back along the arc, or two departures coincide, which is a
/// **curvature** question this crate does not order yet (two tangent circles, M6b's shape). `None`
/// is a zero dot (not tangent after all) or overflow.
///
/// ★ A fillet's or a slot's wall meets its cylinder exactly so, at a corner whose
/// root is `Double`; the angular order reads the pair as a tie (`UnorderedEdges`) and the winding
/// walk as doubling back (`StraightAngle`) because [`antiparallel`] has no line–arc arm — the
/// structure cannot tell, only the geometry can, and this is where it is asked once.
pub(crate) fn tangent_travel_agrees(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    carrier: usize,
    sense: i8,
    arc: &ArcDir,
) -> Option<bool> {
    use nacre_exact::{Orient, quad::QuadVal};
    let d = stored_line_dir(jd, p, carrier)?;
    let (line, s) = &arc.at;
    let mut rel = line.base();
    for (r, c) in rel.iter_mut().zip(arc.centre.iter()) {
        *r = r.checked_sub(*c)?;
    }
    let u = nacre_exact::cross3_rat(&arc.axis, &rel)?;
    let v = nacre_exact::cross3_rat(&arc.axis, &line.dir())?;
    let dot = QuadVal::from_rat(nacre_exact::dot3_rat(&d, &u)?)
        .checked_add(&s.checked_mul_rat(nacre_exact::dot3_rat(&d, &v)?)?)?;
    let way = if arc.ccw { 1i8 } else { -1 };
    match dot.sign() {
        Orient::Positive => Some(sense * way > 0),
        Orient::Negative => Some(sense * way < 0),
        Orient::Zero => None,
    }
}

/// Whether two **departures** from one node are a half turn apart by tangency — a line and an
/// arc tangent at the node, the arc leaving opposite to the line. The structural
/// [`antiparallel`] cannot see it; [`tangent_travel_agrees`] can. `Ok(false)` for any other pair.
pub(crate) fn tangent_pole(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<bool, BoolError> {
    let (line, arc) = match (a, b) {
        (EdgeDir::Line { carrier, sense }, EdgeDir::Arc(arc))
        | (EdgeDir::Arc(arc), EdgeDir::Line { carrier, sense }) => ((*carrier, *sense), arc),
        _ => return Ok(false),
    };
    if turn(jd, p, a, b)? != 0 {
        return Ok(false);
    }
    Ok(tangent_travel_agrees(jd, p, line.0, line.1, arc) == Some(false))
}

/// Whether two **departures** from one node are a half turn apart along one straight line — a
/// plane-pair line and a ruling that are the same line (two planes along a cylinder's axis meeting
/// on it, [`crate::planes::SharedRuling`]), leaving opposite ways: the segment that runs on past
/// the lateral's cap, and the ruling piece it continues. The structural [`antiparallel`] cannot
/// see it (one carrier is a plane class, the other a cylinder); the directions can — the line's
/// is `sense · (n_P × n_carrier)` in the stored frame ([`stored_line_dir`], the frame `sense` is
/// made in), the ruling's `±m`. `Ok(false)` for any other pair.
pub(crate) fn ruling_pole(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<bool, BoolError> {
    Ok(line_along_ruling(jd, p, a, b)? == Some(false))
}

/// For a line and a ruling at one node that lie on **one** line (their turn is `0`): whether they
/// travel the same way. `None` for any other pair.
fn line_along_ruling(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<Option<bool>, BoolError> {
    let ((carrier, sense), r) = match (a, b) {
        (EdgeDir::Line { carrier, sense }, EdgeDir::Ruling(r))
        | (EdgeDir::Ruling(r), EdgeDir::Line { carrier, sense }) => ((*carrier, *sense), r),
        _ => return Ok(None),
    };
    if turn(jd, p, a, b)? != 0 {
        return Ok(None);
    }
    let d =
        stored_line_dir(jd, p, carrier).ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
    let along = match nacre_exact::dot_sign_rat(&d, &r.axis) {
        nacre_exact::Orient::Positive => sense,
        nacre_exact::Orient::Negative => -sense,
        // A zero turn puts the line along the axis, so `d · m` cannot vanish; if it does, the
        // pair is not the shape this answers.
        nacre_exact::Orient::Zero => return Ok(None),
    };
    Ok(Some((along > 0) == r.up))
}

fn arc_side(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    carrier: usize,
    arc: &ArcDir,
) -> Result<i8, BoolError> {
    use nacre_exact::{Orient, Rat};
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
    let d = stored_line_dir(jd, p, carrier).ok_or_else(wide)?;
    let plane = (|| -> Option<[Rat; 4]> {
        Some([
            d[0],
            d[1],
            d[2],
            Rat::from_int(0).checked_sub(nacre_exact::dot3_rat(&d, centre)?)?,
        ])
    })()
    .ok_or_else(wide)?;
    let side = match nacre_exact::quad::plane_side(&plane, &at.0, &at.1) {
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
