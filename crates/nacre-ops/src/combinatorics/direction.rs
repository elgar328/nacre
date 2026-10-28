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
/// It returns the sense **paired with the carrier it belongs to** ([`EdgeDir`]) rather than a
/// bare `i8`, so no caller pairs them by hand. (`EdgeDir::new`
/// still states a pair directly — that is for fixtures, and it is the one place a wrong pairing is
/// still spellable.)
///
/// **Made once.** `order_along` is `sign((V_i − V_j)·d)`, so the direction of travel is either
/// `−order_along(from, to)` (invert the result) or `order_along(to, from)` (swap the arguments) —
/// the same value by the antisymmetry of a difference. Two spellings would agree only because two
/// independent inversions cancel, and a change to `order_along`'s convention would move one and
/// not the other.
///
/// **The zero policy lives here, and that is not a matter of taste**: both callers answered a
/// coincidence the same way (`CoincidentNodes`). Where two consumers want *different* answers —
/// the turn's zero, which is a rejection to one and a bucket to the other — the policy stays with
/// them and only the sign is shared (see [`turn`]).
///
/// **Both ends arrive as a name *and* a pin, because a cylinder-pinned one needs both.** The
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
    // **`None` is `WitnessNotRational`, not `RingNaming`.** A cylinder-pinned end is answered,
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
/// **`node` is not decoration.** A straight edge's tangent is the same at both ends, so for it "the
/// edge's direction" and "the direction *at* this end" are one thing. An arc's two ends differ —
/// and worse, a direction taken at a node the edge does not even touch is a perfectly legal call
/// that returns a confident, wrong sign.
///
/// So the node comes in and is **checked here, at the one place a direction is born** — not at
/// `turn`, which would have to become fallible and drag `Result` through two sorts, and not at the
/// callers, who would each have to remember. `debug_assert` because "a ring is a chain" is a
/// *producer's* invariant, not an input's: the census runs in both profiles, so the debug one is
/// the instrument.
///
/// **The hoist survives.** `angular_order` builds one direction per outgoing half-edge and sorts
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
        // The carried sense first, and only then the endpoints — which can answer a cut end too
        // ([`order_pinned`] orders a pierce point); the better reason is that the sense is the
        // **splitter's own statement** about a piece it made, taken once for the whole segment and
        // handed to every sub-segment, rather than re-derived per piece from two names.
        // The coincident-endpoint check `edge_dir` makes is not lost with it: the split refuses
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
    // **Two causes, two names.** `pierce_meet` folds four causes into one `None` (its doc), and
    // one of them is a *naming* fact: a three-plane node on an arc means the split shipped the
    // two-names case it exists to refuse. Asking the kind first is what keeps `RingNaming`'s
    // sentence true where it is raised.
    //
    // The other three do not reach here: they would have stopped the split that built this arc,
    // since it re-solves the same pair of classes for the same cylinder.
    // A `match` and not a fallback: a third variant must light this up rather than be
    // rejected under a name that is not its own (`winding`'s `coord_key` states the rule).
    match node.kind() {
        NodeKind::Pierce { .. } => {}
        NodeKind::ThreePlane(_) => return Err(reject(RejectReason::RingNaming)),
    }
    let at = pierce_meet(jd, a.cyl, &a.def, node)
        .and_then(PiercePoint::into_exact)
        .ok_or_else(undecided)?;
    let coeffs = jd.world_rat(p).ok_or_else(undecided)?;
    // The circle's centre: where the axis pierces this plane. The canonical sign the coefficients
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

/// **The direction of the line `L = P ∩ Q`, exactly** — `n_p × n_q` over the two classes' stored
/// normals ([`Judge::stored_world_rat`]), which is the direction [`order_along`] orders the line by. One
/// spelling for every reader that needs the line's direction as a vector rather than as a sign
/// against a third plane: the quadratic order's axis, an arc's tangency, an arc's side.
///
/// **Rational, never read off the plane caches.** A component the truth makes exactly `0`
/// stays `0` here; a cross product of the rounded caches can leave it at a few ulps, and a reader
/// that picks "the first nonzero component" as its axis then compares points that all share that
/// coordinate — two distinct points read as one. `None` when a class has no world description or
/// the product overflows.
pub(crate) fn stored_line_dir(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    q: usize,
) -> Option<[nacre_exact::Rat; 3]> {
    let (np, nq) = (jd.stored_world_rat(p)?, jd.stored_world_rat(q)?);
    nacre_exact::cross3_rat(&[np[0], np[1], np[2]], &[nq[0], nq[1], nq[2]])
}

/// The turn at ring node `i`, about the face's **outward** normal: `+1` left, `-1` right.
///
/// **No point is materialized, and no coordinate is read** — the algebra that makes that true
/// lives with the atom, [`turn`], and is not repeated here.
///
/// **It is `0` exactly when the two walls name the same plane.** The node is `P ∩ A ∩ B`, and a
/// point exists there only if the three normals are independent — so with two genuinely different
/// planes the determinant cannot vanish. Two edges through one node whose lines are parallel are
/// two edges on **one** line, and that is the straight stretch a ring can arrive as (the
/// arrangement names a point wherever a feature crosses an edge).
///
/// **"Same plane" is not "same class".** Aliasing lets two classes name one plane, which is why
/// `angular_order` folds them (`Aliases::union_wall`) and still keeps an honest `UnorderedEdges`
/// for the pair aliasing did not reach.
///
/// A ring is not convex, so this is **not** the winding — at a reflex node it is its
/// opposite. [`loop_winding`] asks it at a hull vertex, where the two agree.
// No product caller: `loop_winding` orders the extreme node's departures instead. Its tests ask
// this at a hull vertex, where the turn and the winding agree.
#[cfg(test)]
pub(crate) fn turn_at(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    p: usize,
    ring: &[RingEdge],
    i: usize,
) -> Result<i8, BoolError> {
    let n = ring.len();
    // Both read **at this node** — the one they share. For a straight edge that is today's
    // answer either way; naming it is what makes an arc's two ends tellable apart.
    let arriving = dir_at(jd, cyls, p, &ring[(i + n - 1) % n], ring[i].node)?;
    let leaving = dir_at(jd, cyls, p, &ring[i], ring[i].node)?;
    match turn(jd, p, &arriving, &leaving)? {
        0 => Err(reject(RejectReason::StraightAngle)),
        t => Ok(t),
    }
}

/// **Do two carriers give one direction on `P`?** — `P ∩ a` and `P ∩ b` are parallel.
///
/// **`parallel`, not `collinear`, and the difference is the caller's.** Two parallel lines are
/// the same line only when they share a point; `arrangement`'s wall direction families ask it of
/// walls that share none — their own comment says so: *"Same family ⇒ the two lines are parallel
/// and meet in no point."* The premise belongs to whoever has it.
///
/// **The same primitive answers a different question elsewhere and that is left alone.**
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
/// **It is not "the angle is π"**, and the guard that makes it one stays with the caller:
/// [`angular_order`] asks this only inside the bucket where the turn already came
/// back `0`. Pulled out of that guard the predicate answers a different question, because two
/// directions on *different* carriers can also be collinear (aliasing that `union_wall` did not
/// reach), and those belong in the `0` bucket rather than the `π` one. For two directions at one
/// node it needs no guard to be true when it says yes — one carrier through a shared point is one
/// line or one circle — which is how [`loop_winding`] reads it, through [`one_carrier_pole`]: an
/// outgoing departure that is the reverse of an arrival on one carrier is an antenna.
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
        // The two arcs of one circle meeting at one node are tangent by construction, so the
        // only question left is which way each travels. **They are at one node by construction**
        // too — a fan is built from a single vertex — which is why nothing here has to compare
        // positions (the sign of `s` is no stand-in for "same node": two different nodes can share
        // it).
        (EdgeDir::Arc(x), EdgeDir::Arc(y)) => x.cyl == y.cyl && x.ccw != y.ccw,
        // One ruling, opposite travel — the straight reading of the arc arm above. Two *different*
        // rulings are parallel lines and share no node, so `(cyl, side)` identity is the "one
        // carrier" premise, same as `ca == cb` for lines.
        (EdgeDir::Ruling(x), EdgeDir::Ruling(y)) => {
            x.cyl == y.cyl && x.side == y.side && x.up != y.up
        }
        // **A line and an arc are never the π pole here, and the reason is upstream**: this is
        // asked only where the turn already came back `0`, and a `0` turn against an arc means the
        // segment is *tangent* at that node — while the split cuts only at **transversal**
        // crossings. So `false` is not a shrug: it sends the pair to the `0` bucket, where
        // `angular_order`'s `UnorderedEdges` names the surprise rather than ranking it. The same
        // sentence covers a line against a ruling (their `0` is a parallel-concurrency
        // degeneracy) and a ruling against an arc (they cannot share a node at all).
        _ => false,
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
/// where `sign(det[…])` is `plane_pair_dir_sign`, already exact. The derivation is spelled
/// **here and nowhere else**, and this function is the one place the product computes it.
///
/// **`0` comes back as `0`, on purpose.** [`angular_order`] (through [`departure_turn`])
/// buckets it as the `0`/π pole and only rejects a *pair* of them; the test-only `turn_at` calls it
/// a [`RejectReason::StraightAngle`]. Deciding here would hand one caller a rule that is not its
/// own. By the same test, [`edge_dir`]'s zero *does* live in the
/// atom — both of its callers answer it identically.
///
/// **What arcs will widen, and what they will not.** A circle's tangent rides no plane, so an
/// arc arrives as a direction this signature cannot spell: what grows is the **input type**, not
/// the skeleton around it — the identity above becomes a cross product with one quadratic factor,
/// whose sign closes over `QuadVal::sign` (the tangent `m × (p − o)` is *linear* in the pierce
/// point, so no biquadratic is needed; two tangents at one vertex would be quadratic, and that is
/// two circles meeting, already refused as `CylinderPairContact`).
///
/// **What catches a mistake here, measured.** Negating this product fails **83** tests — but
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
/// **Two more places read the direction's representation, and they are not this atom.** Each is a
/// named function beside this one, so the next widening is a `match` the compiler points at rather
/// than a list a reader has to trust:
/// - [`antiparallel`] — the π pole in `angular_order` ("same wall, opposite sign"), which for arcs
///   becomes "same circle, opposite tangent";
/// - [`ruling_pole`] — a plane-pair line and a ruling that are one line, leaving opposite ways.
///
/// ([`parallel_carriers`] answers `arrangement::split_at_crossings`' wall direction families from
/// plane classes alone, not from this representation.)
///
/// **[`turn`] and [`antiparallel`] answer for arcs, and the arity question is settled where it
/// belongs.** For a straight edge the tangent is the same at both ends, so "the edge's direction"
/// and "the direction *at* this endpoint" coincide and nothing has to tell them apart; an arc's two
/// ends differ. That is [`dir_at`]'s `node`, not a wider argument list here — the direction arrives
/// already bound to the end it was read at, so the sites above compare two directions and nothing
/// else.
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
        // **A segment against an arc needs no new primitive, and the algebra says why.**
        // A circle bound's plane is ⊥ to the axis, so `n_P ∥ m` — see
        // [`crate::planes::class_carries_circle`] for what actually holds that (the gate, not
        // this module). The arc's tangent is
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
        // gate, per face pair, and an operand whose own arcs meet at a corner — two cylinders on
        // one edge — is refused before the arrangement, `CylinderPairEdge`.)
        (EdgeDir::Arc(_), EdgeDir::Arc(_)) => 0,
        // **A segment against a ruling collapses the same way the arc arm did.** The ruling's
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
        // `sense` and the coefficients read the **same** `n_ca` ([`Judge::world_rat`], the
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
/// **The overall sign is measured**: negating it turns the through-boss volume
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
        jd.world_rat(p)
            .map(|n| {
                nacre_exact::dot_sign_rat(&[n[0], n[1], n[2]], &r.axis) == nacre_exact::Orient::Zero
            })
            .unwrap_or(true),
        "a ruling direction on a class its axis does not lie in"
    );
    // The **stored** spelling, like every factor the walk's atoms read (`arc_side`'s
    // canonical→stored turn is measured-locked): `sense` is made against the stored frame, so
    // the carrier coefficients here must be too — the canonical name opposes it on half the
    // classes.
    let Some(n) = jd.stored_world_rat(carrier) else {
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
/// **It is a `Result`, and that is deliberate.** A `0` for both "collinear" and "the arithmetic ran
/// out" reads to a caller as collinear, so a width failure comes back as a *shape* answer. The two
/// are different facts and this says so.
///
/// **Which of its factors the corpus actually locks, measured by red probe** — because "the
/// suite is green" says nothing about a sign no fixture can see:
///
/// | factor | probe | verdict |
/// |---|---|---|
/// | the whole sign | negate the result | **invisible** — this same atom feeds both the cyclic order and the winding, and the walk tries both handednesses, so a *global* flip is absorbed by trying the other one |
/// | `ccw` | drop it | **locked** — the two arcs at a crossing collapse into one bucket, `UnorderedEdges` |
/// | the canonical→stored turn | use [`Judge::world_rat`] | **locked** — one class in the corpus
/// disagrees, and that class's walk merges four cells into one 8-half-edge orbit. The contour
/// count passes it; `arrangement::walk_cells`' Euler condition is what refuses it |
/// | `axis_up` | drop it | **locked** — the unmoved corpus already: the disk-side watcher (`disk_side_probe::record`, *"the disk-side rule and the cell's own corners disagree"*) dies on the first boolean |
/// | `frame_sign` | drop it | **locked** — the commuting oracle's always-on subset turns red on 250 of 396 cells (`t(−4,−4,−2)` puts the plate's caps on a seed plane, `frame_sign = −1`), at the same watcher |
///
/// The last two are derived, not guessed (the algebra is above). Both rows read
/// «unexercised — `+1`/`true` on every class reached» without a fixture made for them: a face
/// lies on a seed plane
/// with its outward along +axis under a translation as ordinary as `t(−4,−4,−2)`, but the
/// corpus does not contain the population a rule needs until a fixture adds it.
///
/// **The `sense` a split carries is locked too.** Flipping the sense a
/// split carries onto its sub-segments (`Carrier::Plane::sense`) attaches the arcs to the wrong
/// cells, and nothing in the walk sees it — the cell count and the contour count both come out
/// right. Its first reader is neither `nest_cells`' root choice nor
/// `label_cells`' keep decision (every order-independent summary of the labels is identical
/// because the two 3-cells *swap* labels). It is `emit_faces`, where the ring comes out the exact
/// reverse, and `bands`' arc fence pins it there through `ClassAudit::outer_rings`.
///
/// The table's first row holds against that lock too: negating the
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
/// the **same** way — the ring doubles back along the arc, or two departures leave together, and
/// then which comes first is curvature's answer ([`tangent_codirected_turn`]). `None` is a zero
/// dot (not tangent after all) or overflow.
///
/// A fillet's or a slot's wall meets its cylinder exactly so, at a corner whose
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

/// **Two departures along one carrier, leaving opposite ways** — one line class, one circle, one
/// ruling ([`antiparallel`]), or a plane-pair line and the ruling it is ([`ruling_pole`]). Asked of
/// a ring's outgoing departure against an incoming one turned back, it says the two are one ray:
/// an antenna, whose contributions cancel.
pub(crate) fn one_carrier_pole(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<bool, BoolError> {
    Ok(antiparallel(a, b) || ruling_pole(jd, p, a, b)?)
}

/// **Two departures a half turn apart** — [`one_carrier_pole`], or a line and an arc tangent at
/// the node with the arc leaving the other way ([`tangent_pole`]). The `π` end of
/// [`angular_order`]'s buckets, asked where [`departure_turn`] already came back `0`. (The two
/// non-structural arms read disjoint pairs — line–ruling, line–arc — so their order is no
/// answer's.)
pub(crate) fn departure_pole(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<bool, BoolError> {
    Ok(one_carrier_pole(jd, p, a, b)? || tangent_pole(jd, p, a, b)?)
}

/// **The turn from one departure to another at their shared node, the curvature tie read.**
/// [`turn`], and where that is `0` because a line and an arc leave tangent and the same way, the
/// turn their bending makes ([`tangent_codirected_turn`]) — the arc bends off the line towards its
/// centre. The one comparator [`angular_order`] sorts a node's departures by.
pub(crate) fn departure_turn(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<i8, BoolError> {
    let t = turn(jd, p, a, b)?;
    if t != 0 {
        return Ok(t);
    }
    Ok(tangent_codirected_turn(jd, p, a, b)?.unwrap_or(0))
}

/// **A line and an arc leaving one node tangent and the same way** — the turn from `a` to `b`
/// that their curvature makes: the line runs straight on and the arc bends towards its centre, so
/// the arc lies counter-clockwise of the line (about the face's outward normal) exactly when its
/// centre is on the line's left. That is the sign [`arc_side`] multiplies a radial side by —
/// `axis_up · ccw · frame_sign`, the arc's travel about the axis read in the face's frame — with
/// the radial factor gone, because at a tangency the line runs along the circle, not across it.
/// A keyhole's cap: the box's wall tangent to the bore's rim at the corner where the box's
/// other wall stands on the rim. `Ok(None)` for any other pair — a turn that is not `0`, or two
/// departures a half turn apart ([`tangent_pole`]).
pub(crate) fn tangent_codirected_turn(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    a: &EdgeDir,
    b: &EdgeDir,
) -> Result<Option<i8>, BoolError> {
    let (line, arc, arc_second) = match (a, b) {
        (EdgeDir::Line { carrier, sense }, EdgeDir::Arc(arc)) => ((*carrier, *sense), arc, true),
        (EdgeDir::Arc(arc), EdgeDir::Line { carrier, sense }) => ((*carrier, *sense), arc, false),
        _ => return Ok(None),
    };
    if turn(jd, p, a, b)? != 0 {
        return Ok(None);
    }
    match tangent_travel_agrees(jd, p, line.0, line.1, arc) {
        Some(true) => {}
        Some(false) => return Ok(None),
        None => return Err(reject(RejectReason::WitnessNotRational)),
    }
    let up = if arc.axis_up { 1i8 } else { -1 };
    let way = if arc.ccw { 1i8 } else { -1 };
    let arc_left = up * way * jd.planes[p].frame_sign;
    Ok(Some(if arc_second { arc_left } else { -arc_left }))
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
    // **Stored-frame normals, not the canonical ones.** `d` is compared against `frame_sign`
    // below, which speaks the stored frame; `Judge::world_rat` speaks the canonical one and points
    // the other way on half the classes. See [`Judge::stored_world_rat`].
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

/// CCW cyclic order of the edges around one arrangement vertex on plane class `w`, **read from no
/// coordinate**. Each edge rides a line `w ∩ fp` and runs in direction `s·(n_w × n_fp)`; it is
/// given as `(fp, s)`. The signed turn between edges i and j is [`departure_turn`].
///
/// The turn sign is transitive only *within an open half-plane* (span < π), where it reproduces
/// the angle order exactly — the textbook Graham-scan fact. So: bucket every edge by its turn
/// against a reference `r = edges[0]` into the two open half-planes, plus the `0`/`π` pole where
/// the turn is 0 (collinear with `r`). The pole splits by same-`fp` opposite direction: same fp,
/// opposite `s` is angle π; otherwise angle 0. Each open bucket is then a real sort. Returned
/// indices are the CCW order starting at `r`.
///
/// **What it relies on, and who supplies it.** Two edges leaving one vertex with a `0` turn are on
/// one line — parallel plus a shared point — and then they have no cyclic order to read. Upstream
/// makes that impossible in two places: `merge_coincident` folds an edge traced twice (its own doc
/// names this bucket as the reason), and `Aliases::union_wall` folds two walls that carry one line
/// so they arrive under a single `fp`. Collinear same-direction overlap is resolved earlier
/// still, by `split_at_crossings`.
///
/// **That reliance is checked, not assumed.** `Aliases::record` only folds where four or more
/// planes meet at a point, and that every same-line wall pair meets such a point is not
/// established — so a `0` turn that survives to here is an honest `UnorderedEdges` reject rather
/// than an order picked arbitrarily. Measured: it fires nowhere in the suite, and disabling
/// `union_wall` makes the four-plane concurrency model reject, which is what says the reliance is
/// real and the check reaches it.
pub(crate) fn angular_order(
    jd: &Judge<'_, WorkingPlane>,
    w: usize,
    edges: &[EdgeDir],
) -> Result<Vec<usize>, BoolError> {
    // The comparator the winding reads too ([`departure_turn`]), which is what keeps the two
    // from drifting; the `0` it returns is *this* function's to interpret — here it is the `0`/π
    // pole.
    let cross = |i: usize, j: usize| departure_turn(jd, w, &edges[i], &edges[j]);
    let (mut zero, mut pos, mut pole, mut neg) = (vec![], vec![], vec![], vec![]);
    for i in 0..edges.len() {
        match cross(0, i)? {
            c if c > 0 => pos.push(i),
            c if c < 0 => neg.push(i),
            // Collinear with the reference: angle 0 (same ray) or π (opposite ray). A line and
            // an arc tangent at this vertex with the arc leaving the other way (a fillet's smooth
            // corner) are π apart too — the structural test cannot see it, the geometric
            // one can ([`tangent_pole`]) — and so are a line and a ruling that are
            // one line, leaving opposite ways ([`ruling_pole`]).
            _ if departure_pole(jd, w, &edges[i], &edges[0])? => pole.push(i),
            _ => zero.push(i),
        }
    }
    // **The π pole is a question about the direction's *representation*, not about the turn**
    // — "same wall class, opposite travel sign", and for arcs "same circle, opposite tangent".
    // Both live in [`antiparallel`]. Measured: the arc/arc arm never fires here,
    // because segments are numbered before arcs and every crossing has one, so `edges[0]` is always
    // a line; the two arcs at a node sit π apart and land in the two *open* buckets instead.
    //
    // **The order below exists only if no two edges share an angle, and that is an
    // assumption about *upstream*, not about this function.** Two edges leaving one vertex whose
    // turn is zero are on one line — parallel plus a shared point — so their cyclic order has no
    // answer, and `next` is built from exactly that order: a guess there is a wrong face, with
    // nothing to say so.
    //
    // Upstream is supposed to make it impossible, in two places, and both are load-bearing:
    // `merge_coincident` folds an edge that was traced twice (its own doc names this bucket as
    // the reason), and `Aliases::union_wall` folds two walls that carry one line. Measured, the
    // second fires 122 times over the coverage suite and this check never does — but **"measured
    // never" is not "cannot"**: `Aliases::record` only folds where four or more planes meet at a
    // point (`s.len() < 4` returns early), and that every same-line wall pair meets such a point
    // is *not* established. So the assumption is checked rather than trusted, and a failure is an
    // honest reject.
    //
    // The pairs are checked before sorting rather than inside the comparator: the sort does not
    // compare every pair, so a tie it happens not to evaluate would slip through — and that is
    // precisely the tie that reorders the result. At the degrees this sees (2-3 measured) a
    // bucket holds at most two entries, so the check costs the one comparison the sort would
    // have made anyway.
    for b in [&pos, &neg] {
        for x in 0..b.len() {
            for y in (x + 1)..b.len() {
                if cross(b[x], b[y])? == 0 {
                    return Err(reject(RejectReason::UnorderedEdges));
                }
            }
        }
    }
    // The reference's own bucket: anything else collinear with it that the pole test did not
    // claim is the same ray, which is an overlap `split_at_crossings` should have resolved.
    if zero.len() > 1 {
        return Err(reject(RejectReason::UnorderedEdges));
    }
    // Within an open half-plane, `a` precedes `b` (smaller angle) iff `d_a × d_b > 0`.
    // `turn` can fail rather than answer — an arc's side is `a + b√c` arithmetic — and the sort
    // stops at it ([`try_sort_by`]). A zero turn is unreachable after the check above.
    let by_turn = |v: &mut Vec<usize>| -> Result<(), BoolError> {
        try_sort_by(v, |&a, &b| cross(a, b).map(|c| 0.cmp(&c)))
    };
    by_turn(&mut pos)?;
    by_turn(&mut neg)?;
    let mut out = zero;
    out.extend(pos);
    out.extend(pole);
    out.extend(neg);
    Ok(out)
}
