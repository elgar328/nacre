use super::*;
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
/// ★★ **`None` where the node is not three planes.** A [`NodeKind::Pierce`] *is* a point, but its
/// coordinates are quadratic-irrational and `orient3d` is the plane-triple judge — so this says
/// "not mine to answer" rather than guessing. Callers turn that into their own vocabulary (the
/// tracer a [`crate::DeclineKind`], the ray caster a reject), which is why it is not a reject here.
///
/// ★★★★★ **That arm used to be unexercised, and the rung that took `plane_ring`'s checks away
/// fired it — with a wrong sign.** The scan road now walks rings whose corners a cylinder made,
/// and the first thing that came back was a pierce corner on the *opposite* side of its own face's
/// plane from its four plane-named neighbours (measured: `sides = [1, -1, -1, -1, -1, 1]`, the two
/// `1`s being the pierce nodes). The scan read those as crossings that are not there and named one
/// with a plane triple that never met, and `orient3d` answered `D = 0`. ★ I reported that panic as
/// a hole upstream in the naming; every ring name measured correct and the fault was here.
pub(crate) fn side_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    n: NodeId,
    q: usize,
) -> Option<i8> {
    match n.kind() {
        NodeKind::ThreePlane(t) => Some(jd.orient3d(t[0], t[1], t[2], q)),
        // ★★ **A pierce point's side of a plane is one `a + b√c` sign.** The point is
        // `line.base() + s·line.dir()`, the plane's coefficients are rational, and
        // `quad::plane_side` is that sign — the same predicate `ruling_side` reads. `None` is a
        // missing description (a class with no world name, a cylinder with no world statement),
        // never a shape this cannot answer.
        NodeKind::Pierce { cyl, .. } => {
            let (line, sv) = pierce_meet(jd, cyl, &cyls.get(cyl)?.def, n)?;
            let co = class_coeffs_rat(jd, q)?;
            let fix = outward_fix(jd, q)?;
            Some(
                fix * match nacre_exact::quad::plane_side(&co, &line, &sv) {
                    nacre_exact::Orient::Positive => 1,
                    nacre_exact::Orient::Negative => -1,
                    nacre_exact::Orient::Zero => 0,
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
/// node* (a road) — and the four callers do different things with each.
///
/// ★ There used to be a fourth outcome, `CurvedDeparture`: a ring leaving the meet along a curved
/// edge between two on-meet nodes, whose side the walk could not name. The caller names it now
/// ([`EdgeMeet::Departs`]), and the walk reads the departure as one more off-line entry.
pub(crate) enum RingWalk {
    /// Where the ring meets the line, in ring order from the first off-`q` node.
    Met(Vec<Feature>),
    /// Every node **and every edge** lies on `q`: a ring in the plane has no flanks to be decided
    /// by.
    ///
    /// ★★ **The second half is new and it is not pedantry.** This used to fire on "every node lies
    /// on `q`", which a ring with a curved edge satisfies while *leaving* the plane — and one
    /// consumer ([`every_ray`]) skips such a ring entirely, dropping its crossings from a parity
    /// count. A ring whose nodes are all on `q` while an edge departs is `Met` instead: the
    /// departure is an off-line entry of its own side ([`EdgeMeet::Departs`]).
    AllOn,
    /// A node whose side this walk cannot answer.
    ///
    /// ★★ **It used to mean "a [`NodeKind::Pierce`]", and it does not any more.** [`side_of`]'s
    /// pierce arm answers, so a cylinder's corner is read like any other; what is left here is
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
    /// a graze over ground it does not bound. The run is cut at each departure now, and the
    /// departure's own side flanks the pieces ([`EdgeMeet::Departs`]), so the sentence above is
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
    /// ★ A stretch that begins where the ring *returned* to the meet is preceded by the departure
    /// itself, and that side is the departure's σ ([`EdgeMeet::Departs`]) — so there is always a
    /// neighbour to read, and the number is never a plausible stand-in.
    Run {
        first: usize,
        len: usize,
        flanks_differ: bool,
        flank: i8,
    },
}

/// Whether ring edge `i` lies on what `q` cuts here, or leaves it — see [`ring_against_plane`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum EdgeMeet {
    /// The edge lies on the meet: two on-`q` nodes joined by a straight edge, or — on the circle
    /// road — an arc carried by `q` itself.
    On,
    /// The edge leaves the meet between two on-`q` nodes (a curved edge) and stays on this side of
    /// `q`, in [`side_of`]'s frame, until it returns — a circle meets a plane in two points, so
    /// the whole excursion lies on one side.
    Departs(i8),
}

/// Read a ring against one plane: where it meets the line, and whether it crosses or only touches.
///
/// ★ **One walk, four consumers.** `arrangement::trace_transversal_face` clips a face's ring
/// against a cut plane, `arrangement::cycle_on_class` reads a lateral's cycle against a circle,
/// [`every_ray`] casts a parity ray along `P ∩ Q_a` and [`segment_meets_face`] alternates a
/// segment against a face; all must answer the same question first — *does the boundary cross
/// this line here?* — and a node sitting **on** the line is the only hard part of it. The tracer
/// had the rule (look at the node's two off-line neighbours: opposite sides is a crossing, equal
/// sides a touch) inlined in its scan, entangled with naming, alias recording and decline kinds;
/// the ray caster had no rule at all and threw such a candidate away. Resilience that lives in one
/// consumer is resilience the other does not have — the same shape the nesting road's retry was in
/// before `nesting::cell_inside` became its one loop.
///
/// What each consumer does *with* a feature stays its own: the tracer turns it into a **named
/// point** (four-plane aliases, `DeclineKind`, occupancy), the ray caster into one bit ("ahead of
/// `v`?"). That is where the sharing stops.
///
/// `AllOn` when every node **and edge** lies on `q` — a ring in the plane has no flanks to be
/// decided by.
///
/// **Features come out in ring order from the first off-`q` entry.** That is the order the tracer's
/// scan produced them in, and its naming step records aliases into a union-find as it goes, so the
/// order is contract, not incident.
///
/// ★★ **`on_meet(i)` answers "does ring edge `i` lie on what `q` cuts here, and if not, which side
/// does it leave to?"** — edge `i` runs from `nodes[i]` to `nodes[i + 1]`, and it is asked only
/// where both ends are on `q`. Two on-`q` nodes do **not** settle it: two points fix a straight
/// line, so a straight edge between them is on it, and a **curved** one leaves and comes back.
/// What counts as "the meet" is the caller's, because it differs by road — a plane class cuts a
/// planar face in a **line** (so a straight edge is `On` and an arc `Departs` to the side its
/// tangent points, [`arc_departure_side`]) and a lateral face in a **circle** (so an arc carried
/// by that very class is `On`). `None` is "no exact description" and answers [`RingWalk::Unnameable`].
///
/// ★★★★★ **A departing edge is read as one more off-line entry, of the departure's side.**
/// The ring's sign sequence is then nodes and departures alike, and the two rules the walk has
/// always had — a sign change between neighbours is a crossing, a maximal run of zeros is an
/// on-line interval flanked by its neighbours — apply unchanged: a run is cut where the ring
/// leaves, each piece's flanks are the departures beside it, and a ring whose every node is on
/// `q` (a half-disk cap: its chord and its arc) is a run flanked by its own arc on both sides.
/// Crossings are still only ever between two *nodes*: a departure sits between two on-`q` nodes,
/// so it is never adjacent to an off-`q` one.
pub(crate) fn ring_against_plane(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    nodes: &[NodeId],
    q: usize,
    on_meet: impl Fn(usize) -> Option<EdgeMeet>,
) -> RingWalk {
    let n = nodes.len();
    let Some(side) = (0..n)
        .map(|i| side_of(jd, cyls, nodes[i], q))
        .collect::<Option<Vec<i8>>>()
    else {
        return RingWalk::Unnameable;
    };
    // The sign sequence: every node, and after node `i` its edge where that edge departs the
    // meet between two on-`q` nodes — `(side, Some(node))` or `(σ, None)`.
    let mut seq: Vec<(i8, Option<usize>)> = Vec::with_capacity(n);
    for i in 0..n {
        seq.push((side[i], Some(i)));
        if side[i] == 0 && side[(i + 1) % n] == 0 {
            match on_meet(i) {
                None => return RingWalk::Unnameable,
                Some(EdgeMeet::On) => {}
                Some(EdgeMeet::Departs(s)) => seq.push((s, None)),
            }
        }
    }
    let m = seq.len();
    let Some(start) = seq.iter().position(|e| e.0 != 0) else {
        return RingWalk::AllOn;
    };
    let mut out = Vec::new();
    let mut j = 0;
    while j < m {
        let i = (start + j) % m;
        if seq[i].0 != 0 {
            let ni = (i + 1) % m;
            if seq[ni].0 != 0 && seq[ni].0 != seq[i].0 {
                let (Some(a), Some(_)) = (seq[i].1, seq[ni].1) else {
                    unreachable!("a departure sits between two on-line nodes")
                };
                out.push(Feature::Crossing {
                    edge: a,
                    from: seq[i].0,
                });
            }
            j += 1;
        } else {
            // A maximal run of on-line vertices, ended by an off-line node or by a departure. Two
            // is the common case, but a vertex whose name had to be taken from its touching
            // planes (`loop_triples`) stays in the ring even when the loop runs straight through
            // it, so a run can be longer — and a run between two departures can be a single node.
            let Some(first) = seq[i].1 else {
                unreachable!("an on-line entry is a node")
            };
            let mut len = 0usize;
            while j < m && seq[(start + j) % m].0 == 0 {
                len += 1;
                j += 1;
            }
            let before = seq[(i + m - 1) % m].0;
            let after = seq[(start + j) % m].0;
            out.push(Feature::Run {
                first,
                len,
                flanks_differ: before != after,
                flank: before,
            });
        }
    }
    RingWalk::Met(out)
}

/// Is the implicit point `v` inside the simple ring `ring`, both on face plane `p`?
///
/// **A ray, cast along a line we already have.** Every ring edge lies on `P ∩ R`, and `v`
/// lies on `P ∩ Q_a` for either of its own two planes. Those two lines meet at
/// `X = {P, Q_a, R}`, which is *itself* a three-plane point — so "is `X` inside the edge"
/// and "is `X` ahead of `v`" are both [`order_along`], the comparator already built
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
