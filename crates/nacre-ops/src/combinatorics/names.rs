use super::*;
/// **A vertex's identity** — the sorted plane triple that names it.
///
/// `Eq`/`Hash` give identity dedup so an A-piece and a B-piece that meet at a seam node share one
/// result vertex/edge; `Ord` gives the deterministic node order replay needs — and it is the bare
/// triple's lexicographic order, so every "smallest name wins" rule reads unchanged.
///
/// This was `boolean::Node`, spoken only by the assembler. It lives here because the arrangement
/// names the same vertices, and the variant is spelled like [`nacre_topo::Vertex::ThreePlane`]
/// so the arrangement, the assembler and the topology store call the thing by one name.
///
/// Read it with `match`, never `let`-`else`: a new variant lights up the first and falls silently
/// into the second — a defect this repository has already had.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub(crate) enum NodeId {
    ThreePlane([usize; 3]), // sorted triple (key into the seam map)
    /// Where two plane classes' meet line crosses a cylinder's lateral surface — the point
    /// [`nacre_topo::Vertex::Pierce`] names, and the point the next rung splits a circle
    /// into arcs at.
    ///
    /// ★ **Two index spaces in one name.** `planes` are plane-class indices and `cyl` is a
    /// **cylinder**-class index (`crate::planes::ClassIx::Cyl`'s payload); they are separate
    /// numberings and a value from one is meaningless in the other. `reuse::canonical` already
    /// carries both kinds in one key, so the precedent is the file's, not this type's.
    Pierce {
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
    /// [`NodeId::pierce`], whose pair carries a root that a re-sort has to restate.
    ///
    /// It does **not** check for a collapsed triple. Two names being equal is a real condition
    /// with *different answers at different callers* — `arrangement::plane_ring` declines
    /// (`CollapsedTriple`), [`loop_triples`] falls back to naming the vertex from every plane
    /// touching it — so making the constructor fallible would copy that fork to all eight minting
    /// sites.
    pub(crate) fn three_planes(t: Canon3) -> NodeId {
        NodeId::ThreePlane(t.planes())
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
    pub(crate) fn pierce(
        first: usize,
        second: usize,
        cyl: usize,
        root: nacre_topo::QuadRoot,
    ) -> NodeId {
        let (planes, root) = nacre_topo::QuadRoot::canonical([first, second], root);
        NodeId::Pierce { planes, cyl, root }
    }
}

/// **The one door out of the identity and into the machinery that assumes every vertex has a
/// three-plane name** — the ray casts, the ring walks, the wall-and-handle derivations, and the
/// comparison keys. They take `[usize; 3]`, and rightly: an in-flight probe point like
/// `point_in_component`'s `{a, b, q}` is three planes without being any vertex of the arrangement,
/// so a name is the wrong type for their parameter.
///
/// **`None` is a pierce point**, and the answer for every caller behind this door is the same one:
/// it has no three-plane name and none of those paths has another to give. What differs is what
/// each does about it, and that splits in two — see [`three_plane_probes`] for the split.
///
/// Sites that *dispatch* instead of declining stay out of it deliberately — [`node_coords_rat`],
/// [`pierce_point`] and [`loop_winding`]'s lexicographic scan, each with its own answer for the
/// other variant. ★ And one is **not** in the ring-and-segment world at all: `reuse::canonical`
/// builds a comparison key, and its answer is that the key's vessel widens (`CanonNode`), not that
/// the question is refused.
///
/// ★★ **The gate that keeps this honest**, and what it exempts — this module and one file:
///
/// ```text
/// rg 'NodeId::(ThreePlane|Pierce)' crates/ \
///   -g '!**/combinatorics/**' -g '!**/reuse.rs' | grep -vE ':\s*//'
/// ```
///
/// It must be empty. Spelling a variant anywhere else means a site went around the door instead of
/// answering — the same gate once found **fifteen** of those with the whole suite
/// already green, and it costs nothing and does not break on a rename.
///
/// ★★★ **[`pierce_name`] is this door's twin.** This door answers only the three-plane half;
/// the sites in `arrangement` that need the `Pierce` half — `arc_split_witness`' `separates`
/// closure, two inside `split_circles`, and the seam table — go through the twin rather than
/// opening the variant directly.
///
/// ★ Not automated, and that is how it rots: three sites can arrive with nothing running the
/// check, and a fourth be nearly added with the suite green. `tests/probes/rotation_sweep.rs`'
/// `side_of` guard is the precedent for making a source scan a test.
pub(crate) fn three_plane_name(n: NodeId) -> Option<[usize; 3]> {
    match n {
        NodeId::ThreePlane(t) => Some(t),
        NodeId::Pierce { .. } => None,
    }
}

/// [`three_plane_name`]'s twin — the payload of a **pierce** name, `None` for a three-plane one.
///
/// The three questions its consumers ask are all payload reads: does this root *separate* (a
/// tangency's `Double` does not), which cylinder is the point on, where does it sort along the
/// meet line. Spelling the variant at those sites instead is what the gate above forbids — the
/// door is total over the enum, so a third variant becomes a compile error here rather than a
/// silent fall-through at four call sites.
pub(crate) fn pierce_name(n: NodeId) -> Option<([usize; 2], usize, nacre_topo::QuadRoot)> {
    match n {
        NodeId::ThreePlane(_) => None,
        NodeId::Pierce { planes, cyl, root } => Some((planes, cyl, root)),
    }
}

/// **The names of a list of *candidates*, pierce points dropped.**
///
/// ★★★ **This is the licence, and its name is where the licence is stated.** Behind the door there
/// are two shapes and only one of them may lose a member:
///
/// - a **probe list** may — its consumers try each member until one decides, and an exhausted list
///   is already a named decline (`nesting::cell_inside` answers `NoClearRay`, `boolean`'s `first_deciding`
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
/// [`NodeId::Pierce`] the edge already carries in its endpoint list. So the pin says **which kind**
/// and the name is read from beside it, rather than a second copy living here.
///
/// ★★ The two arms are two *orders*, not two spellings of one: [`order_along`] reads a class
/// through `orient3d` × `dir_sign` (integer predicates), and a pierce point through the
/// `a + b√c` tower. Naming the kind is what makes the second reachable at all — a `usize` had
/// nowhere to say "not a plane".
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum EndPin {
    /// The third plane class: the point is `P ∩ wall ∩ this`.
    Class(usize),
    /// A cylinder crossing: the point is the [`NodeId::Pierce`] this endpoint is named by.
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
        /// ★★ An arc split cuts a segment at pierce points, and `order_along` speaks three-plane
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

/// **A ring edge's carrier** — the type half of "a line is unordered, a circle is ordered".
///
/// A plane-carried edge rides one wall class, as `Ring.walls` always said. An arc rides a
/// cylinder, and for it the ring additionally remembers **which way around the axis this edge
/// runs**: `ccw` restates `ClassEdges::edge_at`'s own convention (*"`MergedArc::end` runs
/// counter-clockwise about the axis"* — the even half-edge travels that way, its twin the other),
/// carried rather than re-derived. That bit is what will let `edge_for` tell the two
/// complementary arcs between one pair of pierce vertices apart.
///
/// ★ Replacing the `usize::MAX` sentinel with a variant also kills a recorded hazard for free:
/// `dissolve_straight_angles` folds on wall *equality*, and two arcs of different circles — or
/// of one circle in different directions — now compare unequal instead of `MAX == MAX`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Wall {
    Plane(usize),
    Arc {
        cyl: usize,
        ccw: bool,
    },
    /// A ruling piece: straight on the lateral, so like a line its ends
    /// order it — but its carrier is the cylinder, and `(cyl, side)` names which of the two
    /// parallel rulings ([`combinatorics::RulingCarrier::side`]; `0` is a **tangent** wall's
    /// single ruling). `up` restates
    /// `ClassEdges::edge_at`'s convention (`MergedRuling::end` ascends the axis; the even
    /// half-edge travels up, its twin down), carried like `Arc::ccw`.
    Ruling {
        cyl: usize,
        side: i8,
        up: bool,
    },
}

impl Wall {
    /// The same carrier walked the other way: a plane wall is direction-blind, the same
    /// arc walked back runs the other way about the axis, the same ruling piece descends.
    /// The *complementary* arc between the same two nodes keeps its flag instead — which
    /// is what lets a keyed lookup tell "this edge reversed" from "the other arc".
    pub(crate) fn reversed(self) -> Self {
        match self {
            Wall::Plane(p) => Wall::Plane(p),
            Wall::Arc { cyl, ccw } => Wall::Arc { cyl, ccw: !ccw },
            Wall::Ruling { cyl, side, up } => Wall::Ruling { cyl, side, up: !up },
        }
    }
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

/// The **ruling** a straight lateral edge rides: a wall plane parallel
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
    /// Which of the two parallel rulings, by the sign convention above — or `0`, the single
    /// ruling of a **tangent** wall: an identity key like the other two values, never
    /// a sign to multiply by (the sign consumers assert it away).
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
/// ★ **Test-only**: the tracer's crossed-edge
/// wall reads the wall the producer carries (`NamedRing`) instead, so no production path
/// derives ring geometry from names. Kept for hand-built test rings, whose vertices
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
                node: NodeId::three_planes(Canon3::three(a)),
                to: NodeId::three_planes(Canon3::three(b)),
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
/// ★★★★★ **It was written five times before it was written once.** the ring-edge derivation (now `Ring::edges` in `boolean`), the ray
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

/// **A three-plane name that went through the rule** — the only thing
/// [`NodeId::three_planes`] accepts.
///
/// Two ways in, both here: [`canonical_triple`] for a *set* of planes known to pass through a
/// point (the rule picks), and [`Canon3::three`] for a *construction* that produces exactly three
/// — a line pinned by a plane, a solid's corner of three faces — where there is nothing to pick.
/// What the type cannot say is that a caller passed **every** plane it knew: three of four is the
/// defect in a new coat, and the operand-vertex audit (`operand_vertex_audit`) is what measures
/// that, corpus-wide. What it does say is that no site spells the choice a seventh time.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub(crate) struct Canon3([usize; 3]);

impl Canon3 {
    /// Exactly three plane classes through the point, **by construction** — sorted here, so the
    /// order given does not matter. Not for a set the caller cut down to three: that is
    /// [`canonical_triple`]'s question.
    pub(crate) fn three(mut t: [usize; 3]) -> Canon3 {
        t.sort_unstable();
        Canon3(t)
    }
    /// The three classes, ascending.
    pub(crate) fn planes(self) -> [usize; 3] {
        self.0
    }
}

/// **The one rule that names a point from the planes through it.**
///
/// `s` is every plane class known to pass through one point, sorted and deduplicated. The name is
/// the lexicographically first triple of `s` that is **independent** — three planes sharing a line
/// name no point — and `None` when fewer than three classes are given.
///
/// Three classes need no judgement: they are the only candidate, and every producer that meets a
/// three-plane point derives the same three (the ordinary vertex; this arm leaves its name what it
/// always was, and asks nothing of the judge). Four or more is a concurrency, and the choice among
/// the candidates is what has to be **one rule**: the alias table's representative
/// (`arrangement::Aliases`) is the minimum of its union-find, i.e. this very triple, so a vertex
/// an operand names here and a point the arrangement discovers fold onto the same name.
///
/// ★ **Why one function.** The rule has many users —
/// the arrangement's run-vertex discovery, its alias representative, the result vertex's
/// definition, a line's wall family, an edge's end triples. Spelled once per user, with
/// the operand's own vertices keeping a rule of their own (a
/// triple per face loop), a concurrency arrives under four names, one of them dependent,
/// and the judge reads the dependent one as lying on every class.
///
/// ★ **Completeness** — why an operand's set plus the arrangement's `{wc} ∪ t` records make one
/// component: with `|F| = 4` planes through a point every record *is* `F`, so the union-find has
/// one component whose minimum is this triple. Five or more is not in the corpus
/// (`concurrent_vertices_are_four_planes_and_the_trace_sees_all_of_them`) and is the recorded
/// stop condition.
pub(crate) fn canonical_triple(jd: &Judge<'_, WorkingPlane>, s: &[usize]) -> Option<Canon3> {
    debug_assert!(
        s.windows(2).all(|w| w[0] < w[1]),
        "sorted, deduplicated: {s:?}"
    );
    match *s {
        [] | [_] | [_, _] => None,
        [a, b, c] => Some(Canon3::three([a, b, c])),
        _ => {
            for i in 0..s.len() {
                for j in (i + 1)..s.len() {
                    for k in (j + 1)..s.len() {
                        if jd.plane_pair_dir_sign(s[i], s[j], s[k]) != 0 {
                            return Some(Canon3::three([s[i], s[j], s[k]]));
                        }
                    }
                }
            }
            None
        }
    }
}

/// **What pins a named point on the line `a ∩ b`** — the one rule, read from the name.
///
/// A three-plane name is pinned by whichever of its planes cuts the line ([`pin_on_line`]). A
/// pierce name whose pair *is* the line's is pinned by its cylinder ([`EndPin::Cylinder`] — the
/// quadric's root along that very line); a pierce name with another pair is a point the line
/// passes through by coincidence, pinned by whichever of that pair's planes cuts the line.
///
/// ★ Why it exists: a segment's endpoint travels as a (name, pin) pair, and the alias table can
/// fold the name onto a representative of another **variant** — a tangent corner's `Pierce` onto
/// the `ThreePlane` of its planes with the class through it. The pin is a fact about the
/// representative, so it is derived again from it, here, rather than carried across the fold.
pub(crate) fn pin_for(
    jd: &Judge<'_, WorkingPlane>,
    a: usize,
    b: usize,
    name: NodeId,
) -> Option<EndPin> {
    match name {
        NodeId::ThreePlane(t) => pin_on_line(jd, a, b, t).map(EndPin::Class),
        NodeId::Pierce { planes, .. } => {
            let mut line = [a, b];
            line.sort_unstable();
            if planes == line {
                return Some(EndPin::Cylinder);
            }
            planes
                .iter()
                .copied()
                .find(|&c| c != a && c != b && jd.plane_pair_dir_sign(a, b, c) != 0)
                .map(EndPin::Class)
        }
    }
}
