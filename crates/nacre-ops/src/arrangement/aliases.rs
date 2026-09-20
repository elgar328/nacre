use super::*;
/// Sort a triple and reject it if two of its planes coincide.
///
/// The triples already carry dense plane ids — `loop_triples` maps face indices through `plane_ix`
/// at the source — so there is nothing to canonize here; the point is the collapse check. `None`
/// when two names are equal: such a triple defines no point (`three_planes` answers `None`, and the
/// exact predicates need `D ≠ 0`), so the face is declined rather than fed a degenerate meet. A
/// face→plane map still collapses — two faces of one solid meeting a vertex on one plane — so this
/// guard outlives the old face/plane ambiguity it was born with.
/// The reject a declined face reports to the caller.
///
/// Most kinds *are* the answer — the tracer says what it could not do, and `face` says where.
/// [`DeclineKind::FourPlane`] is the exception: the naming failure is a symptom, and reporting it
/// as one would hide the substrate limit that caused it, so it is raised as the cause. The face
/// handle is dropped there because [`RejectReason::FourPlane`] carries no payload; per-face detail
/// remains in the class audit, which is where `TraceDeclined`'s own docs put it.
///
/// **One function for both consumers.** The boolean's error and the audit's `failed_at` must agree
/// — two copies of this mapping would let them drift.
pub(super) fn decline_to_reject(kind: DeclineKind, face: Option<Handle<Face>>) -> RejectReason {
    match kind {
        DeclineKind::FourPlane => RejectReason::FourPlane,
        kind => RejectReason::TraceDeclined { kind, face },
    }
}

/// **The alias table learns from the operands first**: every concurrency an operand's
/// own topology knows — a vertex with four or more incident plane classes, carried by its rings as
/// `NamedRing::concurrencies` — is recorded before any class is traced. The representative the
/// ring already named the vertex by (`canonical_triple`) is then the representative the
/// arrangement's own discoveries (`{wc} ∪ t` at a run vertex, a wall family's line) fold onto.
pub(super) fn seed_from_operands(
    aliases: &mut Aliases,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    trace_in: &combinatorics::TraceInput,
) {
    let mut corners: Vec<NodeId> = Vec::new();
    for side in &trace_in.faces {
        for (_, loops) in side {
            let rings = loops
                .outer
                .iter()
                .chain(loops.holes.iter().flatten())
                .chain(loops.cycles.iter().flatten().map(|(_, r)| r));
            for lr in rings {
                if let Some(nr) = lr.poly() {
                    for s in &nr.concurrencies {
                        aliases.record(jd, s);
                    }
                    corners.extend(
                        nr.triples
                            .iter()
                            .copied()
                            .filter(|&n| combinatorics::pierce_name(n).is_some()),
                    );
                }
            }
        }
    }
    // ★ **The operands' cylinder corners, against every class.** A pierce corner
    // (`Pierce{[p0, p1], cyl, root}` — a fillet's tangent corner, a boss's foot) is a point some
    // *other* class may pass through: the gusset's side plane through the fillet axis contains
    // the tangent ruling, and so the corner. Then the point has three names — its own, the
    // three-plane `[p0, p1, wc]`, and the class's crossing of that ruling — and
    // [`Aliases::record_on_cylinder`] joins them. Asked here, once, of the operands' own
    // topology and the class table, so every round and every class trace starts knowing it:
    // asked during a trace instead (where `third_on_l` meets the corner on `wc`), a class traced
    // in the same round could refuse the coincidence before the round that learnt it — a round
    // that declines returns, and there is no next. `side_of` is exact (the quad tower), and the
    // question is corners × classes, both small.
    corners.sort_unstable();
    corners.dedup();
    for &corner in &corners {
        let Some((planes, _, _)) = combinatorics::pierce_name(corner) else {
            continue;
        };
        for c in 0..jd.planes.len() {
            if planes.contains(&c) {
                continue;
            }
            if combinatorics::side_of(jd, cyls, corner, c) == Some(0) {
                aliases.record_on_cylinder(jd, cyls, corner, c);
            }
        }
    }
}

/// The names that turned out to denote **one** feature, learned while tracing.
///
/// Both kinds come from the same discovery. When a producer learns the full set `S` of planes
/// through a point:
///
/// - **A point** with `|S| > 3` has one valid name per 3-subset, so the tables would enter it once
///   per name. [`Aliases::point`] folds them onto the lexicographically first subset that actually
///   names a point (three planes sharing a line name none, so such a subset cannot be the winner).
/// - **A line** shows up as a 3-subset of `S` whose planes share a line rather than meeting at the
///   point. On each of those three classes the other two are *the same wall* — the engine calls a
///   line by the plane it rides besides the cut plane, so that is two names for one line, and
///   [`Aliases::wall`] folds them onto the smallest.
///
/// Neither fold is optional if the other happens: `merge_coincident` keys an edge by
/// `(wall, endpoints)`, so a duplicate edge merges only when **both** its wall and its endpoint
/// names agree. That is why the two live in one table and are applied together.
///
/// ★ The table has **two sources**. The operands seed it before any class is traced
/// ([`seed_from_operands`] — every vertex with four or more incident plane classes, which the
/// operand's own topology knows), and the tracer adds what it discovers (`{wc} ∪ t` at a run
/// vertex, a wall family's line). Among plane names the representative is
/// [`combinatorics::canonical_triple`]'s answer, which is also the name the operand's ring
/// already gave the point; with four planes through a point every record is the whole set, so the
/// two sources land in one component.
///
/// ★ A point of the table need not be a plane name at all. The seed also learns every
/// **pierce corner** an operand carries that some further class passes through
/// ([`Aliases::record_on_cylinder`]), so one component can hold a corner, a three-plane name and
/// a class's ruling crossing; the representative among those is a **pierce** name
/// ([`Aliases::rep_rank`] — a point on a cylinder is represented on the cylinder, because the
/// lateral chart reads the cylinder off the name).
// `Clone` so a round can hand every class the table as it stood when the round began, and
// merge their discoveries afterwards — see `trace_result_faces`. In the ordinary model the
// maps are empty, so the copy costs nothing.
#[derive(Default, Debug, Clone)]
pub(crate) struct Aliases {
    /// Union-find over vertex names.
    point: HashMap<NodeId, NodeId>,
    /// Union-find over `(class, wall)` — walls of one class that carry the same line.
    wall: HashMap<(usize, usize), usize>,
}

impl Aliases {
    /// Record that every plane in `s` (sorted, deduped) passes through one point.
    pub(super) fn record(&mut self, jd: &Judge<'_, WorkingPlane>, s: &[usize]) {
        if s.len() < 4 {
            return; // three planes meeting at a point is the ordinary case, and names nothing new
        }
        // ★ The representative is [`combinatorics::canonical_triple`]'s answer — the one
        // rule every producer of a point's name calls — so a name an operand's ring already gave
        // the point is the representative it folds onto here.
        let rep = combinatorics::canonical_triple(jd, s).map(NodeId::three_planes);
        for i in 0..s.len() {
            for j in (i + 1)..s.len() {
                for k in (j + 1)..s.len() {
                    let t = [s[i], s[j], s[k]];
                    if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                        // Shares a line: names no point, and tells us two walls are one line.
                        self.union_wall(t[0], t[1], t[2]);
                        self.union_wall(t[1], t[0], t[2]);
                        self.union_wall(t[2], t[0], t[1]);
                    } else if let Some(rep) = rep {
                        self.union_point(rep, NodeId::three_planes(Canon3::three(t)));
                    }
                }
            }
        }
    }

    /// **A pierce corner found on a further plane class** — the cylinder twin of
    /// [`Aliases::record`]. The corner's planes `p0, p1` and the class `wc` all pass through the
    /// point, and so does the cylinder; every name that set can produce denotes it: the corner's
    /// own, the three-plane name of `{p0, p1, wc}` (when independent), and `wc`'s ruling crossing
    /// with the corner's cap at the root that lies on the corner's side of `wc`. All of them are
    /// folded here, in one place, so the ruling sweep and the plane roads — which mint the names
    /// on their own — never have to decide "same point" themselves.
    ///
    /// The crossing's side is the corner's own: [`ruling_side`] of the corner's meet against
    /// `wc`'s stored normal, the predicate the sweep uses to tell `wc`'s two rulings apart.
    pub(crate) fn record_on_cylinder(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        corner: NodeId,
        wc: usize,
    ) {
        let Some((planes, cyl, _)) = combinatorics::pierce_name(corner) else {
            return;
        };
        let mut s = vec![planes[0], planes[1], wc];
        s.sort_unstable();
        s.dedup();
        if s.len() < 3 {
            return; // `wc` is one of the corner's own planes: the ordinary corner, nothing to fold
        }
        if let Some(t) = combinatorics::canonical_triple(jd, &s) {
            self.union_point(corner, NodeId::three_planes(t));
        }
        let Some(wcy) = cyls.get(cyl) else { return };
        let def = &wcy.def;
        let Some(w) = combinatorics::class_coeffs_rat(jd, wc) else {
            return;
        };
        let Some(meet) = combinatorics::pierce_meet(jd, cyl, def, corner) else {
            return;
        };
        let Some(side) = ruling_side(&w, def, (&meet.0, &meet.1)) else {
            return; // the corner sits on `wc`'s axis plane's own line: no ruling to name
        };
        for fc in planes {
            if let Ok(id) = crossing_on_ruling(jd, def, fc, wc, cyl, side) {
                self.union_point(corner, id);
            }
        }
    }

    fn find_point(&self, t: NodeId) -> NodeId {
        let mut x = t;
        while let Some(&p) = self.point.get(&x) {
            if p == x {
                break;
            }
            x = p;
        }
        x
    }

    /// **The representative of a point's names** — the one every key uses (`canon_point`).
    ///
    /// Among three-plane names it is the least, which is [`combinatorics::canonical_triple`]'s
    /// own answer (the operand's name and the arrangement's discoveries meet there).
    /// ★ **A point on a cylinder is represented on the cylinder.** A `Pierce` name
    /// locates the point exactly *and* says which cylinder it lies on — the lateral chart's
    /// geometry (θ about the axis, `pierce_meet`) reads that from the name — while a three-plane
    /// name of the same point (a class through a tangent corner: `ThreePlane([cap, t, wc])`) is
    /// a key only. So a pierce name outranks a three-plane name; among pierce names the least.
    /// Either way the representative is a function of the class alone.
    fn rep_rank(n: NodeId) -> (u8, NodeId) {
        (u8::from(matches!(n, NodeId::ThreePlane(_))), n)
    }

    fn union_point(&mut self, a: NodeId, b: NodeId) {
        let (ra, rb) = (self.find_point(a), self.find_point(b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if Self::rep_rank(ra) < Self::rep_rank(rb) {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.point.insert(hi, lo);
        self.point.entry(lo).or_insert(lo);
    }

    fn find_wall(&self, class: usize, w: usize) -> usize {
        let mut x = w;
        while let Some(&p) = self.wall.get(&(class, x)) {
            if p == x {
                break;
            }
            x = p;
        }
        x
    }

    fn union_wall(&mut self, class: usize, a: usize, b: usize) {
        let (ra, rb) = (self.find_wall(class, a), self.find_wall(class, b));
        if ra == rb {
            return;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.wall.insert((class, hi), lo);
        self.wall.entry((class, lo)).or_insert(lo);
    }

    /// The identity to key a vertex on — itself when nothing was merged with it.
    pub(crate) fn canon_point(&self, t: NodeId) -> NodeId {
        if self.point.is_empty() {
            return t; // the overwhelmingly common case pays nothing
        }
        self.find_point(t)
    }

    /// The name to call a line by on `class` — `w` itself when no other wall carries it.
    pub(crate) fn canon_wall(&self, class: usize, w: usize) -> usize {
        if self.wall.is_empty() {
            return w;
        }
        self.find_wall(class, w)
    }

    /// The walls that carry **one line** with `w` on `class` — `w` alone when none does.
    ///
    /// The wall fold says two classes name one line; this reads that back out, because a point on
    /// that line lies on **every** plane in the family, and the point fold needs to hear about it
    /// (see the report in [`split_at_crossings`]).
    pub(super) fn wall_family(&self, class: usize, w: usize) -> Vec<usize> {
        if self.wall.is_empty() {
            return vec![w];
        }
        let rep = self.find_wall(class, w);
        let mut fam: Vec<usize> = self
            .wall
            .keys()
            .filter(|(c, _)| *c == class)
            .map(|&(_, x)| x)
            .filter(|&x| self.find_wall(class, x) == rep)
            .collect();
        fam.push(w);
        fam.sort_unstable();
        fam.dedup();
        fam
    }

    pub(super) fn absorb(&mut self, other: &Aliases) {
        for (&k, &v) in &other.point {
            self.union_point(k, v);
        }
        for (&(c, w), &v) in &other.wall {
            self.union_wall(c, w, v);
        }
    }

    /// How much has been learned — the fixed-point loop's measure.
    pub(super) fn len(&self) -> usize {
        self.point.len() + self.wall.len()
    }
}
