use super::*;
/// Sort a triple and reject it if two of its planes coincide.
///
/// The triples already carry dense plane ids — `loop_triples` maps face indices through `plane_ix`
/// at the source — so there is nothing to canonize here; the point is the collapse check. `None`
/// when two names are equal: such a triple defines no point (`three_planes` answers `None`, and the
/// exact predicates need `D ≠ 0`), so the face is declined rather than fed a degenerate meet. A
/// face→plane map still collapses — two faces of one solid meeting a vertex on one plane — so this
/// guard is needed.
/// The reject a declined face reports to the caller.
///
/// Most kinds *are* the answer — the tracer says what it could not do, and `face` says where.
/// Three are exceptions, raised as their cause. [`DeclineKind::FourPlane`]: the naming failure is a
/// symptom, and reporting it as one would hide the substrate limit that caused it.
/// [`DeclineKind::ObliqueSection`]: an oblique section the trace does not read is the oblique
/// coverage's end, which [`RejectReason::ObliqueCylinderCut`] names wherever else it ends; so is
/// [`DeclineKind::CylinderPair`], two cylinders meeting at a point a trace tied, which
/// [`RejectReason::CylinderPairContact`] names wherever else it is met. The face
/// handle is dropped there because [`RejectReason::FourPlane`] carries no payload; per-face detail
/// remains in the class audit, which is where `TraceDeclined`'s own docs put it.
///
/// **One function for both consumers.** The boolean's error and the audit's `failed_at` must agree
/// — two copies of this mapping would let them drift.
pub(super) fn decline_to_reject(kind: DeclineKind, face: Option<Handle<Face>>) -> RejectReason {
    match kind {
        DeclineKind::FourPlane => RejectReason::FourPlane,
        DeclineKind::ObliqueSection => RejectReason::ObliqueCylinderCut,
        DeclineKind::CylinderPair => RejectReason::CylinderPairContact,
        kind => RejectReason::TraceDeclined { kind, face },
    }
}

/// **The alias table learns from the operands first**: every concurrency an operand's own topology
/// knows — a vertex no cylinder made with four or more incident plane classes, carried by its rings
/// as `NamedRing::concurrencies` — is recorded before any class is traced. The representative the
/// ring already named the vertex by (`canonical_triple`) is then the representative the
/// arrangement's own discoveries (`{wc} ∪ t` at a run vertex, a wall family's line) fold onto.
pub(super) fn seed_from_operands(
    aliases: &mut Aliases,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    trace_in: &combinatorics::TraceInput,
) {
    aliases.nameless = (0..jd.planes.len())
        .map(|c| jd.world_rat(c).is_none())
        .collect();
    let mut corners: Vec<NodeId> = Vec::new();
    let mut three: Vec<NodeId> = Vec::new();
    for side in &trace_in.faces {
        for (_, loops) in side {
            let rings = loops
                .outer
                .iter()
                .chain(loops.holes.iter().flatten())
                .chain(loops.cycles.iter().flatten());
            for lr in rings {
                if let Some(nr) = lr.poly() {
                    for s in &nr.concurrencies {
                        let _ = aliases.record(jd, s);
                    }
                    corners.extend(
                        nr.triples
                            .iter()
                            .copied()
                            .filter(|&n| combinatorics::pierce_name(n).is_some()),
                    );
                    three.extend(
                        nr.triples
                            .iter()
                            .copied()
                            .filter(|&n| combinatorics::three_plane_name(n).is_some()),
                    );
                }
            }
        }
    }
    // **The operands' cylinder corners, against every class.** A pierce corner
    // (`Pierce{[p0, p1], cyl, root}` — a fillet's tangent corner, a boss's foot) is a point some
    // *other* class may pass through: the gusset's side plane through the fillet axis contains
    // the tangent ruling, and so the corner. Then the point has every name that set gives it —
    // its own, the three-plane `[p0, p1, wc]`, and each of its planes' crossings with the class —
    // and [`Aliases::record_on_cylinder`] joins them. Asked here, once, of the operands' own
    // topology and the class table, so every round and every class trace starts knowing it:
    // asked during a trace instead (where `third_on_l` meets the corner on `wc`), a class traced
    // in the same round could decline for want of it before the round that learnt it — a round
    // that declines returns, and there is no next. `side_of` is exact (the quad tower), and the
    // question is corners × classes, both small.
    // **And the dual: a line of two classes lying on a cylinder** — the gate's record, one
    // entry per such line, folded before any class is traced for the same reason as the corners
    // below.
    for (ci, cyl) in cyls.iter().enumerate() {
        for sr in &cyl.shared {
            aliases.record_shared_ruling(jd, cyls, ci, sr);
        }
        for gl in &cyl.grazed {
            aliases.record_graze_line(jd, cyls, ci, gl);
        }
    }
    // **The operands' three-plane corners on a cylinder.** A corner whose three planes meet exactly
    // on a cylinder's surface — a gusset's foot on a bore's rim — is one point with its pairs'
    // crossings there, folded before any class is traced for the same reason as the pierce corners
    // below. Asked of the surface, not of a lateral face: a name is a surface's (design «점의
    // 이름»), and whether a face is there is the result's question — the result vertex names only
    // the surfaces of the faces through it. A corner with no rational coordinates (a class with no
    // world statement among its planes) is not asked: no exact test puts it on the cylinder.
    three.sort_unstable();
    three.dedup();
    for &n in &three {
        let Some(t) = combinatorics::three_plane_name(n) else {
            continue;
        };
        // The coordinates are the cylinders' question: asked once, when the first one asks.
        let mut coords = None;
        for (ci, cy) in cyls.iter().enumerate() {
            let Some(p) = coords.get_or_insert_with(|| combinatorics::node_coords_rat(jd, n))
            else {
                break;
            };
            let def = &cy.def;
            if nacre_exact::cylinder_radial_side(p, &def.origin(), &def.dir(), def.r2())
                == nacre_exact::Orient::Zero
            {
                let _ = aliases.record_cylinder_point(jd, cyls, ci, n, &t);
            }
        }
    }
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
            #[expect(
                clippy::disallowed_methods,
                reason = "a point question: does this pierce corner lie on class `c`"
            )]
            let on_class = combinatorics::side_of(jd, cyls, corner, c) == Some(0);
            if on_class {
                let _ = aliases.record_on_cylinder(jd, cyls, corner, c);
            }
        }
    }
}

/// **A union the table refuses: the two components stand on two cylinders** — their point is
/// where two cylinders meet ([`Aliases::join_points`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct TwoCylinders;

impl From<TwoCylinders> for RejectReason {
    fn from(_: TwoCylinders) -> Self {
        RejectReason::CylinderPairContact
    }
}

/// The names that denote **one** feature, learned while tracing.
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
///   [`Aliases::canon_wall`] folds them onto one — the least by the line's order ([`Line::key`]),
///   read off the one table of lines ([`Lines`]), so every class of the line answers alike.
///
/// Neither fold is optional if the other happens: `merge_coincident` keys an edge by
/// `(wall, endpoints)`, so a duplicate edge merges only when **both** its wall and its endpoint
/// names agree. That is why the two live in one table and are applied together.
///
/// The table has **two sources**. The operands seed it before any class is traced
/// ([`seed_from_operands`] — every vertex with four or more incident plane classes, which the
/// operand's own topology knows), and the tracer adds what it discovers (`{wc} ∪ t` at a run
/// vertex, a wall family's line). Among plane names the representative is
/// [`combinatorics::canonical_triple`]'s answer, which is also the name the operand's ring
/// already gave the point; with four planes through a point every record is the whole set, so the
/// two sources land in one component.
///
/// A point of the table need not be a plane name at all. The seed also learns every
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
    /// The lines three or more plane classes hold ([`Lines`]).
    lines: Lines,
    /// Which classes have no world statement — a pierce name whose pair holds one is placed by an
    /// interval alone ([`Aliases::rep_rank`]). Learnt before any record ([`seed_from_operands`]).
    nameless: Vec<bool>,
    /// How many unions joined two components — the fixed-point loop's measure ([`Aliases::learned`]).
    joined: usize,
}

impl Aliases {
    /// Record that every plane in `s` (sorted, deduped) passes through one point. Every union is
    /// made; `Err` when one of them is refused ([`Aliases::join_points`]).
    pub(super) fn record(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        s: &[usize],
    ) -> Result<(), TwoCylinders> {
        // Three independent planes meeting at a point is the ordinary case and names nothing new;
        // three that share a line are a pencil, which is news for the line table.
        if s.len() < 3 || (s.len() == 3 && jd.plane_pair_dir_sign(s[0], s[1], s[2]) != 0) {
            return Ok(());
        }
        let mut out = Ok(());
        // The representative is [`combinatorics::canonical_triple`]'s answer — the one
        // rule every producer of a point's name calls — so a name an operand's ring already gave
        // the point is the representative it folds onto here.
        let rep = combinatorics::canonical_triple(jd, s).map(NodeId::three_planes);
        for i in 0..s.len() {
            for j in (i + 1)..s.len() {
                for k in (j + 1)..s.len() {
                    let t = [s[i], s[j], s[k]];
                    if jd.plane_pair_dir_sign(t[0], t[1], t[2]) == 0 {
                        // Shares a line: names no point, and tells us two walls are one line. A
                        // class of the set cuts the line where the set holds an independent three;
                        // a set that is the pencil alone asks every class.
                        if let Some(cut) = s
                            .iter()
                            .copied()
                            .chain(0..jd.planes.len())
                            .find(|&d| jd.plane_pair_dir_sign(t[0], t[1], d) != 0)
                        {
                            self.learn_line(jd, &t, cut);
                        }
                    } else if let Some(rep) = rep {
                        if let Err(e) =
                            self.union_point(rep, NodeId::three_planes(Canon3::three(t)))
                        {
                            out = Err(e);
                        }
                    }
                }
            }
        }
        out
    }

    /// **A pierce corner found on a further plane class** — the cylinder twin of
    /// [`Aliases::record`]. The corner's planes `p0, p1` and the class `wc` all pass through the
    /// point, and so does the cylinder; every name that set produces denotes it — the corner's own,
    /// the three-plane name of `{p0, p1, wc}` (when independent), and each of the corner's planes'
    /// crossing with `wc` at the root on the corner's other plane ([`Aliases::record_cylinder_point`]).
    /// All of them are folded here, before any trace, so the ruling sweep and the plane roads —
    /// which mint the names on their own — find them already one point.
    pub(crate) fn record_on_cylinder(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        corner: NodeId,
        wc: usize,
    ) -> Result<(), TwoCylinders> {
        let Some((planes, cyl, _)) = combinatorics::pierce_name(corner) else {
            return Ok(());
        };
        let mut s = vec![planes[0], planes[1], wc];
        s.sort_unstable();
        s.dedup();
        if s.len() < 3 {
            return Ok(()); // `wc` is one of the corner's own planes: the ordinary corner
        }
        self.record_cylinder_point(jd, cyls, cyl, corner, &s)
    }

    /// **Two names a round proved one point, recorded** — the record a tie inside the fixed-point
    /// loop makes instead of a refusal: the tie was an exact proof (an interval never orders two
    /// names equal), and whether two names are one point is the converged table's question, not a
    /// round's. Every plane of either name passes through the point, so the set is recorded as a
    /// concurrency ([`Aliases::record`] — a pencil among them is learnt as a line), a point on a
    /// cylinder gets every name that set gives it there ([`Aliases::record_cylinder_point`]), and
    /// the two names are joined. A union the table refuses — the point stands on two cylinders
    /// ([`Aliases::join_points`]) — is a point where two cylinders meet, which no road names
    /// (`CylinderPairContact`); two roots of one pair are two points, and a tie between them is the
    /// arithmetic's defect (`CoincidentNodes`).
    pub(super) fn record_same_point(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        a: NodeId,
        b: NodeId,
    ) -> Result<(), RejectReason> {
        // Read as the table reads them: each name restated onto its line's name.
        let (a, b) = (self.lines.restate(a), self.lines.restate(b));
        let (pa, pb) = (combinatorics::pierce_name(a), combinatorics::pierce_name(b));
        if let (Some((qa, ca, ra)), Some((qb, cb, rb))) = (pa, pb)
            && qa == qb
            && ca == cb
            && ra != rb
        {
            return Err(RejectReason::CoincidentNodes);
        }
        self.union_point(a, b)?;
        let mut s: Vec<usize> = combinatorics::name_planes(a);
        s.extend(combinatorics::name_planes(b));
        s.sort_unstable();
        s.dedup();
        self.record(jd, &s)?;
        if let Some((_, cyl, _)) = pa.or(pb) {
            let at = if pa.is_some() { a } else { b };
            self.record_cylinder_point(jd, cyls, cyl, at, &s)?;
        }
        Ok(())
    }

    /// [`Aliases::record_same_point`] for a trace, whose refusals are declines: two cylinders'
    /// crossings decline as [`DeclineKind::CylinderPair`], two roots of one pair as
    /// [`DeclineKind::CoincidentFeatures`].
    pub(super) fn record_same_point_traced(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        a: NodeId,
        b: NodeId,
    ) -> Result<(), DeclineKind> {
        self.record_same_point(jd, cyls, a, b).map_err(|r| match r {
            RejectReason::CylinderPairContact => DeclineKind::CylinderPair,
            _ => DeclineKind::CoincidentFeatures,
        })
    }

    /// **Every name of one point on cylinder `cyl`, folded** — `at` is one, `s` the plane classes
    /// through it ([`combinatorics::names_on_cylinder`], the one rule). Every seed of a point on a
    /// cylinder ends here. Every union is made; `Err` when one of them is refused — a name already
    /// standing on another cylinder ([`Aliases::join_points`]).
    fn record_cylinder_point(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        cyl: usize,
        at: NodeId,
        s: &[usize],
    ) -> Result<(), TwoCylinders> {
        let mut out = Ok(());
        for n in combinatorics::names_on_cylinder(jd, cyls, cyl, at, s) {
            if let Err(e) = self.union_point(at, n) {
                out = Err(e);
            }
        }
        out
    }

    /// **A plane-pair line lying on a cylinder's lateral face** — the other twin of
    /// [`Aliases::record`], read off the gate's record ([`crate::planes::SharedRuling`]) rather
    /// than off an operand's corners. Every named class that cuts the line `t1 ∩ t2` — across the
    /// axis or slanted to it — meets it in one point, and every name that point has is folded here
    /// ([`Aliases::record_cylinder_point`]): the three-plane `{c, t1, t2}`, each secant's crossing
    /// with `c` (the root on the recorded line), a tangent partner's double root, and a sectioned
    /// partner's root, the one the intervals leave. Asked of every class cutting the line, not of
    /// the operands' corners: where the prism reaches past the cylinder's cap, the point on the cap
    /// is no corner of either operand, and it is the same point all the same.
    pub(super) fn record_shared_ruling(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        cyl: usize,
        sr: &crate::planes::SharedRuling,
    ) {
        let (t1, t2) = sr.line;
        for c in 0..jd.planes.len() {
            if jd.world_rat(c).is_none() || jd.plane_pair_dir_sign(t1, t2, c) == 0 {
                continue; // only a named class that cuts the line meets it in a point it can state
            }
            let mut s = [c, t1, t2];
            s.sort_unstable();
            let Some(t) = combinatorics::canonical_triple(jd, &s) else {
                continue;
            };
            // The point by its most exact name: a named partner's crossing with `c`, the root the
            // record holds on the line (the record knows the line exactly where a partner with no
            // world statement leaves the three-plane name to the intervals).
            let at = [t1, t2]
                .into_iter()
                .find_map(|x| record_root(jd, cyls, cyl, sr, c, x))
                .unwrap_or(NodeId::three_planes(t));
            // A seed: a point a refused union leaves on another cylinder keeps its names apart.
            let _ = self.record_cylinder_point(jd, cyls, cyl, at, &s);
        }
    }

    /// **A ruling a sectioned class draws and a grazing class touches along** — the gate's record
    /// ([`crate::planes::GrazeLine`]). Every class `k` across the axis whose station falls on a lateral
    /// face meets the line in one point: the three-plane `{k, section, graze}` and the root of
    /// `k ∩ section` the intervals leave ([`combinatorics::root_through`] — the record is the word that
    /// the point lies on the cylinder), folded as every seed's point is
    /// ([`Aliases::record_cylinder_point`]). A station off the faces is left alone: the point is on the
    /// infinite cylinder only, and folding it would hand a corner there (a box's, standing past the
    /// lateral) a pierce name a class with no world statement cannot build into a vertex — unlike
    /// an operand's corner of three named planes on a cylinder's surface, whose names are all exact
    /// ([`seed_from_operands`]).
    ///
    /// The fold puts the grazing class into the name of a point off its faces, which the disk
    /// shortcut in `nesting` would misread («a corner a silent class names lies on its faces»). It does
    /// not, because the representative is the pierce name (`rep_rank` puts an interval pierce before
    /// a three-plane one, and rings carry representatives) — asserted here, where the fold is made.
    pub(super) fn record_graze_line(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        cyls: &[crate::planes::WorkingCyl],
        cyl: usize,
        gl: &crate::planes::GrazeLine,
    ) {
        let Some(wcy) = cyls.get(cyl) else { return };
        let def = &wcy.def;
        let m = def.dir();
        for k in 0..jd.planes.len() {
            let Some(wk) = jd.world_rat(k) else {
                continue;
            };
            if !class_carries_circle(&[wk[0], wk[1], wk[2]], &m) {
                continue; // only a class across the axis meets the line in a point
            }
            let Some(t) = crate::planes::axis_param_of_plane(&wk, def) else {
                continue;
            };
            if !gl.spans.iter().any(|[lo, hi]| *lo <= t && t <= *hi) {
                continue;
            }
            let mut s = [k, gl.section, gl.graze];
            s.sort_unstable();
            let Some(tri) = combinatorics::canonical_triple(jd, &s) else {
                continue;
            };
            let point = NodeId::three_planes(tri);
            // A seed: a point a refused union leaves on another cylinder keeps its names apart.
            let _ = self.record_cylinder_point(jd, cyls, cyl, point, &s);
            debug_assert!(
                !combinatorics::names_on_cylinder(jd, cyls, cyl, point, &s)
                    .iter()
                    .any(|&n| combinatorics::pierce_name(n).is_some())
                    || combinatorics::pierce_name(self.canon_point(point)).is_some(),
                "a point folded onto a grazed ruling is represented by its pierce name"
            );
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

    /// **Whether two names are one point by this table** ([`combinatorics::same_name`] with
    /// [`Aliases::canon_point`]).
    pub(crate) fn one_point(&self, a: NodeId, b: NodeId) -> bool {
        combinatorics::same_name(&|n| self.canon_point(n), a, b)
    }

    /// **The representative of a point's names** — the one every key uses (`canon_point`).
    ///
    /// Among three-plane names it is the least, which is [`combinatorics::canonical_triple`]'s
    /// own answer (the operand's name and the arrangement's discoveries meet there).
    /// **A point on a cylinder is represented on the cylinder.** A `Pierce` name
    /// locates the point exactly *and* says which cylinder it lies on — the lateral chart's
    /// geometry (θ about the axis, `pierce_meet`) reads that from the name — while a three-plane
    /// name of the same point (a class through a tangent corner: `ThreePlane([cap, t, wc])`) is
    /// a key only. So a pierce name outranks a three-plane name, and an exact pierce name one an
    /// interval alone can place: a road that reads the representative on a line of two named
    /// classes reads an exact name as that line's crossing, and an interval name there as a second
    /// point it cannot order against the first. An interval name still outranks a three-plane one:
    /// it says which cylinder the point is on. Among names of one rank the least. Either way the
    /// representative is a function of the class alone.
    fn rep_rank(&self, n: NodeId) -> (u8, NodeId) {
        let rank = match n.kind() {
            NodeKind::Pierce { planes, .. }
                if planes.iter().any(|&p| self.nameless.get(p) == Some(&true)) =>
            {
                1
            }
            NodeKind::Pierce { .. } => 0,
            NodeKind::ThreePlane(_) => 2,
        };
        (rank, n)
    }

    fn union_point(&mut self, a: NodeId, b: NodeId) -> Result<(), TwoCylinders> {
        self.join_points(self.lines.restate(a), self.lines.restate(b))
    }

    /// [`Aliases::union_point`] of the names as given — the table's keys are restated names
    /// ([`Lines::restate`]), and this is what puts a key learnt before its line was in the
    /// component of its restatement.
    ///
    /// **Every union of the table is made here, and here a point on two cylinders is refused.** A
    /// component holding a pierce name is represented by one (`rep_rank` puts every pierce name
    /// before every three-plane one), and every pierce name of a component stands on its
    /// representative's cylinder — this refusal is what keeps it so — so two components stand on
    /// two cylinders exactly when their representatives do. Such a point is where two cylinders
    /// meet, which no road names: the union is not made and does not count as learnt, and the caller
    /// says what that means — a seed keeps the names apart, a road that proved them one point
    /// refuses the boolean.
    fn join_points(&mut self, a: NodeId, b: NodeId) -> Result<(), TwoCylinders> {
        let (ra, rb) = (self.find_point(a), self.find_point(b));
        if ra == rb {
            return Ok(());
        }
        if let (Some((_, ca, _)), Some((_, cb, _))) = (
            combinatorics::pierce_name(ra),
            combinatorics::pierce_name(rb),
        ) && ca != cb
        {
            return Err(TwoCylinders);
        }
        self.joined += 1;
        let (lo, hi) = if self.rep_rank(ra) < self.rep_rank(rb) {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.point.insert(hi, lo);
        self.point.entry(lo).or_insert(lo);
        Ok(())
    }

    /// The identity to key a vertex on — itself when nothing was merged with it.
    pub(crate) fn canon_point(&self, t: NodeId) -> NodeId {
        if self.point.is_empty() && self.lines.is_empty() {
            return t; // the overwhelmingly common case pays nothing
        }
        // Restated before the find, as every key is; and after it, since a key learnt before its
        // line names the component by a pair the line no longer goes by.
        self.lines.restate(self.find_point(self.lines.restate(t)))
    }

    /// The name to call a line by on `class` — `w` itself when no other wall carries it.
    pub(crate) fn canon_wall(&self, class: usize, w: usize) -> usize {
        self.lines.wall(class, w)
    }

    /// The walls that carry **one line** with `w` on `class` — `w` alone when none does.
    ///
    /// The wall fold says two classes name one line; this reads that back out, because a point on
    /// that line lies on **every** plane in the family, and the point fold needs to hear about it
    /// (see the report in [`split_at_crossings`]).
    pub(super) fn wall_family(&self, class: usize, w: usize) -> Vec<usize> {
        self.lines.family(class, w)
    }

    /// Learn everything `other` holds. A union refused here is two tables' knowledge putting one
    /// point on two cylinders, `CylinderPairContact` ([`Aliases::join_points`]).
    pub(super) fn absorb(
        &mut self,
        jd: &Judge<'_, WorkingPlane>,
        other: &Aliases,
    ) -> Result<(), RejectReason> {
        if self.nameless.is_empty() {
            self.nameless.clone_from(&other.nameless);
        }
        for line in other.lines.learnt() {
            self.learn_line(jd, &line.held, line.cut);
        }
        for (&k, &v) in &other.point {
            self.union_point(k, v)?;
        }
        Ok(())
    }

    /// Learn that every class in `held` (ascending, three or more) holds one line `cut` cuts. A
    /// pierce name already in the table under another pair of that line is joined with its
    /// restatement, so the table keeps one component per point whichever pair a key was learnt by.
    fn learn_line(&mut self, jd: &Judge<'_, WorkingPlane>, held: &[usize], cut: usize) {
        let named = |c: usize| jd.world_rat(c).is_some();
        let sense = |p, q, cut| combinatorics::line_sense(jd, p, q, cut);
        self.learn_line_by(held, cut, named, sense);
    }

    /// [`Aliases::learn_line`] with the planes' facts handed in ([`Lines::learn`]).
    fn learn_line_by(
        &mut self,
        held: &[usize],
        cut: usize,
        named: impl Fn(usize) -> bool,
        sense: impl Fn([usize; 2], [usize; 2], usize) -> i8,
    ) {
        if !self.lines.learn(held, cut, named, sense) {
            return;
        }
        self.joined += 1;
        let stale: Vec<NodeId> = self
            .point
            .keys()
            .copied()
            .filter(|&k| self.lines.restate(k) != k)
            .collect();
        for k in stale {
            // A restatement keeps the name's cylinder, so this union is never refused.
            let joined = self.join_points(self.lines.restate(k), k);
            debug_assert!(
                joined.is_ok(),
                "a key and its restatement stand on one cylinder"
            );
        }
    }

    /// **How much has been learned** — the fixed-point loop's measure: the unions that joined two
    /// components, of names and of lines. Each one changes the partition and nothing else does, so the
    /// loop has converged exactly when a round adds none. The number of keys is no such measure: a
    /// union of two components whose roots are both keys already re-points one of them and leaves the
    /// count where it was.
    pub(super) fn learned(&self) -> usize {
        self.joined
    }
}

/// **The crossing of named classes `c` and `x` on the line a [`crate::planes::SharedRuling`] records**
/// — `Pierce{[c, x]}` at the root the record holds on the line (a tangency's `Double` likewise);
/// `None` where `x` has no world statement or no root of the pair is on the line.
fn record_root(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    cyl: usize,
    sr: &crate::planes::SharedRuling,
    c: usize,
    x: usize,
) -> Option<NodeId> {
    use crate::combinatorics::QuadRoot;
    use nacre_exact::quad::{CylinderMeet, QuadVal};
    let def = &cyls.get(cyl)?.def;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let (wc, wx) = (jd.world_rat(c)?, jd.world_rat(x)?);
    let on_line = |line: &nacre_exact::quad::MeetLine, sv: &QuadVal| {
        sr.holds(&wx, &m, line, sv) == Some(true)
    };
    match nacre_exact::quad::plane_plane_cylinder(&wc, &wx, &o, &m, r2)? {
        CylinderMeet::Pair { line, s } => [(QuadRoot::Lo, &s[0]), (QuadRoot::Hi, &s[1])]
            .into_iter()
            .find(|(_, sv)| on_line(&line, sv))
            .map(|(root, _)| NodeId::pierce(c, x, cyl, root)),
        CylinderMeet::Tangent { line, s } if on_line(&line, &QuadVal::from_rat(s)) => {
            Some(NodeId::pierce(c, x, cyl, QuadRoot::Double))
        }
        _ => None,
    }
}

/// **The lines plane classes hold, one name each.** A line held by three or more classes (a pencil)
/// is one entry: the set of every class that holds it, and every pair of that set maps to it — so a
/// class asks the same set whichever of its partners it names the line by, and two pencils learnt
/// apart that share two classes are one line. A line held by two classes alone has no entry; its
/// pair is its name.
///
/// **The line's name is its least pair in the representative's order** ([`Line::key`] — a pair of
/// two classes with a world statement before any other, then by number; `Aliases::rep_rank` reads
/// a pierce name the same way, an interval-placed name being one whose pair holds a class with no
/// world statement). The pierce names on the line are restated onto that pair
/// ([`Lines::restate`]), so a restated name is never greater than the name it restates: a key the
/// table learnt before its line can never be a component's representative, and the representative
/// does not depend on the order the unions and the lines were learnt in. A named pair keeps the
/// names exact wherever the line has one.
#[derive(Default, Debug, Clone)]
struct Lines {
    /// Every pair of classes of a learnt line → that line's index in `lines`, and which way the
    /// pair runs along the line's name ([`combinatorics::line_sense`]).
    of: HashMap<[usize; 2], (usize, i8)>,
    /// The learnt lines. A line merged into another is left empty.
    lines: Vec<Line>,
}

/// One learnt line ([`Lines`]).
#[derive(Default, Debug, Clone)]
struct Line {
    /// Its classes, ascending.
    held: Vec<usize>,
    /// Which of them have a world statement, beside `held`.
    named: Vec<bool>,
    /// A class whose plane cuts it — what [`combinatorics::line_sense`] reads the pairs' senses
    /// against.
    cut: usize,
}

impl Line {
    /// A pair of the line's classes in the representative's order: both with a world statement
    /// first, then the ascending pair by number.
    fn key(&self, a: usize, b: usize) -> (bool, [usize; 2]) {
        let named = |c: usize| self.held.binary_search(&c).is_ok_and(|i| self.named[i]);
        (!(named(a) && named(b)), [a.min(b), a.max(b)])
    }

    /// The pair the line is named by — its least pair ([`Line::key`]).
    fn name(&self) -> [usize; 2] {
        let mut best = (true, [usize::MAX; 2]);
        for i in 0..self.held.len() {
            for j in (i + 1)..self.held.len() {
                best = best.min(self.key(self.held[i], self.held[j]));
            }
        }
        best.1
    }
}

impl Lines {
    fn line(&self, a: usize, b: usize) -> Option<&Line> {
        let pair = if a < b { [a, b] } else { [b, a] };
        self.of.get(&pair).map(|&(i, _)| &self.lines[i])
    }

    fn is_empty(&self) -> bool {
        self.of.is_empty()
    }

    /// The wall `class` calls the line `class ∩ w` by: the other class whose pair with `class` is
    /// least ([`Line::key`]) — on a class of the line's name, the name's other class.
    fn wall(&self, class: usize, w: usize) -> usize {
        self.line(class, w)
            .and_then(|line| {
                line.held
                    .iter()
                    .copied()
                    .filter(|&x| x != class)
                    .min_by_key(|&x| line.key(class, x))
            })
            .unwrap_or(w)
    }

    /// The walls that carry the line `class ∩ w` on `class` — every other class that holds it.
    fn family(&self, class: usize, w: usize) -> Vec<usize> {
        match self.line(class, w) {
            Some(line) => line.held.iter().copied().filter(|&x| x != class).collect(),
            None => vec![w],
        }
    }

    /// Every learnt line.
    fn learnt(&self) -> impl Iterator<Item = &Line> {
        self.lines.iter().filter(|l| !l.held.is_empty())
    }

    /// **A pierce name restated onto its line's name** — the pair a learnt line is called by, the
    /// root carried along it ([`combinatorics::QuadRoot::along`]). Any other name as it is.
    fn restate(&self, n: NodeId) -> NodeId {
        let Some((planes, cyl, root)) = combinatorics::pierce_name(n) else {
            return n;
        };
        let Some(&(i, sense)) = self.of.get(&planes) else {
            return n;
        };
        let [a, b] = self.lines[i].name();
        if [a, b] == planes || sense == 0 {
            return n; // the line's own name — or a sense the predicate did not decide: kept apart
        }
        NodeId::pierce(a, b, cyl, root.along(sense > 0))
    }

    /// Learn that every class in `held` (ascending, three or more) holds one line that `cut` cuts,
    /// joining it with every learnt line that shares a pair with it. `named` says which classes
    /// have a world statement (the planes' order), `sense` which way one pair runs along another
    /// against `cut` ([`combinatorics::line_sense`]). `true` when that changed the table.
    fn learn(
        &mut self,
        held: &[usize],
        cut: usize,
        named: impl Fn(usize) -> bool,
        sense: impl Fn([usize; 2], [usize; 2], usize) -> i8,
    ) -> bool {
        let mut set: Vec<usize> = held.to_vec();
        let mut joined: Vec<usize> = Vec::new();
        loop {
            let before = set.len();
            for i in 0..set.len() {
                for j in (i + 1)..set.len() {
                    if let Some(&(k, _)) = self.of.get(&[set[i], set[j]])
                        && !joined.contains(&k)
                    {
                        joined.push(k);
                        set.extend(self.lines[k].held.iter().copied());
                    }
                }
            }
            set.sort_unstable();
            set.dedup();
            if set.len() == before {
                break;
            }
        }
        if joined.len() == 1 && self.lines[joined[0]].held == set {
            return false;
        }
        let at = joined.iter().copied().min().unwrap_or(self.lines.len());
        if at == self.lines.len() {
            self.lines.push(Line::default());
        }
        for &k in &joined {
            self.lines[k] = Line::default();
        }
        let line = Line {
            named: set.iter().map(|&c| named(c)).collect(),
            held: set,
            cut,
        };
        let name = line.name();
        for i in 0..line.held.len() {
            for j in (i + 1)..line.held.len() {
                let pair = [line.held[i], line.held[j]];
                let along = if pair == name {
                    1
                } else {
                    sense(pair, name, cut)
                };
                self.of.insert(pair, (at, along));
            }
        }
        self.lines[at] = line;
        true
    }
}

#[cfg(test)]
#[path = "../tests/arrangement/aliases.rs"]
mod tests;
