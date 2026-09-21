use super::*;
/// [`trace_result_faces`] for the band pass's tests, with the reuse mode fixed at `Proved` so
/// the cylinder guard inside is what turns it off.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn trace_result_faces_full_for_test(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    n_a: usize,
    class_owner: &[Option<SolidSide>],
    trace_in: &combinatorics::TraceInput,
) -> Result<(Vec<LocalFace>, Curved, Option<BoolError>), BoolError> {
    trace_result_faces(
        model,
        kind,
        a,
        b,
        jd,
        faces,
        plane_ix,
        cyls,
        n_a,
        class_owner,
        crate::reuse::ClassReuse::Proved,
        trace_in,
    )
}

/// Per `(cylinder class, plane class)`, the four bits of that circle's disk cell — see the band
/// pass. Keyed rather than positional because a class carries a circle only when a cylinder cuts
/// it, and a cylinder cuts only some classes.
pub(crate) type DiskLabels = HashMap<(usize, usize), Label>;

/// **A cut circle's per-arc disk-side labels** — the sector sibling of [`DiskLabels`] (the
/// rulings ladder): a cut circle bounds no whole disk, so "the material immediately above and
/// below this plane inside the circle" is answered **per arc**, by the label of the cell on the
/// arc's disk side. Keyed like `DiskLabels`; each entry lists one [`ArcLabel`] per arc in the
/// split's own arc order. The frame is the disk labels' (the class's stored normal).
pub(crate) type ArcLabels = HashMap<(usize, usize), Vec<ArcLabel>>;

/// One arc's row in [`ArcLabels`] — what the arc bounds, and **who traced it**.
///
/// ★★★★★ **Two different questions live here and they are not one question.** `label` answers
/// *membership* — "is this side material" — and that is all a label can say. Whether this lateral
/// face is even **here** to bound it is a different proposition, and the trace answers that one:
/// `marks` is the contribution list that covered this arc ([`MergedArc::merged`]), where a face
/// running through says `Transversal` and one whose boundary stops at the arc says `Graze`.
/// A **holed** lateral re-entering a boolean makes the difference visible — its two sectors can
/// carry a literally identical label, and only the kind says which one is a face at all.
///
/// ★ The two are set together, in one pass over one arc, for the reason `tri_pt3` and `rotated`
/// are: split into two maps they could disagree, and then nothing could say which was the truth.
#[derive(Clone, Debug)]
pub(crate) struct ArcLabel {
    /// The arc's end pair, in the split's own arc order.
    pub(crate) ends: [NodeId; 2],
    /// The four bits of the cell on the arc's disk side, in the class's stored frame.
    pub(crate) label: Label,
    /// The `(solid, kind)` contributions covering this arc — [`MergedArc::merged`] verbatim.
    ///
    /// ★★★ **A lateral face's mark here is never `Seated`, and the type already says so.** The
    /// circle arm of [`trace_one`] builds its kind from [`CylOnClass`], whose only two answers are
    /// `Crosses` and `Grazes` — there is no `Seated` arm to take, because a cylinder's lateral
    /// surface cannot lie *in* a plane. `Seated` circle contributions come from the seated arm one
    /// branch further down, which a face enters only when its own class **is** this class: a
    /// planar face's rim. So "skip `Seated`" is not a guess about producers — it reads the one
    /// thing on this list that a lateral could not have written.
    pub(crate) marks: Vec<(SolidSide, SegKind)>,
}

/// **One ruling piece, as the cylinder's own chart needs it** — the wall it rides and the
/// **axis interval** it spans.
///
/// ★★★★★ **Carried out rather than recomputed.** `rulings_on_class` already decides this extent
/// (the rulings ladder), and a chart that worked it out again would be a **second source of one
/// fact** — the defect shape this repository names first. So the value leaves by the road
/// [`CutRims`] already takes: made per class, folded once, keyed in the global class space.
///
/// ★ `end` is the piece's own two pierce nodes and `z` their axis parameters
/// ([`node_axis_param`]). Both are carried because the chart needs the **order** (from the node's
/// `(MeetLine, QuadVal)` name, via `circular_order_about_seam`) and the **position** (`z`), and
/// deriving one from the other twice is how a name loses a sign.
/// ☑ Measured over the whole suite: **826 ruling pieces, every one with both ends named** — so a
/// piece whose ends have no ⊥ partner (which would have no `z` at all) is not a population today.
///
/// ★★★★ **There is no `side` here, and that is the point.** The first spelling carried one, and
/// when the cells arrived (D1b) nothing read it: a ruling's identity is `(wall class, root)`, and
/// which of the two parallel rulings that is derives from it through the one production spelling
/// ([`crate::combinatorics::ruling_side`], which is how `cyl_chart::emit_lateral` gets it at
/// emission time). A carried copy of a
/// derived value is the second spelling this ladder keeps being bitten by.
///
/// ★★★★★ **`cfg_attr(not(test), ...)`, not a blanket allow.** The only reader is `cyl_chart`,
/// which is `#[cfg(test)]`, so outside a test build these fields have none — but *inside* one they
/// must genuinely be read, and that is what caught `side`. A plain `#[allow(dead_code)]` would
/// have kept carrying it silently, which is how the field survived a whole rung.
#[derive(Clone, Debug)]
pub(crate) struct RulingExtent {
    pub(crate) wall: usize,
    /// Which of the wall's two rulings — or `0`, a **tangent** wall's single one, the
    /// one station that carries no label: nothing changes across it, the face ends there.
    pub(crate) side: i8,
    pub(crate) end: [NodeId; 2],
    pub(crate) z: [nacre_exact::Rat; 2],
    /// **The chart's vertical answer** (capability D, third rung): the label of the cell this
    /// ruling borders **inside** the cylinder, on the wall class this ruling lies on.
    ///
    /// ★★★★★ **This is the half `disk_labels`/`arc_labels` do not carry.** A [`Label`] holds the
    /// material on **both** sides of its plane, so a ⊥ class's disk label already answers the
    /// chart's *horizontal* crossings entire — which is why `cyl_chart::Chart::read_cell` reads two
    /// of them and asks them to agree. Crossing a *vertical* line is crossing the **wall**, and the
    /// cell inside the strip is where that plane's two sides are stated at the lateral.
    ///
    /// `None` when the wall's rational name and its stored normal cannot be related
    /// ([`world_rat_sense`]) — the ruling is then left without an answer and **counted**, never
    /// guessed at.
    ///
    /// ★ It ships, not an instrument — the emitter does not read the horizontal lines only: a
    /// cut end
    /// the reader cannot pair with its rim leaves a cell with no horizontal answer at all, and
    /// this is what answers it ([`crate::cyl_chart::Chart::read_cell`]).
    pub(crate) label: Option<Label>,
    /// The `(solid, kind)` contributions that covered this piece — [`MergedRuling::merged`], the
    /// same list [`ArcLabel::marks`] carries for an arc.
    ///
    /// ★★★★★ **A label answers *membership*; this answers *existence*.** The same split
    /// holds on the ⊥ side: the chart collects **every** perpendicular class, and a ruling's
    /// extent is
    /// set by its *wall*, not by this lateral face — so a vertical line can be in the chart while
    /// the face is not there at all. `cyl_chart::read_cell`'s `face_spans` already states the
    /// rule for exactly this list, and it is carried beside the label for the reason
    /// `ArcLabel` carries both: split into two maps they could disagree, and then nothing
    /// could say which was the truth.
    #[cfg(test)]
    pub(crate) marks: Vec<(crate::planes::SolidSide, SegKind)>,
}

/// **Which of a ruling's two half-edges borders the cell inside the cylinder** — `true` for the
/// even one (the piece's own `end[0] → end[1]`, ascending the axis).
///
/// ★★★★★ **Derived from the definitions, not re-spelled.** Three sentences already in this crate
/// compose to the answer, and none of them is written a second time here:
///
/// * [`combinatorics::RulingCarrier::side`] is `sign((x − o) · (m̂ × n̂))` against the class's
///   **canonical rational** name, so the strip's interior lies along `−side · (m̂ × n̂_r)`; lifted
///   to the stored normal by [`world_rat_sense`] (`κ`) that is `−side·κ·(m̂ × n̂_P)` — the
///   product [`plus_theta_is_above`] already spells.
/// * The walk keeps a cell **on the left of its travel in the root face's outward frame**,
///   `n_out = frame_sign · n̂_P` ([`crate::planes::WorkingPlane::frame_sign`] — the sentence lives
///   there), and the even half-edge travels `+m̂` ([`MergedRuling::end`]) — so its cell lies
///   along `n_out × m̂ = −frame_sign · (m̂ × n̂_P)`.
///
/// ```text
///   even half-edge's cell is interior  ⟺  −frame_sign·(m̂ × n̂_P) ∥₊ −side·κ·(m̂ × n̂_P)
///                                      ⟺  side · κ · frame_sign = +1
///                                      ⟺  plus_theta_is_above ≠ (frame_sign > 0)
/// ```
///
/// — letter for letter the ⊥ road's rule for a cut circle's disk side (`axis_up ≠ (frame_sign >
/// 0)`, `emit_faces`), which is the other half-edge rule the chart frame enters.
///
/// ★★★★★ **`frame_sign` is a factor, not an ornament.** `side · κ` alone is frame-free —
/// true of *which ruling* it names,
/// false of *which half-edge's cell*: the chart's left is the root face's, not the stored
/// normal's. An unmoved corpus never tells them apart because a wall class is `frame_sign = −1`
/// only when its face lies on a **seed plane** with its outward along +axis (`Model::new` plants
/// x = 0, y = 0, z = 0 with cache direction −axis), which the commuting
/// oracle reaches by putting the plate's max faces there by a translation and turning walls onto
/// them by a rotation — 105 cells, every one caught by `ruling_probe::SIDE_CHECK` at the fact.
pub(super) fn ruling_interior_is_even(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    side: i8,
) -> Option<bool> {
    // The strip's interior is the side `+θ̂` does **not** enter, read in the chart's frame.
    Some(plus_theta_is_above(jd, wc, side)? != (jd.planes[wc].frame_sign > 0))
}

/// What the arrangement learned about the curved boundary, bundled: the band pass reads
/// `disk_labels`, the assembly reads `cut_rims`. One struct so the trace's return does not grow
/// element by element (it was widened once already, for `deferred`).
pub(crate) struct Curved {
    pub(crate) disk_labels: DiskLabels,
    pub(crate) arc_labels: ArcLabels,
    pub(crate) cut_rims: CutRims,
    /// ★ The alias table the class world settled on — every name in the labels, rims
    /// and rulings above is its representative, and the lateral chart, which mints station
    /// names of its own, asks it (`canon_point`) so a station on a corner *is* the corner.
    pub(crate) aliases: Aliases,
    /// Per cylinder class, the ruling pieces on it — the **vertical** lines of that cylinder's
    /// chart, where `disk_labels`/`arc_labels` carry the horizontal ones' answers.
    ///
    /// ★ **Its consumer is `cyl_chart`** — production: `emit_lateral` names
    /// every ruling from these (`Chart::node_on`, `Chart::ruling_name`), and the census reads their
    /// `label`/`marks`.
    pub(crate) rulings: HashMap<usize, Vec<RulingExtent>>,
}

/// **The seam table — every node the result faces reference, realized to a coordinate and a
/// measured tolerance.** The weld table `assemble_fuse_cut` reads; built directly from the
/// emitted rings (no `build_seam`: that is raw-index and pierce-only), rejecting rather than
/// panicking on a degenerate meet.
///
/// ★ A named function rather than a block for the same reason `per_class` is one: the arc fence
/// calls it on the very faces production feeds it. The deferred stopper intercepts the whole
/// stretch this runs in, so a failure *here* never reaches an arc population's caller — which
/// means no reject name can testify that the pierce arm works, and only a direct second consumer
/// can (measured: with the arm disabled wholesale, every boolean-level fence stays green).
pub(crate) fn seam_table(
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    jd: &Judge<'_, WorkingPlane>,
) -> Result<Vec<SeamVertex>, BoolError> {
    watch!(SEAM);
    let geom = jd.planes;
    let mut seam: Vec<SeamVertex> = Vec::new();
    let mut seen: HashMap<NodeId, ()> = HashMap::new();
    for f in faces {
        for loop_ in f.poly_rings() {
            for &node in loop_.iter() {
                if seen.insert(node, ()).is_some() {
                    continue;
                }
                // ★★ **A pierce node is realized from its name, like everything else
                // here: the truth is the definition, the coordinate its cache.** The
                // coordinate is `a + b√c` — `pierce_point` re-solves it from the name's
                // `(line, s)`; a rational road cannot hold it (`node_coords_rat`'s doc
                // calls that a type fact, not a width decline). The tolerance is the same
                // rule as the three-plane arm below: how far the realized point sits from
                // each surface that defines it, plus the closed-form pairwise meet
                // (`pierce_vertex_tol` carries the argument for which pairwise curves are
                // in and out).
                //
                // ★ The **vertex minting** past this table is `boolean`'s
                // `def_triple`/`node_handle`, which names a pierce node's vertex as
                // `Vertex::Pierce` — a cut rim's node and the scan's crossing on a
                // ruling both travel that road.
                if let Some(([p0, p1], cyl, _)) = combinatorics::pierce_name(node) {
                    let wcy = &cyls[cyl];
                    let arr = combinatorics::pierce_point(jd, cyl, &wcy.def, node)
                        .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                    let point = nacre_math::Point3::from_array(arr);
                    seam.push(SeamVertex {
                        point,
                        triple: node,
                        tol: crate::planes::pierce_vertex_tol(
                            point,
                            &geom[p0].plane,
                            &geom[p1].plane,
                            &wcy.realized,
                        ),
                    });
                    continue;
                }
                // Declining, never `continue`: a skipped seam entry surfaces downstream as
                // `MissingSeam`, whose class is `SuspectedDefect` and whose sentence is
                // "a reconstruction dropped a crossing" — a wrong diagnosis for an input
                // the kernel simply does not build yet.
                let t = three_plane_name(node)
                    .ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?;
                let point = three_planes(&geom[t[0]].plane, &geom[t[1]].plane, &geom[t[2]].plane)
                    .ok_or_else(|| reject(RejectReason::ThreePlanes))?;
                seam.push(SeamVertex {
                    point,
                    triple: node,
                    tol: vertex_tol(
                        point,
                        &geom[t[0]].plane,
                        &geom[t[1]].plane,
                        &geom[t[2]].plane,
                    ),
                });
            }
        }
    }
    // **Two names, one point.** Every arrangement vertex is a distinct plane triple, and the
    // materialized coordinate is only its cache — so two *different* triples landing on the same
    // coordinate means the exact substrate and the f64 cache disagree about how many vertices
    // there are. Downstream that becomes a zero-length edge, so catch it here, where both triples
    // are still in hand, instead of letting `assemble_fuse_cut` discover it as a degenerate line.
    //
    // The usual cause is a **split plane table**: one geometric plane carried by two classes, whose
    // triples then name one point twice (measured — two `add_cuboid` walls at the same
    // x that `planes_coplanar` could not prove coplanar because their un-normalized coefficients
    // are not exactly proportional). A genuine 4-plane concurrency does the same.
    for (i, u) in seam.iter().enumerate() {
        for v in &seam[i + 1..] {
            if u.point == v.point {
                return Err(reject(RejectReason::SeamAlias));
            }
        }
    }
    Ok(seam)
}

/// Every result face across all plane classes, before assembly (the driver's risky half, testable
/// by face count without mutating the model). A declining class aborts the whole boolean.
///
/// Classes are visited in order, so the faces are produced in a sequence that is a function of the
/// input — which is what keeps `assemble_fuse_cut`'s handle minting replayable.
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_result_faces(
    model: &Model,
    kind: BoolKind,
    a: Handle<Solid>,
    b: Handle<Solid>,
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    cyls: &[crate::planes::WorkingCyl],
    n_a: usize,
    class_owner: &[Option<SolidSide>],
    reuse: crate::reuse::ClassReuse,
    trace_in: &combinatorics::TraceInput,
) -> Result<(Vec<LocalFace>, Curved, Option<BoolError>), BoolError> {
    let planes = jd.planes;
    let mut local_faces: Vec<LocalFace> = Vec::new();

    // **What each class contributes, for the classes the other operand cannot reach.** Everything
    // else stays `Arrange`, which is the whole engine as it was.
    // ★ **A cylinder in the operands turns the reuse shortcut off entirely**.
    // `reuse::pass_through` moves faces **per plane class** (`plane_ix[fi].plane() != wc →
    // continue`), and a lateral face belongs to no plane class — so a class decided
    // `PassThrough` would carry the planes across and leave the cylinder's own faces behind,
    // silently. The guard also keeps `VertexClasses::of` (which walks face slots and asks each
    // for its plane row) away from cylinder rows, so it does double duty.
    let has_cyl = plane_ix.iter().any(|c| matches!(c, ClassIx::Cyl(_)));
    let reuse = if has_cyl {
        crate::reuse::ClassReuse::Off
    } else {
        reuse
    };
    let plans = crate::reuse::class_plans(model, reuse, kind, a, b, planes, class_owner);

    // ★ Two passes, because an identity must not depend on the order classes happen to be visited.
    // Pass A traces and splits every class, learning aliases as it goes; pass B builds the cells.
    // Doing both in one loop would key an early class's tables before a later class had reported
    // the alias that renames one of its vertices — the same point would then be assembled under two
    // names and rejected as a `SeamAlias`, for a reason that is really "we asked too early".
    //
    // `split_at_crossings`' output does not depend on the table (it only writes to it), so pass A
    // is a pure prefix of the old work rather than extra work — except that a *new* alias found
    // during a split leaves the merge that ran before it stale, so pass A repeats until the table
    // stops growing. Discoveries only accumulate and are bounded, so this terminates; a model with
    // no concurrency at all makes exactly one round.
    // ★ **Every class that has anything to arrange, in a fixed order.** Only the classes an operand
    // face lies on: a result face on `W` is part of `∂A` or `∂B`, so `W` carries an operand face.
    let work: Vec<usize> = {
        // ★ **Plane classes only.** A cylinder row rides in `trace_in` too (its lateral face is
        // what leaves circles on the ⊥ classes), and a lateral surface *is* not a plane class —
        // asking which one it is has no answer, which is exactly what `ClassIx::plane`'s panic
        // says. Filtering here is the upstream filter that panic is a detector for; without it
        // the first cylinder boolean past the C2 stopper aborts the kernel.
        let mut c: Vec<usize> = trace_in
            .faces
            .iter()
            .flatten()
            .filter_map(|(fp, _)| match plane_ix[*fp] {
                ClassIx::Plane(i) => Some(i),
                ClassIx::Cyl(_) => None,
            })
            .collect();
        c.sort_unstable();
        c.dedup();
        c
    };
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, jd, cyls, trace_in);
    #[allow(clippy::type_complexity)]
    let mut splits: Vec<(Vec<MergedSeg>, Vec<MergedCircle>, Vec<MergedRuling>)> = Vec::new();
    loop {
        let before = aliases.len();
        // **Every class in the round sees the table as it stood when the round began**, and
        // its own discoveries on top — where the sequential loop also showed it whatever the
        // lower-numbered classes had found meanwhile. The fixed point is the same, and for
        // two reasons that are both properties of `Aliases` rather than of the schedule:
        // a merged class's representative is its **minimum** element, so the final partition
        // does not depend on the order the unions happened in; and discoveries only
        // accumulate, so a round that learns something later than it used to just costs one
        // more round. The last round — the one whose splits are kept — runs on a table that
        // has stopped growing either way.
        let snapshot = aliases.clone();
        // Class-major, so the reject a decline raises is still the lowest-numbered class's.
        let round = crate::par::try_map_range(work.len(), |k| {
            let wc = work[k];
            // **Pass A runs for every class, including the ones pass B will not arrange.**
            //
            // It used to skip them, and that was unsound: pass B's reuse can *decline* — a vertex
            // where four planes meet has four possible names and only the alias table settles
            // which, so the class falls back to arranging — and the fallback reads `splits[wc]`,
            // which skipping pass A leaves empty. The class's faces would then vanish.
            //
            // Nothing in the code stopped that; it simply needed a model with a concurrency in a
            // region the other operand cannot reach, and the corpus has none. The cheap repair is
            // to keep the fallback a real one, which is what this does.
            let mut tr = timed!(
                TRACE_ON,
                trace_on_class(trace_in, wc, jd, cyls, faces, plane_ix, &snapshot)
            );
            let mut local = snapshot.clone();
            local.absorb(&std::mem::take(&mut tr.aliases));
            // An incomplete trace ⇒ honest reject, naming what the tracer could not do and on
            // which operand face. A class can decline several faces; the first is the one
            // reported, and `try_map_range` picks the lowest-numbered class, which is the one
            // the sequential loop returned at.
            if let Some(&(fp, kind)) = tr.declined.first() {
                // The witness is the face handle, read **kind-agnostically**: a lateral face can
                // decline too (`DeclineKind::CylSpan` is exactly that), and asking it for its
                // plane row would abort where an honest reject belongs.
                return Err(reject(decline_to_reject(kind, faces[fp].face())));
            }
            // ★ A disk cap's chord is one of these segments (both ends pierce-pinned,
            // from the same parity sweep as every polygon's section).
            let merged = timed!(MERGE, merge_coincident(jd, &tr.segs, wc, &local));
            let split = timed!(SPLIT, split_at_crossings(jd, cyls, wc, &merged, &mut local))?;
            let split = drop_newsless(split)?;
            let circles = merge_circles(&tr.circles, cyls, &local)?;
            let rulings = merge_rulings(&tr.rulings, cyls, &local);
            Ok(((split, circles, rulings), local))
        })?;
        splits.clear();
        for (split, local) in round {
            splits.push(split);
            // Absorbing the snapshot back is a no-op; only the round's discoveries are new.
            aliases.absorb(&local);
        }
        if aliases.len() == before {
            break;
        }
    }

    // Pass B is independent per class — it reads `splits[wc]` and the judging context, and
    // returns owned faces — so it is evaluated across cores. `try_map_range` is what keeps
    // that from being observable: the faces are consumed in class order, so the handles
    // `assemble_fuse_cut` mints are the ones a single thread would have minted, and a class
    // that declines surfaces the same rejection the sequential loop returned (the lowest
    // index, not whichever worker got there first).
    // The vertex→class map of each operand, built once and only if some class needs it.
    let vc_a = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::A))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, 0..n_a));
    let vc_b = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::B))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, n_a..faces.len()));

    let per_class = crate::par::try_map_range(splits.len(), |k| {
        let wc = work[k];
        let (split, circles, rulings) = &splits[k];
        // The per-class product: the faces, the disk labels the band pass reads (empty for a
        // class with no circles — and for a reused class, which is why cylinders switch reuse
        // off), the cut rims, and the **stopper socket's** deferred reject (empty
        // — see the socket note below).
        type ClassOut = (
            Vec<LocalFace>,
            Vec<(usize, Label)>,
            Vec<(usize, ArcLabel)>,
            Vec<(usize, CutRim)>,
            Vec<(usize, RulingExtent)>,
            Option<BoolError>,
        );
        let arrange = |wc: usize| -> Result<ClassOut, BoolError> {
            watch!(CELLS);
            // ★★ **One `ClassEdges` for the whole pipeline.** The split renumbers half-edges, so
            // everything below must read the *same* edges the walk did — building it here is what
            // makes that structural instead of a promise. (`frame_audit` runs its own copy of this
            // pipeline and must build it the same way; the arc fence locks that they agree.)
            let edges = timed!(
                C_SPLIT,
                ClassEdges::of(jd, cyls, wc, split, circles, rulings, &aliases)
            )?;
            // ★★★ **A stopper is *made* here — and only made.** A stopper lives in
            // this socket: built where the class-level fact lives, raised inside `reconstruct`
            // at the assembly's very end, intercepting every stage between via the
            // `deferred.unwrap_or(e)` below — so a refused population carries one name out
            // however far the pipeline gets. The arc population is supported and the socket holds
            // `None`; the ladder stays for the next out-of-coverage class (ellipses).
            //
            // ★★ `the_audit_and_the_boolean_agree_about_an_arc_class` compares the audit
            // *with* the boolean, so it stays green when both slide to the same wrong name; asking
            // it about interception measures the proposition next door.
            let staged = per_class(jd, cyls, kind, wc, &edges);
            // ★ **The stopper socket.** A stopper is made here — per class, raised at
            // the assembly's very end, intercepting every stage between (seven layers of
            // `deferred.unwrap_or` down the whole pipeline). The arc population is supported,
            // so the socket holds `None`; the next out-of-coverage class (ellipses,
            // say) plugs its own deferred reject in here and inherits the entire interception
            // ladder instead of re-plumbing it.
            let deferred: Option<BoolError> = None;
            let s = match staged {
                Ok(s) => s,
                // The socket is empty, so clippy sees a literal `None` being unwrapped — the
                // yield's *shape* is the point (a plugged stopper wins here), kept as is.
                #[allow(clippy::unnecessary_literal_unwrap)]
                Err(e) => return Err(deferred.unwrap_or(e)),
            };
            Ok((
                s.faces,
                s.disk_labels,
                s.arc_labels,
                edges.cut_rims.clone(),
                s.ruling_extents,
                deferred,
            ))
        };
        // **The plan decides, and only ever downwards.** A `PassThrough` that cannot name one of
        // its vertices falls back to arranging, so this can lose the shortcut but never the answer.
        let reused = match plans[wc] {
            crate::reuse::ClassPlan::Arrange => None,
            crate::reuse::ClassPlan::Empty => Some(Vec::new()),
            crate::reuse::ClassPlan::PassThrough(side) => {
                watch!(REUSE);
                let (vc, range) = match side {
                    SolidSide::A => (vc_a.as_ref(), 0..n_a),
                    SolidSide::B => (vc_b.as_ref(), n_a..faces.len()),
                };
                vc.and_then(|vc| {
                    crate::reuse::pass_through(
                        model,
                        wc,
                        &planes[wc],
                        faces,
                        plane_ix,
                        range,
                        vc,
                        |t| aliases.canon_point(t),
                    )
                })
            }
        };
        match reused {
            Some(f) => Ok((f, Vec::new(), Vec::new(), Vec::new(), Vec::new(), None)),
            None => arrange(wc),
        }
    })?;
    let mut curved = Curved {
        disk_labels: HashMap::new(),
        arc_labels: HashMap::new(),
        cut_rims: HashMap::new(),
        rulings: HashMap::new(),
        aliases,
    };
    // The first arc class's deferred reject, in `work` order — the map above may run its classes
    // in parallel, but this fold reads the vec in order, so the choice is deterministic.
    let mut deferred: Option<BoolError> = None;
    for (k, (faces, labels, arcs, rims, ruls, d)) in per_class.into_iter().enumerate() {
        local_faces.extend(faces);
        // ★ `work[k]` is the translation from this arrangement's k-th class to the global plane
        // class index — the curved maps are keyed in the global space the assembly speaks.
        for (cyl, label) in labels {
            curved.disk_labels.insert((cyl, work[k]), label);
        }
        for (cyl, al) in arcs {
            curved
                .arc_labels
                .entry((cyl, work[k]))
                .or_default()
                .push(al);
        }
        for (cyl, rim) in rims {
            curved.cut_rims.insert((cyl, work[k]), rim);
        }
        // ★ `wall` was this arrangement's k-th class; restate it in the global space the rest of
        // the map already speaks, exactly as the three keys above do.
        for (cyl, mut r) in ruls {
            r.wall = work[k];
            curved.rulings.entry(cyl).or_default().push(r);
        }
        deferred = deferred.or(d);
    }
    Ok((local_faces, curved, deferred))
}
