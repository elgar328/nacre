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
) -> Result<(Vec<LocalFace>, Curved), BoolError> {
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

/// **A cut circle's per-arc disk-side labels** — the sector sibling of [`DiskLabels`]: a cut circle
/// bounds no whole disk, so "the material immediately above and below this plane inside the circle"
/// is answered **per arc**, by the label of the cell on the arc's disk side. Keyed like
/// `DiskLabels`; each entry lists one [`ArcLabel`] per arc in the split's own arc order. The frame
/// is the disk labels' (the class's stored normal).
pub(crate) type ArcLabels = HashMap<(usize, usize), Vec<ArcLabel>>;

/// One arc's row in [`ArcLabels`] — what the arc bounds, and **who traced it**.
///
/// **Two different questions live here and they are not one question.** `label` answers
/// *membership* — "is this side material" — and that is all a label can say. Whether this lateral
/// face is even **here** to bound it is a different proposition, and the trace answers that one:
/// `marks` is the contribution list that covered this arc ([`Piece::merged`]), where a face
/// running through says `Transversal` and one whose boundary stops at the arc says `Graze`.
/// A **holed** lateral re-entering a boolean makes the difference visible — its two sectors can
/// carry a literally identical label, and only the kind says which one is a face at all.
///
/// The two are set together, in one pass over one arc, for the reason `tri_pt3` and `rotated`
/// are: split into two maps they could disagree, and then nothing could say which was the truth.
#[derive(Clone, Debug)]
pub(crate) struct ArcLabel {
    /// The arc's end pair, in the split's own arc order.
    pub(crate) ends: [NodeId; 2],
    /// The four bits of the cell on the arc's disk side, in the class's stored frame.
    pub(crate) label: Label,
    /// The `(solid, kind)` contributions covering this arc — [`Piece::merged`] verbatim.
    ///
    /// **A lateral face's mark here is never `Seated`, and the type already says so.** The
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
/// **Carried out rather than recomputed.** `rulings_on_class` already decides this extent, and a
/// chart that worked it out again would be a **second source of one fact** — the defect shape this
/// repository names first. So the value leaves by the road [`CutRims`] already takes: made per
/// class, folded once, keyed in the global class space.
///
/// `end` is the piece's own two pierce nodes: the chart reads both its θ (the node's
/// `(MeetLine, QuadVal)` name, via `circular_order_about_seam`) and the lines it stops on (the
/// planes of the names) off them.
#[derive(Clone, Debug)]
pub(crate) struct RulingExtent {
    pub(crate) wall: usize,
    /// Which of the wall's two rulings — or `0`, a **tangent** wall's single one, the
    /// one station that carries no label: nothing changes across it, the face ends there.
    pub(crate) side: i8,
    pub(crate) end: [NodeId; 2],
    /// **The chart's vertical answer**: the label of the cell this
    /// ruling borders **inside** the cylinder, on the wall class this ruling lies on.
    ///
    /// **This is the half `disk_labels`/`arc_labels` do not carry.** A [`Label`] holds the
    /// material on **both** sides of its plane, so a ⊥ class's disk label already answers the
    /// chart's *horizontal* crossings entire — which is why `cyl_chart::Chart::read_cell` reads two
    /// of them and asks them to agree. Crossing a *vertical* line is crossing the **wall**, and the
    /// cell inside the strip is where that plane's two sides are stated at the lateral.
    ///
    /// `None` for a tangent ruling (`side == 0`) only: a class with no world statement reads its
    /// stored facing ([`crate::arrangement::plus_theta_is_above`]), and every half-edge of the
    /// class's walk bounds a cell — the chart's census asserts the rest are formed.
    ///
    /// It ships, not an instrument — the cell reader reads it beside the rims' labels
    /// ([`crate::arrangement::cyl_chart::Chart::read_cell`]): it answers a cell no rim speaks
    /// for, and where both speak the two must agree or the cell has no chamber.
    pub(crate) label: Option<Label>,
    /// The label of the cell this ruling borders **outside** the cylinder, on the same wall — the
    /// other half-edge's. Where secant walls run through one ruling and the lateral keeps the
    /// material outside the cylinder, the chambers that decide whether the lateral runs on across
    /// the line are the ones between neighbouring walls out there, and only these cells state them
    /// ([`crate::arrangement::cyl_chart::Chart::slit_at`]). `None` exactly where `label` is.
    pub(crate) outer: Option<Label>,
    /// The `(solid, kind)` contributions that covered this piece — [`MergedRuling::merged`], the
    /// same list [`ArcLabel::marks`] carries for an arc.
    ///
    /// **A label answers *membership*; this answers *existence*.** The same split
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
/// **Derived from the definitions, not re-spelled.** Three sentences already in this crate
/// compose to the answer, and none of them is written a second time here:
///
/// * [`combinatorics::RulingCarrier::side`] is `sign((x − o) · (m̂ × n̂))` against the class's
///   **canonical rational** name, so the strip's interior lies along `−side · (m̂ × n̂_r)`; lifted
///   to the plane's facing by the world name's sense (`κ`) that is `−side·κ·(m̂ × n̂_P)` — the
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
/// **`frame_sign` is a factor, not an ornament.** `side · κ` alone is frame-free —
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

/// What the arrangement learned about the curved boundary, bundled: the lateral chart reads the
/// labels against `split_rims` and cuts its rims there, and the assembly's rims are derived from
/// the cleaned faces ([`crate::draft::held_rims`]). One struct so the trace's return does not grow element by
/// element.
pub(crate) struct Curved {
    pub(crate) disk_labels: DiskLabels,
    pub(crate) arc_labels: ArcLabels,
    /// How the arrangement split each circle it cut ([`CutRims`]) — what the arc labels are
    /// indexed by and the lateral's rims are cut at. Not the result's rim: the cleaning pass can
    /// dissolve its nodes, from the cap and the lateral together.
    pub(crate) split_rims: CutRims,
    /// Per cylinder class, the ruling pieces on it — the **vertical** lines of that cylinder's
    /// chart, where `disk_labels`/`arc_labels` carry the horizontal ones' answers.
    ///
    /// **Its consumer is `cyl_chart`** — production: `emit_lateral` names
    /// every ruling from these (`Chart::ruling_name`), and the census reads their `label`/`marks`.
    pub(crate) rulings: HashMap<usize, Vec<RulingExtent>>,
}

/// **The seam table — every node the result faces reference, realized to the cache its vertex
/// receives.** The weld table `assemble_fuse_cut` reads; built directly from the
/// emitted rings, rejecting rather than panicking on a degenerate meet.
///
/// **A node's coordinate is the cache its vertex will receive.** Each node is realized from
/// the definition the minting pushes ([`Def::vertex`] through [`crate::realize::point_cache`] —
/// the push funnel's own question), so the arrangement's point and the model's cache are one
/// answer. The classes' `f64` plane caches are not read for it: two different planes whose caches
/// round to the same bits are parallel in `f64`, and solving on them refused a corner the planes
/// really have. They supply only the construction figure where the realization has no road, as
/// every producer's figure does.
///
/// A named function rather than a block for the same reason `per_class` is one: the arc fence
/// calls it on the very faces production feeds it — a direct second consumer, which is what
/// testifies that the pierce arm works (measured: with the arm disabled wholesale, every
/// boolean-level fence stays green).
pub(crate) fn seam_table(
    model: &Model,
    faces: &[LocalFace],
    cyls: &[crate::planes::WorkingCyl],
    jd: &Judge<'_, WorkingPlane>,
    planes: &mut crate::realize::PlaneMemo,
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
                // Declining, never `continue`: a skipped seam entry surfaces downstream as
                // `MissingSeam`, whose class is `SuspectedDefect` and whose sentence is
                // "a reconstruction dropped a crossing" — a wrong diagnosis for an input
                // the kernel simply does not build yet.
                let def =
                    Def::of_name(node).ok_or_else(|| reject(RejectReason::PierceVertexUnnamed))?;
                // The construction figure, asked only where the realization has no road. A
                // pierce node's is re-solved from its name — the coordinate is `a + b√c` and
                // `pierce_point` reads it off the name's `(line, s)`; a rational road cannot hold
                // it (`node_coords_rat`'s doc calls that a type fact, not a width decline).
                let figure = || {
                    match def {
                        Def::Pierce { cyl, .. } => {
                            combinatorics::pierce_point(jd, cyl, &cyls[cyl].def, node)
                                .map(nacre_math::Point3::from_array)
                        }
                        Def::Three(t) => {
                            let t = t.planes();
                            three_planes(&geom[t[0]].plane, &geom[t[1]].plane, &geom[t[2]].plane)
                        }
                    }
                    .ok_or_else(|| reject(RejectReason::ThreePlanes))
                };
                let cache = crate::realize::point_cache(
                    model,
                    &def.vertex(jd, geom, cyls)?,
                    planes,
                    figure,
                )?;
                seam.push(SeamVertex {
                    cache,
                    triple: node,
                });
            }
        }
    }
    // **Two names, one coordinate.** Every arrangement vertex is a distinct name, and where the
    // realization answers its coordinate is the nearest `f64` of the truth — so two names on one
    // coordinate are either **one point named twice** (the defect this is here for) or two points
    // closer than an `f64` tells apart, which no output can hold as two vertices either.
    // Downstream both become a zero-length edge, so catch it here, where both names are still in
    // hand, instead of letting `assemble_fuse_cut` discover it as a degenerate line. (The two are
    // not told apart yet — todo 「`SeamAlias` 는 두 원인을 한 딱지로 거절한다」.)
    //
    // One cause of a point named twice is a **split plane table**: one geometric plane carried by two classes, whose
    // triples then name one point twice. Two named planes cannot split — the class merge compares
    // their names — but planes without a world name (a turn off the quarters) are merged by a
    // judgement, and a judgement that fails to prove them one splits them. A genuine 4-plane
    // concurrency does the same.
    for (i, u) in seam.iter().enumerate() {
        for v in &seam[i + 1..] {
            if u.cache.coord() == v.cache.coord() {
                return Err(reject(RejectReason::SeamAlias));
            }
        }
    }
    Ok(seam)
}

/// **The whole result's faces, cleaned once** — the plane arrangement's faces merged per coplanar
/// group ([`crate::assembly::unify_coplanar_faces`]), every cylinder class's lateral faces
/// appended ([`cyl_chart::emit_lateral`]), then one straight-angle pass over the whole list
/// ([`crate::assembly::dissolve_straight_angles`]), and the result's rims read off what that pass
/// left ([`crate::draft::held_rims`]).
///
/// **The pass runs after the laterals exist** because its rule — a node goes only where every ring
/// through it runs straight on — is a sentence about every ring. Run over the plane faces alone it
/// drops a node a lateral's edge still holds: where a flat ends on a lateral along a ruling, and
/// another class's crossings split that ruling, the flat runs straight through the split's nodes
/// and the lateral's ruling edge keeps them — an edge used once. So the lateral cuts its rims and
/// rulings at the split's nodes, and the pass drops a node from the cap or the flat and the lateral
/// together.
///
/// One function rather than a stretch of `boolean`, because tests that assemble a result by hand
/// run this same sequence, and a copy of the order would silently keep the old one.
#[allow(clippy::too_many_arguments)]
pub(crate) fn whole_result_faces(
    kind: BoolKind,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: Vec<LocalFace>,
    curved: &Curved,
    faces_tab: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
) -> Result<(Vec<LocalFace>, HeldRims), BoolError> {
    // Merge coplanar, same-normal faces that share a full edge (e.g. the split side walls a fused
    // coincident interface leaves) so the result is a minimal, chainable solid — a second boolean
    // on it then sees no redundant coplanar planes.
    let mut faces = timed!(
        UNIFY,
        crate::assembly::unify_coplanar_faces(faces, jd, cyls)
    )?;
    if !cyls.is_empty() {
        let rows = crate::bands::cyl_rows(faces_tab, plane_ix, n_a)?;
        // **The lateral faces come from the chart**: every cylinder class's cells, read off the
        // plane arrangement's own labels and emitted as `LocalFace`s, the type the plane side
        // emits.
        let lateral = cyl_chart::emit_lateral(kind, jd, cyls, curved, &rows);
        // The census sits *before* the cleaning pass on purpose: the cleaning dissolves nodes,
        // which would blur the face-by-face question being asked. The census holds the chart
        // against its own rules (and the emitter against the census's independent count), in
        // test builds; the emitter's refusal is handed to it before `?` decides, so a refused
        // class is still recorded.
        #[cfg(test)]
        cyl_chart::census(jd, cyls, kind, &faces, curved, &rows, &lateral);
        // Merged, they would not be: each lateral face is a region of its chart already, so there
        // is no phantom seam to erase (measured: a merging pass over them merges nothing over the
        // suite).
        faces.extend(lateral?);
    }
    timed!(DISSOLVE, {
        let all: Vec<usize> = (0..faces.len()).collect();
        crate::assembly::dissolve_straight_angles(
            &mut faces,
            &all,
            &combinatorics::SharedLines::of(jd, cyls),
        );
    });
    let rims = crate::draft::held_rims(&faces, &curved.split_rims);
    Ok((faces, rims))
}

/// One class's segments, circles and ruling pieces, before the circles and rulings are split.
pub(super) type Skeleton = (Vec<MergedSeg>, Vec<MergedCircle>, Vec<MergedRuling>);

/// **One class's traced edges, split** — pass A's step for one class: the skeleton
/// ([`class_skeleton`]), then its circles and rulings split against its segments
/// ([`split_class_edges`]). Every step learns into `local`, so every name the class's edges carry
/// is born inside the round.
pub(super) fn class_split(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    tr: &Trace,
    local: &mut Aliases,
) -> Result<ClassSplit, BoolError> {
    let (segs, circles, rulings) = class_skeleton(jd, cyls, wc, tr, local)?;
    timed!(
        C_SPLIT,
        split_class_edges(jd, cyls, faces, wc, &segs, &circles, &rulings, local)
    )
}

/// **One class's traced edges, merged and split** — pass A's skeleton step: the segments merged
/// where coincident ([`merge_coincident`]) and cut at their crossings ([`split_at_crossings`],
/// which learns aliases into `local`), which of them carry news asked ([`carries_news`] — its
/// conflicts first, nearest their cause), the circles and rulings merged, the ruling pieces cut at
/// the ends of **every** segment along them ([`cut_rulings_along_segments`] — a newsless knife
/// edge's ends are vertices of the line it lies on), and then the newsless segments dropped. One
/// function because the audit and tests that run a class by hand take the same step; a copy of
/// the order could drift from it.
pub(super) fn class_skeleton(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    tr: &Trace,
    local: &mut Aliases,
) -> Result<Skeleton, BoolError> {
    // A disk cap's chord is one of these segments (both ends pierce-pinned,
    // from the same parity sweep as every polygon's section).
    let merged = timed!(MERGE, merge_coincident(jd, &tr.segs, wc, local))?;
    let split = timed!(SPLIT, split_at_crossings(jd, cyls, wc, &merged, local))?;
    let news = carries_news(&split)?;
    let circles = merge_circles(&tr.circles, local)?;
    let rulings = merge_rulings(&tr.rulings, local);
    let rulings = cut_rulings_along_segments(jd, cyls, wc, &split, rulings, local)?;
    let split = split
        .into_iter()
        .zip(news)
        .filter_map(|(s, keep)| keep.then_some(s))
        .collect();
    Ok((split, circles, rulings))
}

/// **Every class that has anything to arrange, in a fixed order.** Only the classes an operand
/// face lies on: a result face on `W` is part of `∂A` or `∂B`, so `W` carries an operand face.
pub(super) fn work_classes(
    trace_in: &combinatorics::TraceInput,
    plane_ix: &[ClassIx],
) -> Vec<usize> {
    // **Plane classes only.** A cylinder row rides in `trace_in` too (its lateral face is
    // what leaves circles on the ⊥ classes), and a lateral surface *is* not a plane class —
    // asking which one it is has no answer, which is exactly what `ClassIx::plane`'s panic
    // says. Filtering here is the upstream filter that panic is a detector for; without it
    // the first cylinder boolean aborts the kernel.
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
}

/// **The alias rounds: pass A, run until the table stops growing** — the settled table, and each
/// class of `work` traced and split on it ([`class_split`]), or the refusal it met there.
///
/// Two passes, because an identity must not depend on the order classes happen to be visited.
/// Pass A traces and splits every class, learning aliases as it goes; pass B builds the cells.
/// Doing both in one loop would key an early class's tables before a later class had reported
/// the alias that renames one of its vertices — the same point would then be assembled under two
/// names and rejected as a `SeamAlias`, for a reason that is really "we asked too early".
///
/// A *new* alias found during a class's steps leaves what ran before it stale, so pass A repeats
/// until the table stops growing. Discoveries only accumulate and are bounded, so this terminates;
/// a model with no concurrency at all makes exactly one round. The last round — the one whose
/// skeletons are kept — ran on a table that did not grow, so each skeleton is the one the settled
/// table gives.
pub(super) fn settle_rounds(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    trace_in: &combinatorics::TraceInput,
    work: &[usize],
) -> Result<(Aliases, Vec<Result<ClassSplit, BoolError>>), BoolError> {
    let mut aliases = Aliases::default();
    seed_from_operands(&mut aliases, jd, cyls, trace_in);
    let absorb = |table: &mut Aliases, local: &Aliases| table.absorb(jd, local).map_err(reject);
    fixed_point(aliases, absorb, |snapshot| {
        // Every class, refusals included, in class order — the lowest-numbered refusal is the one
        // a settled round reports, in both builds.
        crate::par::map_every(work.len(), |k| {
            let wc = work[k];
            // **Pass A runs for every class, including the ones pass B will not arrange.**
            //
            // Skipping them is unsound: pass B's reuse can *decline* — a vertex where four planes
            // meet has four possible names and only the alias table settles which, so the class
            // falls back to arranging — and the fallback reads the class's skeleton, which skipping
            // pass A leaves empty, so the class's faces would vanish. It takes a model with a
            // concurrency in a region the other operand cannot reach, and the corpus has none;
            // running pass A everywhere keeps the fallback a real one.
            let mut tr = timed!(
                TRACE_ON,
                trace_on_class(trace_in, wc, jd, cyls, faces, plane_ix, snapshot)
            );
            // An incomplete trace ⇒ honest reject, naming what the tracer could not do and on
            // which operand face. A class can decline several faces; the first is the one
            // reported.
            if let Some(&(fp, kind)) = tr.declined.first() {
                // The witness is the face handle, read **kind-agnostically**: a lateral face can
                // decline too (`DeclineKind::CylSpan` is exactly that), and asking it for its
                // plane row would abort where an honest reject belongs.
                return Err(reject(decline_to_reject(kind, faces[fp].face())));
            }
            let mut local = snapshot.clone();
            local
                .absorb(jd, &std::mem::take(&mut tr.aliases))
                .map_err(reject)?;
            let skeleton = class_split(jd, cyls, faces, wc, &tr, &mut local)?;
            Ok((skeleton, local))
        })
    })
}

/// **The fixed point of the alias rounds**, whatever a round computes: `round` maps the table a
/// round begins with to each class's product and the table it learnt on top, or its refusal; the
/// products are absorbed and the rounds repeat until the table stops growing, and that round's
/// products and refusals are the answer.
///
/// **Every class in a round sees the table as it stood when the round began**, and its own
/// discoveries on top — not whatever lower-numbered classes found meanwhile, as a sequential loop
/// would show it. The fixed point is the same, and for two reasons that are both properties of
/// `Aliases` rather than of the schedule: a merged class's representative is its **least** name in
/// one order fixed by the names alone (`Aliases::rep_rank`) — and a name restated onto its line's
/// name is never greater than the name it restates — so neither the partition nor its
/// representatives depend on the order the unions and the lines were learnt in; and discoveries
/// only accumulate, so learning something a round later just costs one more round.
///
/// **A round's refusal is read only on the settled table.** A class that refuses in a round whose
/// table grew is asked again in the next one: what refused it may be a name another class reports
/// in the same round. A refusal carries no table, so only what succeeded is absorbed — a record can
/// fail half-made ([`Aliases::record_same_point`] joins the names before the concurrency that
/// refuses a point on two cylinders), and absorbing the half would grow the table, close the record
/// behind `one_point` next round, and lose the refusal. A class loses only what its own refused
/// steps had learnt, which it learns again. `absorb`'s own refusal ends the rounds at once: two
/// cylinders' components joined is a fact about the geometry that no round undoes.
fn fixed_point<S>(
    mut aliases: Aliases,
    absorb: impl Fn(&mut Aliases, &Aliases) -> Result<(), BoolError>,
    round: impl Fn(&Aliases) -> Vec<Result<(S, Aliases), BoolError>>,
) -> Result<(Aliases, Vec<Result<S, BoolError>>), BoolError> {
    loop {
        let before = aliases.learned();
        let snapshot = aliases.clone();
        let mut classes = Vec::new();
        for r in round(&snapshot) {
            classes.push(match r {
                Ok((product, local)) => {
                    // Absorbing the snapshot back is a no-op; only the round's discoveries are new.
                    absorb(&mut aliases, &local)?;
                    Ok(product)
                }
                Err(e) => Err(e),
            });
        }
        if aliases.learned() == before {
            return Ok((aliases, classes));
        }
    }
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
) -> Result<(Vec<LocalFace>, Curved), BoolError> {
    let planes = jd.planes;
    let mut local_faces: Vec<LocalFace> = Vec::new();

    // **What each class contributes, for the classes the other operand cannot reach.** Everything
    // else stays `Arrange`, which is the whole engine.
    // **A cylinder in the operands turns the reuse shortcut off entirely**.
    // `reuse::pass_through` moves faces **per plane class** (`plane_ix[fi].plane() != wc →
    // continue`), and a lateral face belongs to no plane class — so a class decided
    // `PassThrough` would carry the planes across and leave the cylinder's own faces behind,
    // silently. The guard also keeps `VertexClasses::of` (which walks face slots and asks each
    // for its plane row) away from cylinder rows, so it does double duty.
    // Two more stand on it. `reuse::class_plans` reads a plane that misses every vertex of a solid
    // as missing the solid, which holds only while every face is planar; and a reused class makes
    // no cells and no labels, while a band reads its membership off the disk labels an arranged
    // class makes (pass B below). Taking reuse to inputs with cylinders is therefore a redesign of
    // how a band reads its membership, not a performance step.
    let has_cyl = plane_ix.iter().any(|c| matches!(c, ClassIx::Cyl(_)));
    let reuse = if has_cyl {
        crate::reuse::ClassReuse::Off
    } else {
        reuse
    };
    let plans = crate::reuse::class_plans(model, reuse, kind, a, b, planes, class_owner);

    let work = work_classes(trace_in, plane_ix);
    let (aliases, classes) = settle_rounds(jd, cyls, faces, plane_ix, trace_in, &work)?;
    // The lowest class's refusal, read on the settled table.
    let splits: Vec<ClassSplit> = classes.into_iter().collect::<Result<_, _>>()?;

    // Pass B is independent per class — it reads `splits[k]` and the judging context, and
    // returns owned faces — so it is evaluated across cores. `try_map_range` is what keeps
    // that from being observable: the faces are consumed in class order, so the handles
    // `assemble_fuse_cut` mints are the ones a single thread would have minted, and the refusal is
    // the lowest-index class's in both builds, not whichever worker got there first.
    // The vertex→class map of each operand, built once and only if some class needs it.
    let vc_a = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::A))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, 0..n_a));
    let vc_b = plans
        .contains(&crate::reuse::ClassPlan::PassThrough(SolidSide::B))
        .then(|| crate::reuse::VertexClasses::of(model, faces, plane_ix, n_a..faces.len()));

    let one_class = |k: usize| {
        let wc = work[k];
        // The per-class product: the faces, the disk labels the band pass reads (empty for a
        // class with no circles — and for a reused class, which is why cylinders switch reuse
        // off), the cut rims and the rulings' extents.
        type ClassOut = (
            Vec<LocalFace>,
            Vec<(usize, Label)>,
            Vec<(usize, ArcLabel)>,
            Vec<(usize, CutRim)>,
            Vec<(usize, RulingExtent)>,
        );
        let arrange = |wc: usize| -> Result<ClassOut, BoolError> {
            watch!(CELLS);
            // **One `ClassEdges` for the whole pipeline.** The split renumbers half-edges, so
            // everything below must read the *same* edges the walk did — building it here is what
            // makes that structural instead of a promise. (`frame_audit` reads the same split; the
            // arc fence locks that they agree.)
            let edges = ClassEdges::of(cyls, wc, &splits[k]);
            let s = per_class(jd, cyls, kind, wc, &edges)?;
            Ok((
                s.faces,
                s.disk_labels,
                s.arc_labels,
                edges.cut_rims.to_vec(),
                s.ruling_extents,
            ))
        };
        // **The plan decides, and only ever downwards.** A `PassThrough` that cannot name one of
        // its vertices falls back to arranging, so this can lose the shortcut but never the answer.
        // Decided here, inside the closure, never by shrinking the range: a class that is not
        // arranged keeps its index, so `try_map_range` still reports the lowest-index reject and
        // the faces still come back in class order — the order `assemble_fuse_cut` mints handles in.
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
                    crate::reuse::pass_through(model, jd, wc, faces, plane_ix, range, vc, |t| {
                        aliases.canon_point(t)
                    })
                })
            }
        };
        match reused {
            Some(f) => Ok((f, Vec::new(), Vec::new(), Vec::new(), Vec::new())),
            None => arrange(wc),
        }
    };
    // The lowest-index failure, read in `work` order, so the parallel and serial builds choose
    // alike.
    let per_class = crate::par::try_map_range(splits.len(), one_class)?;
    let mut curved = Curved {
        disk_labels: HashMap::new(),
        arc_labels: HashMap::new(),
        split_rims: HashMap::new(),
        rulings: HashMap::new(),
    };
    for (k, (faces, labels, arcs, rims, ruls)) in per_class.into_iter().enumerate() {
        local_faces.extend(faces);
        // `work[k]` is the translation from this arrangement's k-th class to the global plane
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
            curved.split_rims.insert((cyl, work[k]), rim);
        }
        // `wall` was this arrangement's k-th class; restate it in the global space the rest of
        // the map already speaks, exactly as the three keys above do.
        for (cyl, mut r) in ruls {
            r.wall = work[k];
            curved.rulings.entry(cyl).or_default().push(r);
        }
    }
    Ok((local_faces, curved))
}

#[cfg(test)]
#[path = "../tests/arrangement/result.rs"]
mod tests;
