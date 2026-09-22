use super::*;
/// **Which arcs of the circle around an edge a wall fills**, over one sub-interval of that edge.
///
/// An arrangement edge on plane class `W` is the meeting of `W` with one other plane, so the little
/// circle around it is cut into just two arcs — above `W` and below. A contributing face either
/// fills both (it passes through `W` there) or exactly one (its own boundary lies along the edge
/// there, and it hangs off to one side). That is the whole content of this enum, and `edge_mask`
/// turns it into the label flip. It is the BRep spelling of a Nef local pyramid restricted to two
/// planes, and the flip itself is binary winding propagation (Zhou et al. 2016).
///
/// ★ **The occupancy is a fact about the sub-interval, not about the face.** One face can fill both
/// arcs over one stretch of `W ∩ fp` and one arc over the next — a notched or holed face does
/// exactly that where its boundary rides the line. So a segment must be *homogeneous* in this kind:
/// [`split_at_crossings`] subdivides segments later and every piece inherits the kind verbatim, so
/// a mixed segment silently mislabels the pieces that disagree with it. Phase C therefore closes and
/// reopens wherever the kind changes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SegKind {
    /// `L` runs through the face's **interior** here, so the face is on both sides of `W` and
    /// crossing this segment flips the above-label and the below-label **together**. `mat` is which
    /// side of the directed line `W ∩ fp` the solid's material lies on (`+1`/`-1`, relative to
    /// `d = n_W × n_fp`).
    Transversal { mat: i8 },
    /// A boundary edge of a face lying **in** `W`. The face is the boundary between material and
    /// void on one side of `W`, so crossing it flips **one** label — the side the body occupies
    /// (`body_above`, from the face's own outward normal against `W`'s stored normal).
    Seated { body_above: bool },
    /// `L` runs along the face's **boundary** here, so the face fills one arc only and crossing this
    /// segment flips **one** label — `body_above`, the side it fills, from [`run_body_above`].
    ///
    /// Not "the face merely touches `W`": a face may pass clean through `W` elsewhere and still fill
    /// one arc over *this* stretch. Kept distinct
    /// from `Seated` so `edge_mask` can let it override a coincident seated rim at a reflex
    /// dihedral, where the two disagree on which side flips.
    Graze { body_above: bool },
    /// The **tangent ruling** of a face lying in `W`: the face ends where its cylinder
    /// begins, tangent to `W`. Beyond the line what `W` separates is decided by which side the
    /// cylinder bends toward — its **axis** side: a convex fillet's axis is on the body side and
    /// that side turns void past the line (the cusp under the arc); a hole's axis is on the void
    /// side and that side turns material. So crossing this edge flips the axis side's label, one
    /// bit, whichever side the body was on.
    Tangent { axis_above: bool },
}

/// One segment of a solid's trace on a plane class, named entirely in plane-class triples.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Seg {
    /// Canon plane class the segment rides besides the cut plane (its line is `W ∩ wall`).
    pub wall: usize,
    /// Endpoint identities, canonized.
    pub end: [NodeId; 2],
    /// The same two endpoints as **handles on this segment's own line**: the third plane pinning
    /// each on `wc ∩ wall`. Carried, not recovered from `end`, because a canonical name need not
    /// mention either of this line's planes.
    pub end_h: [combinatorics::EndPin; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
}

/// A **full circle** of a solid's trace on a plane class: a cylinder's mark, closed —
/// no endpoints, no wall, no place in the segment machinery. Seated circles come from disk
/// faces and circular holes lying in the class; transversal circles from a lateral surface
/// crossing it, and circles join the arrangement only at the cell stage.
///
/// ★★ **That it is *full* is checked, not inherited.** [`split_circles`] asks it of the segments
/// themselves, so a crossed circle is **split into arcs** rather than treated as the closed cell
/// it is not. Inheriting it from the population gate's wall rule would make a promise about
/// *circles* rest on a rule about *walls*.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CircleTrace {
    /// The cylinder class ([`ClassIx::Cyl`] payload) whose circle this is.
    pub cyl: usize,
    pub solid: SolidSide,
    pub kind: SegKind,
    /// **The angular extent this contribution covers**, counter-clockwise `[from, to]` about the
    /// cylinder's axis — `None` for the whole circle.
    ///
    /// ★★ A lateral face with a **hole** contributes over part of its circle and nothing over the
    /// rest, and a rim of that hole contributes a *different* kind from the band around it. The
    /// straight sibling [`RulingTrace`] has carried ends since it was written, for the same
    /// reason; a circle carried none because before holes there was no partial circle to state.
    /// [`split_circles`] cuts at these ends and hands each arc only the contributions that cover
    /// it.
    pub arc: Option<[NodeId; 2]>,
}

/// One circle of the class's arrangement after merging: both operands' contributions on one
/// cylinder class, plus the cylinder's exact statement (the containment predicates read it).
/// The circle twin of [`MergedSeg`].
#[derive(Clone, Debug)]
pub(crate) struct MergedCircle {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    /// Each contribution with the **angular extent** it covers — see [`CircleTrace::arc`]. `None`
    /// is the whole circle, which is what every contribution was before a lateral face could have
    /// a hole.
    pub merged: Vec<(SolidSide, SegKind, Option<[NodeId; 2]>)>,
}

impl MergedCircle {
    /// The contributions as `edge_mask` reads them.
    ///
    /// ★ A circle only reaches a mask **uncut**, and a partial contribution forces its own cut
    /// ([`split_circles`] feeds every extent end into the split), so every extent here is `None`.
    /// A `Some` would mean a whole-circle mask was about to be built from a partial trace — the
    /// silently wrong answer this vessel exists to prevent — so it is refused by name.
    pub(super) fn whole_marks(&self) -> Result<Vec<(SolidSide, SegKind)>, BoolError> {
        self.merged
            .iter()
            .map(|&(s, k, arc)| match arc {
                None => Ok((s, k)),
                Some(_) => Err(reject(RejectReason::PartialCircleUncut)),
            })
            .collect()
    }
}

/// **One arc of a circle a segment cut** — the DCEL edge a crossed circle becomes.
///
/// ★★ **A circle leaves the parallel road here.** An uncut circle is a one-edge closed cell that
/// the orbit walk cannot express (orbits need ≥ 2 half-edges and a vertex to turn at), which is why
/// [`walk_cells`] appends it *after* the walk on pseudo-half-edges. Cut, it has endpoints and is
/// an ordinary edge of the walk.
///
/// ★ `end` is in the arc's **own travel order**, `end[0]` → `end[1]` counter-clockwise about the
/// circle's normal. That is not a convention this crate invents: STEP's `EDGE_CURVE` runs from
/// `edge_start` to `edge_end` along increasing parameter, and `nacre-step` already writes a
/// circular edge that way — so the exporter needs nothing, and the STEP round-trip measures it.
#[derive(Clone, Debug)]
pub(crate) struct MergedArc {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    pub end: [NodeId; 2],
    /// The contributions that **cover this arc** — not the circle's list entire.
    ///
    /// ★★ Not inherited whole from the circle: a face with a **hole** is not there at all over
    /// the hole, and along the hole's rim it grazes where the band around it crosses.
    /// [`split_circles`] selects per arc; `edge_mask` reads the result unchanged.
    ///
    /// ★ Read by `label_cells`' `mask_of`: crossing an arc flips the same bits crossing its circle
    /// would, which is what "a piece of the same trace" means. (It was carried before that consumer
    /// existed, because the split is the only place that knows which circle an arc came from.)
    pub merged: Vec<(SolidSide, SegKind)>,
}

/// **One ruling piece of the class's arrangement** — the straight sibling of [`MergedArc`], for a
/// class **parallel** to a cylinder's axis: the wall plane meets the
/// lateral surface in up to two axis-parallel lines, and a piece of one is an ordinary two-ended
/// edge whose carrier is the cylinder, not a plane pair.
///
/// ★ `end` is in the piece's **own travel order**: `end[0]` → `end[1]` runs along `+m` (the
/// cylinder's axis direction) — the straight reading of [`MergedArc`]'s CCW convention, read the
/// same way (`edge_at`: the even half-edge travels the stated way, its twin the other).
///
/// `side` names which of the two parallel rulings ([`combinatorics::RulingCarrier::side`] — the
/// sign of `(x − o) · (m × n̂)` against the class's canonical coefficients).
#[derive(Clone, Debug)]
pub(crate) struct MergedRuling {
    pub cyl: usize,
    pub def: nacre_topo::CylinderDef,
    pub side: i8,
    pub end: [NodeId; 2],
    /// The `(solid, kind)` contributions, like a segment's — `edge_mask` reads it unchanged.
    pub merged: Vec<(SolidSide, SegKind)>,
    /// The lateral's `orient_sign` (`+1` material inside the cylinder, `−1` outside — a boss or a
    /// bore), for the vertical answer's content check; an instrument's field, like the label.
    #[cfg(test)]
    pub orient: i8,
}

/// Group a class's circle traces by cylinder class, in ascending class order (deterministic —
/// replay mints handles from this order). The def comes from the class table, which a
/// `ClassIx::Cyl` index indexes directly.
pub(super) fn merge_circles(
    circles: &[CircleTrace],
    cyls: &[crate::planes::WorkingCyl],
    aliases: &Aliases,
) -> Result<Vec<MergedCircle>, BoolError> {
    // ★ The merge layer is where names become canonical — `merge_coincident` for a
    // segment's ends, here for an extent's, `merge_rulings` for a ruling's — so everything after
    // it (splits, the walk, the chart, assembly) sees one name per point.
    let canon = |e: [NodeId; 2]| [aliases.canon_point(e[0]), aliases.canon_point(e[1])];
    let mut out: Vec<MergedCircle> = Vec::new();
    let mut sorted: Vec<&CircleTrace> = circles.iter().collect();
    sorted.sort_by_key(|c| c.cyl);
    for c in sorted {
        if let Some(last) = out.last_mut() {
            if last.cyl == c.cyl {
                last.merged.push((c.solid, c.kind, c.arc.map(canon)));
                continue;
            }
        }
        // ★ The class table answers directly. This used to scan the face rows for a lateral face
        // of the same class and `expect` one — a runtime claim where the truth is structural: a
        // `ClassIx::Cyl` index is an index *into this table*, built from the same map that issued
        // it ("rebuilt from the same map so the two cannot drift", `dense_planes`).
        let def = cyls[c.cyl].def.clone();
        out.push(MergedCircle {
            cyl: c.cyl,
            def,
            merged: vec![(c.solid, c.kind, c.arc.map(canon))],
        });
    }
    Ok(out)
}

/// **One piece of a solid's ruling trace on a ∥ class**: a lateral face's
/// mark on a plane class that contains the cylinder's axis, along one of the two axis-parallel
/// lines. The ruling sibling of [`CircleTrace`], with ends because a ruling is not closed:
/// `end[0]` → `end[1]` ascends the axis, and both are Pierce names (`{wc, ⊥ plane, cyl, root}` —
/// real classes, found by the tracer).
///
/// ★★ **Not "spanning the face's own rims" — a hole makes that false.** Where the
/// face is buried in the other body it does not *cross* the wall its hole's vertical edges lie on
/// — it ends at it — so one ruling comes in pieces of different [`SegKind`]s and each is its own
/// trace. See [`ruling_sweep`]. A face with no hole still yields exactly one piece per ruling,
/// rim to rim.
#[derive(Clone, Debug)]
pub(crate) struct RulingTrace {
    /// The cylinder class, with `side` the ruling's identity
    /// ([`combinatorics::RulingCarrier::side`]).
    pub cyl: usize,
    pub side: i8,
    pub end: [NodeId; 2],
    pub solid: SolidSide,
    pub kind: SegKind,
    /// The lateral's `orient_sign`, carried to [`MergedRuling::orient`].
    #[cfg(test)]
    pub orient: i8,
}

/// Group the class's ruling traces into [`MergedRuling`]s — the ruling twin of
/// [`merge_circles`]: dedupe identical rulings (two operands stating one), collect their
/// contributions, def from the class table. Deterministic order: `(cyl, side descending, extent)` —
/// a fixed rule, ascending cyl with the `+1` ruling first.
///
/// ★★★★★ **The extent is part of the key, and leaving it out is a silent duplicate waiting for a
/// population.** The fold below is *adjacency*-based: it merges a trace into the previous one when
/// the two name the same edge. While each `(cyl, side)` carries exactly one trace per operand, the
/// sort puts those two side by side and the fold works. The moment one ruling comes in **pieces** —
/// a lateral face with a hole grazes part of its own ruling instead of crossing it — a stable sort
/// on `(cyl, side)` alone leaves operand A's pieces before operand B's, so the two statements of one
/// piece are **not adjacent** and each becomes its own `MergedRuling`: one edge appearing twice,
/// with half its contributions each. Sorting by the extent too puts identical pieces together
/// whatever order the tracer emitted them in.
///
/// ☑ Measured over the suite: **no `(cyl, side)` group carries more than one trace today**, so the
/// fold below never actually fires. Both facts are written for
/// the population that has not arrived yet — a cylinder belongs to one operand, so two statements
/// of one ruling need a shape nothing builds so far.
pub(super) fn merge_rulings(
    rulings: &[RulingTrace],
    cyls: &[crate::planes::WorkingCyl],
    aliases: &Aliases,
) -> Vec<MergedRuling> {
    let mut out: Vec<MergedRuling> = Vec::new();
    // Canonical ends (see `merge_circles`) — the key below is by name.
    let mut sorted: Vec<(&RulingTrace, [NodeId; 2])> = rulings
        .iter()
        .map(|r| {
            (
                r,
                [aliases.canon_point(r.end[0]), aliases.canon_point(r.end[1])],
            )
        })
        .collect();
    sorted.sort_by_key(|(r, end)| (r.cyl, -r.side, *end));
    for (r, end) in sorted {
        if let Some(last) = out.last_mut() {
            if last.cyl == r.cyl && last.side == r.side && last.end == end {
                last.merged.push((r.solid, r.kind));
                continue;
            }
        }
        out.push(MergedRuling {
            cyl: r.cyl,
            def: cyls[r.cyl].def.clone(),
            side: r.side,
            end,
            merged: vec![(r.solid, r.kind)],
            #[cfg(test)]
            orient: r.orient,
        });
    }
    out
}

/// A solid's trace on one plane class. `declined` non-empty ⇒ incomplete: a consumer must not
/// read "no segments" as "the plane misses the solid".
#[derive(Default, Debug)]
pub(crate) struct Trace {
    pub segs: Vec<Seg>,
    /// The closed circle elements beside the segments — see [`CircleTrace`].
    pub circles: Vec<CircleTrace>,
    /// The rulings beside them — see [`RulingTrace`].
    pub rulings: Vec<RulingTrace>,
    /// Names this trace found to denote one feature — see [`Aliases`].
    pub aliases: Aliases,
    /// Single-point tangential contacts — real arrangement vertices, but not segments (a
    /// zero-length chord would abort at `Line::through_points`).
    pub touches: Vec<NodeId>,
    /// `(face index, reason)` for every face this brick could not trace.
    pub declined: Vec<(usize, DeclineKind)>,
}
