//! **Which parts of a cylinder's lateral surface survive a boolean** (M6-2a C4b).
//!
//! ★★ **Since the D2b cutover (2026-08-30) production does not read this file's decision roads.**
//! The lateral faces are emitted by `cyl_chart::emit_lateral` from the chart's cells; `band_faces`,
//! `bands_of`, `chamber` and `panel_faces` are `#[cfg(test)]` — the **reference road** the census
//! compares that emission against, face by face, until D3 deletes them. What production still
//! reads here: `cyl_rows` (a face's span is the existence truth at an uncut end), `face_spans`,
//! `read_bits`, `keep_for` and `axis_param` — one spelling each, called from the chart.
//!
//! The plane arrangement decides one plane class at a time; a cylinder's wall is not a plane, so
//! it is decided here instead — and it is decided *coarsely*, because in this population it can
//! be. The M6-2a gate admits only ⊥ cuts and ∥ walls whose **faces** stand clear of the lateral
//! surface, and those two facts together give the **uniform-slab theorem**:
//!
//! > Between two consecutive ⊥ cuts, the other operand's boundary does not meet the open
//! > cylinder slab at all — the ⊥ faces are the cuts themselves and every ∥ wall **face** misses
//! > the lateral. So membership in the other operand is **uniform** over that slab, and one
//! > witness decides the whole band.
//!
//! ★ **Faces, because that is what a boundary is made of.** The gate long tested each ∥ wall's
//! infinite *plane*, which is a cheaper sufficient condition and reads almost the same — but it
//! refused bodies standing well clear of the cylinder whose plane, extended, crossed it. The
//! theorem never asked for that; it asks about the boundary, and the gate now does too.
//!
//! ★★ **A wall face is asked about a rectangle, and the "consecutive cuts" above are per face.**
//! What a ∥ plane meets of the solid cylinder is a rectangle of that plane: the strip across, a
//! lateral face's axis-parameter span along. So the gate reads two separating axes, and a wall
//! clear on either one misses the rectangle (`planes.rs`' `face_clears_footprint`). With several
//! lateral faces there are several rectangles, and that is exactly the granularity the theorem
//! needs: `bands_of` clips every cut to its own row's span, so no band is ever built in the gap
//! between two of them, and a wall sitting in such a gap breaks no premise — there is none there
//! to break. Reading the spans as one `min..max` would invent both the band and the refusal.
//!
//! ★ The theorem replaced a wrong one. "The band is in the other solid iff its z-range is" reads
//! plausibly and is **false**: an L-notch's inner corner is clear of every wall by more than `r`
//! while sitting outside the material. Asking the arrangement is what fixes it, and the
//! counterexample is a fixture of this pass (`a_cylinder_in_the_notch_keeps_no_band`, with its
//! end-to-end twin: cutting with a drill standing in the notch removes nothing).
//!
//! ★★ **The chamber is read off the arrangement, not measured.** The plane arrangement already
//! makes a **disk cell** for every circle a cylinder leaves on a class, and `label_cells` writes
//! four bits on it: which solid's material lies immediately above and below that plane *inside the
//! circle*. That is the band's chamber, stated by the engine that decided it.
//!
//! ★ This replaced a witness ray (C4a's `point_in_faces_rat`): a rational point on the axis, cast
//! through the counterpart's faces, counting crossings. It gave the same answers, but it answered
//! a **3D containment** question — the shape the *component* probe asks — when the band's question
//! is the same shape as "does this face survive": two chambers either side of a boundary. Reading
//! the label needs no coordinates, no ray direction, no abstention retry, and has no width
//! ceiling, and it is why a cylinder may now stand on either side of the boolean.

#[cfg(test)]
use crate::boolean::{Bound, LocalFace};
use crate::combinatorics;
use crate::planes::SolidSide;
use crate::planes::{ClassIx, FaceRow, WorkingPlane, axis_param_of_plane};
use crate::tolerant::Judge;
use crate::{BoolError, BoolKind, RejectReason, reject};
use nacre_scalar::Rat;

/// **One lateral face**, as the band pass reads it.
///
/// ★★ **The row is a face, not a class.** A cylinder *class* is one lateral `Surface` handle
/// (`planes.rs`' `cyl_ix.entry(cf.surf)`), and one surface can carry **several faces of one
/// solid** — bore a plate through and then cut away the middle of the bore, and its wall becomes
/// two disjoint bands. "Which parts survive" is a question about a face; the class only says which
/// surface that face lies on.
///
/// This used to be one row per class, with the first face speaking for the rest. The others' spans
/// were then clipped away in [`bands_of`] and their bands never emitted, which left the result
/// shell open — the boolean came back `OpenResultShell`, naming a symptom of our own omission.
///
/// ★ **Merging the spans into one `min..max` is not the fix either**: it would invent a band across
/// the gap where the solid has no face at all.
pub(crate) struct CylRow {
    /// The [`ClassIx::Cyl`] payload — which lateral surface this face lies on. Several rows may
    /// share it.
    pub(crate) class: usize,
    /// Read by the reference road (test builds); the chart reads `cyls[k].def`, the same value.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) def: nacre_topo::CylinderDef,
    /// **This face's** own extent in the axis parameter (its two rims).
    pub(crate) span: [Rat; 2],
    /// Which operand this face belongs to.
    pub(crate) side: SolidSide,
}

/// One row per lateral **face**, ordered by `(class, lower t)`.
///
/// ★ That order is [`band_faces`]' replay contract stated at the row level. A class's faces have
/// disjoint spans, so sorting by `(class, span[0])` keeps each class's bands in ascending `t` —
/// the contract generalizes rather than bends. With one face per class it *is* the old class-index
/// order, which is why existing results do not move.
///
/// ★★ **That disjointness holds because every lateral face is a full 2π band today.** M6-2b's
/// θ-partial faces will put two faces at the *same* `t`, and the tie then falls to the stable
/// sort's face order — still deterministic, but the sentence above stops being the reason. Naming
/// the premise here so the day it expires is a thing a reader can check, not a surprise.
///
/// A face whose rim span cannot be stated declines by name (`DeclineKind::CylSpan`) — the tracer
/// declines it for the same reason.
pub(crate) fn cyl_rows(
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
) -> Result<Vec<CylRow>, BoolError> {
    let mut out: Vec<CylRow> = Vec::new();
    let mut n_class = 0usize;
    for (i, row) in faces.iter().enumerate() {
        let (ClassIx::Cyl(k), FaceRow::Cylinder(cf)) = (plane_ix[i], row) else {
            continue;
        };
        n_class = n_class.max(k + 1);
        let Some(span) = cf.span else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        // A lateral with no world statement declines the same way a missing span does — the
        // band pass compares against world planes throughout.
        let Some(def) = cf.def.clone() else {
            return Err(reject(RejectReason::TraceDeclined {
                kind: crate::DeclineKind::CylSpan,
                face: cf.face,
            }));
        };
        out.push(CylRow {
            class: k,
            def,
            span,
            side: if i < n_a { SolidSide::A } else { SolidSide::B },
        });
    }
    out.sort_by(|a, b| (a.class, a.span[0]).cmp(&(b.class, b.span[0])));
    // ★ The guard the old per-class table carried, kept as its own sentence: a class exists
    // because a face made it, so a class with no row is a wiring failure rather than an input.
    // It has never fired; an unfired *guard* stays (an unfired *name* would not).
    for k in 0..n_class {
        if !out.iter().any(|r| r.class == k) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
    }
    Ok(out)
}

/// **The surviving lateral bands, as result faces.**
///
/// `plane_faces` is the plane arrangement's output: its circle bounds are where the cuts are, and
/// taking the boundaries from *there* (rather than from every ⊥ class) is what makes each rim
/// shared by exactly the cap face and the band that meet on it.
///
/// Emitted in `(cylinder class, lower t)` order — `reconstruct` mints handles in face order, so
/// this ordering is the replay contract.
/// `labels` is the arrangement's own answer per `(cylinder class, plane class)` — see the module
/// docs. Cylinders may sit on **both** operands: nothing here asks a solid to describe itself.
#[cfg(test)]
pub(crate) fn band_faces(
    kind: BoolKind,
    plane_faces: &[LocalFace],
    rows: &[CylRow],
    jd: &Judge<'_, WorkingPlane>,
    labels: &crate::arrangement::DiskLabels,
    cut_rims: &crate::arrangement::CutRims,
    arc_labels: &crate::arrangement::ArcLabels,
) -> Result<Vec<LocalFace>, BoolError> {
    let mut out = Vec::new();
    for row in rows {
        for (lo, hi) in bands_of(row, plane_faces, jd, cut_rims)? {
            // ★ A **both-cut** interval has a disk cell at neither end — its chambers are per
            // θ-sector, read from the cut circles' arc labels: the panel road (rulings ladder).
            if let (Some(lo_rim), Some(hi_rim)) = (
                cut_rims.get(&(row.class, lo)),
                cut_rims.get(&(row.class, hi)),
            ) {
                panel_faces(
                    kind,
                    row,
                    (lo, lo_rim),
                    (hi, hi_rim),
                    jd,
                    arc_labels,
                    &mut out,
                )?;
                continue;
            }
            let (in_own_inside, in_other) = chamber(jd, row, lo, hi, labels)?;
            // ★ **The wall is a boundary face of its own solid, so that solid's membership flips
            // across it** — read inside from the label, and outside is its negation. The
            // counterpart does *not* flip: the gate keeps its boundary faces clear of the lateral, so
            // the slab theorem covers both sides. The face survives exactly when the two chambers
            // disagree under `keep`.
            let keep_side = |in_own: bool| keep_for(kind, row.side, in_own, in_other);
            let (keep_in, keep_out) = (keep_side(in_own_inside), keep_side(!in_own_inside));
            if keep_in == keep_out {
                continue;
            }
            // The stored cylinder normal points **away from the axis**. It is the result face's
            // outward normal when the kept material is inside the wall (a boss), and must be
            // flipped when the material is outside it (a hole). `keep_in` is that question
            // directly, whichever solid the wall came from.
            out.push(LocalFace {
                surf: ClassIx::Cyl(row.class),
                outer: Bound::Band { lo, hi },
                inner: Vec::new(),
                flip: !keep_in,
            });
        }
    }
    Ok(out)
}

/// The plane classes bounding **this face's** bands, as consecutive `(lo, hi)` pairs ordered by
/// axis parameter.
///
/// The boundaries are the classes that **emitted a circle** for this face's cylinder, plus the
/// classes this face's own rims sit on (found by matching the span's parameters — two ⊥ planes
/// with the same axis parameter *are* one plane, so the match is exact, not a tolerance).
///
/// ★ Circles are gathered per *class*, so a sibling face's boundaries are collected here too — and
/// then clipped away by this face's span. That clip is what keeps two faces of one surface from
/// borrowing each other's cuts.
#[cfg(test)]
pub(crate) fn bands_of(
    row: &CylRow,
    plane_faces: &[LocalFace],
    jd: &Judge<'_, WorkingPlane>,
    cut_rims: &crate::arrangement::CutRims,
) -> Result<Vec<(usize, usize)>, BoolError> {
    let k = row.class;
    let mut classes: Vec<usize> = Vec::new();
    let push = |c: usize, classes: &mut Vec<usize>| {
        if !classes.contains(&c) {
            classes.push(c);
        }
    };
    for lf in plane_faces {
        let ClassIx::Plane(c) = lf.surf else { continue };
        let carries = std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, Bound::Circle { cyl } if *cyl == k));
        if carries {
            push(c, &mut classes);
        }
    }
    // ★ A **cut** circle is a band boundary too (the rulings ladder): it emits no `Circle`
    // bound — its boundary is arc pieces — so the walk above cannot see it, and before this arm
    // the band pass emitted one full-height band straight across it (measured on the lifted
    // through-boss). The arrangement's own record of every cut circle is `cut_rims`, keyed by
    // `(cylinder class, plane class)`; a rim-cut circle (straddle/hung) names the same class the
    // span-end rule below finds, and `push`'s dedup folds them. The final interval list is
    // sorted by axis parameter, so the map's iteration order decides nothing.
    for (&(kk, c), _) in cut_rims.iter() {
        if kk == k {
            push(c, &mut classes);
        }
    }
    for c in 0..jd.planes.len() {
        // ★ A class ⊥ to the axis is a **potential band boundary**, so failing to place it is a
        // refusal, not a skip: silently missing one would merge two bands whose membership
        // differs and hand back a closed, wrong solid. Classes that are not ⊥ cannot bound a
        // band at all — a ∥ wall meets the lateral along **rulings**, which cut bands in θ, not
        // in the axis (the panel road's job; while a pair sits behind the gate's `crossings`
        // record the wall rule has proven it clear of every band) — so those are skipped for
        // cause.
        let Some(coeffs) = combinatorics::class_coeffs_rat(jd, c) else {
            continue;
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        if !nacre_scalar::parallel_rat(&n, &row.def.dir()) {
            continue;
        }
        let t = param(jd, c, row)?;
        if t == row.span[0] || t == row.span[1] {
            push(c, &mut classes);
        }
    }
    let mut keyed: Vec<(Rat, usize)> = classes
        .into_iter()
        .map(|c| Ok((param(jd, c, row)?, c)))
        .collect::<Result<_, BoolError>>()?;
    keyed.retain(|(t, _)| *t >= row.span[0] && *t <= row.span[1]);
    keyed.sort_by_key(|(t, _)| *t);
    Ok(keyed.windows(2).map(|w| (w[0].1, w[1].1)).collect())
}

/// The axis parameter of a plane class, or the honest refusal earned by a class with no exact
/// description — and by one whose parameter is a **value** too wide for `Rat` (the road's own
/// name: the gate's questions are signs and were made total, this one is not).
#[cfg(test)]
fn param(jd: &Judge<'_, WorkingPlane>, c: usize, row: &CylRow) -> Result<Rat, BoolError> {
    axis_param(jd, c, &row.def)
}

/// The same question asked of a **cylinder** rather than one of its faces — the chart road
/// (`crate::cyl_chart`) has a class and a def but no row, and this rule may not be spelled twice.
pub(crate) fn axis_param(
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    def: &nacre_topo::CylinderDef,
) -> Result<Rat, BoolError> {
    param_opt(jd, c, def).ok_or_else(|| reject(RejectReason::WitnessNotRational))
}

fn param_opt(jd: &Judge<'_, WorkingPlane>, c: usize, def: &nacre_topo::CylinderDef) -> Option<Rat> {
    let coeffs = combinatorics::class_coeffs_rat(jd, c)?;
    axis_param_of_plane(&coeffs, def)
}

/// **Which chamber the band sits in, read off the arrangement.**
///
/// The disk cell on a boundary class carries `[inA_above, inA_below, inB_above, inB_below]` — the
/// material immediately above and below that plane *inside the circle*. The band leaves `lo`
/// going up in the axis parameter, so which of "above"/"below" faces it is decided by the class
/// normal against the axis (`sign(n · m)`); at `hi` the band is on the other side.
///
/// ★ **Both ends are read, and they must agree.** The uniform-slab theorem says the counterpart's
/// membership is constant across the open slab; two ends disagreeing means the theorem's premise
/// broke (a ∥ wall crossed the slab — the gate's clearance proof failing), and that is refused
/// rather than resolved by picking one. Free, because both labels are already in hand.
///
/// ★★ **Both solids' bits are data, not a control — and assuming one of them was a real bug.**
/// The first spelling asserted "the band lies inside its own cylinder, so its own solid's bit is
/// true". That holds for a *tool* (a drill solid fills its own cylinder) and is **false for a wall
/// inherited from an earlier bore**: inside that circle the plate has no material — the bore is a
/// hole. Both are ordinary inputs the moment a cylinder may sit on either operand, so the own-solid
/// membership is read from the label too, and the wall's other side follows from what the wall
/// *is*: a boundary face of that solid, so its own membership flips across it while the
/// counterpart's (by the slab theorem) does not.
/// **The θ-panels of a both-cut interval** — the rulings ladder's sector road. Both bounding
/// circles are cut, so the interval's lateral is divided by the rulings at the cut nodes'
/// angles into sectors, and each sector is its own chamber question, answered by the cut
/// circles' per-arc labels ([`crate::arrangement::ArcLabels`]).
///
/// The two rims' nodes must correspond one-to-one by `(wall class, root)` — the same wall
/// cutting both circles puts a branch node at the same angle on each. A geometry whose upper
/// and lower circles are cut by *different* walls has no such pairing and is refused by the
/// ladder's name (its cell has not arrived).
///
/// A kept sector is emitted as a **ring** of the ladder's own vocabulary — `[arc, ruling, arc,
/// ruling]` (`Wall::Arc`/`Wall::Ruling`) — not a `Band`: a band's boundary is two whole
/// circles, a panel's is four pieces.
///
/// ★ **The ring's absolute winding is measured** (cell 4): flipping the emitted orientation
/// (`flip`) turns the through-boss volume oracles red — the watcher this note used to say was
/// still to come.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn panel_faces(
    kind: BoolKind,
    row: &CylRow,
    (lo, lo_rim): (usize, &crate::arrangement::CutRim),
    (hi, hi_rim): (usize, &crate::arrangement::CutRim),
    jd: &Judge<'_, WorkingPlane>,
    arc_labels: &crate::arrangement::ArcLabels,
    out: &mut Vec<LocalFace>,
) -> Result<(), BoolError> {
    use crate::boolean::{Ring, Wall};
    let ladder = || reject(RejectReason::RulingBoundNotYet);
    let k = row.class;
    // A node's `(wall class, root)` — the identity that pairs the two rims' nodes. The wall is
    // the Branch name's plane that is not the rim's own class.
    let name_on = |rim_class: usize,
                   n: combinatorics::NodeId|
     -> Result<(usize, nacre_topo::QuadRoot), BoolError> {
        let (pair, _, root) = combinatorics::branch_name(n).ok_or_else(ladder)?;
        let wall = if pair[0] == rim_class {
            pair[1]
        } else if pair[1] == rim_class {
            pair[0]
        } else {
            return Err(ladder());
        };
        Ok((wall, root))
    };
    let m = lo_rim.nodes.len();
    if m != hi_rim.nodes.len() || m < 2 {
        return Err(ladder());
    }
    let mut hi_of: std::collections::HashMap<(usize, nacre_topo::QuadRoot), combinatorics::NodeId> =
        std::collections::HashMap::new();
    for &n in &hi_rim.nodes {
        if hi_of.insert(name_on(hi, n)?, n).is_some() {
            return Err(ladder()); // two nodes with one name: not this ladder's geometry
        }
    }
    let arc_at = |c: usize, a: combinatorics::NodeId, b: combinatorics::NodeId| {
        arc_labels
            .get(&(k, c))
            .and_then(|v| v.iter().find(|r| r.ends == [a, b]))
            .ok_or_else(ladder)
    };
    let toward_hi = |c: usize| -> bool { crate::planes::plus_t_is_above(&jd.planes[c], &row.def) };
    for j in 0..m {
        let (a_lo, b_lo) = (lo_rim.nodes[j], lo_rim.nodes[(j + 1) % m]);
        let a_hi = *hi_of.get(&name_on(lo, a_lo)?).ok_or_else(ladder)?;
        let b_hi = *hi_of.get(&name_on(lo, b_lo)?).ok_or_else(ladder)?;
        // The sector's chamber, from both rims — the same both-ends-agree discipline as
        // `chamber`'s, on the per-arc labels. The hi arc runs `[a_hi, b_hi]` too: the pairing
        // preserves angle, so the CCW adjacency is the same on both circles.
        let r_lo = arc_at(lo, a_lo, b_lo)?;
        let r_hi = arc_at(hi, a_hi, b_hi)?;
        // ★★★★★ **Existence, before membership.** A label says where *material* is; it does not
        // say whether this lateral face is here to bound it, and a holed lateral coming back in
        // as an operand makes the two sectors' labels literally identical. So the trace is asked
        // first — and where the face is absent there is no membership question to ask.
        let spans_lo = face_spans(r_lo, row.side, toward_hi(lo))?;
        let spans_hi = face_spans(r_hi, row.side, !toward_hi(hi))?;
        #[cfg(test)]
        for r in [r_lo, r_hi] {
            panel_probe::marks(r, row.side);
        }
        if spans_lo != spans_hi {
            // The rims of a hole are themselves band boundaries (`bands_of` pushes every cut
            // rim), so within one interval a sector either exists over its whole height or not
            // at all. Two ends disagreeing means that is false here, and picking an end to
            // believe is the guess `chamber` refuses to make one screenful down.
            return Err(reject(RejectReason::CylinderFaceUndecided));
        }
        if !spans_lo {
            #[cfg(test)]
            panel_probe::sector(false, false);
            continue;
        }
        let (own_lo, other_lo) = read_bits(&r_lo.label, row.side, toward_hi(lo));
        let (own_hi, other_hi) = read_bits(&r_hi.label, row.side, !toward_hi(hi));
        if (own_lo, other_lo) != (own_hi, other_hi) {
            return Err(ladder());
        }
        let keep_side = |in_own: bool| keep_for(kind, row.side, in_own, other_lo);
        let (keep_in, keep_out) = (keep_side(own_lo), keep_side(!own_lo));
        if keep_in == keep_out {
            #[cfg(test)]
            panel_probe::sector(true, false);
            continue;
        }
        #[cfg(test)]
        panel_probe::sector(true, true);
        // The ruling identity at each cut angle, from the one spelling
        // (`arrangement::ruling_side` against the wall class's canonical coefficients).
        let side_at = |n: combinatorics::NodeId| -> Result<i8, BoolError> {
            let (wall, _) = name_on(lo, n)?;
            let w = combinatorics::class_coeffs_rat(jd, wall).ok_or_else(ladder)?;
            let (line, s) = combinatorics::branch_meet(jd, k, &row.def, n).ok_or_else(ladder)?;
            crate::arrangement::ruling_side(&w, &row.def, (&line, &s)).ok_or_else(ladder)
        };
        out.push(LocalFace {
            surf: ClassIx::Cyl(k),
            outer: Bound::Ring(Ring::new(
                vec![a_lo, b_lo, b_hi, a_hi],
                vec![
                    Wall::Arc { cyl: k, ccw: true },
                    Wall::Ruling {
                        cyl: k,
                        side: side_at(b_lo)?,
                        up: true,
                    },
                    Wall::Arc { cyl: k, ccw: false },
                    Wall::Ruling {
                        cyl: k,
                        side: side_at(a_lo)?,
                        up: false,
                    },
                ],
            )),
            inner: Vec::new(),
            flip: !keep_in,
        });
    }
    Ok(())
}

/// **Does this row's lateral face reach into the interval, at this arc?** — the *existence*
/// question, which no label answers.
///
/// A label states where material is. It is written about a **cell**, and it stays the same whether
/// the cylinder's own face bounds that cell or some other face does. That is enough while a
/// lateral marks a class over its whole circle; it stops being enough the moment a **holed**
/// lateral comes back in as an operand, because then two sectors of one rim can carry a literally
/// identical label and differ only in whether the face is there at all.
///
/// The trace already said which. One rule, no shapes counted:
///
/// | this face's mark on the arc | reaches into the interval |
/// |---|---|
/// | [`SegKind::Transversal`] | **yes** — the face passes through the plane here |
/// | [`SegKind::Graze`] | only when `body_above` names the interval's side |
/// | nothing | **no** — the face stops short of this arc |
///
/// `band_is_above` is the interval's side of this rim's plane in that plane's **stored** frame —
/// the very argument [`read_bits`] takes, and produced by the same `toward_hi` the caller already
/// holds. There is no second derivation of "which way is the band" here, deliberately: this
/// ladder has been bitten four times by a sign re-derived one call away from its twin.
///
/// ★ [`SegKind::Seated`] is skipped because it is a *planar* face's word — see [`ArcLabel::marks`],
/// where the type is what rules a lateral out, not a convention.
///
/// ★★ **The rule is not about holes.** An outer rim grazes too (`circle_on_class` says
/// `Grazes { body_above: up }` at the face's own end), and there the interval on the face's side
/// gets `true` from this same test — one sentence about every rim. What is out of reach is only
/// an **uncut** outer rim: it never becomes an [`ArcLabels`] entry at all (a whole circle goes to
/// `DiskLabels` and the band road), so today only cut rims ask here.
///
/// ☑ **How often each row actually fires** (whole binary, production calls only): `Transversal`
/// **265** · `Graze` **13**, of which **1** reaches and **12** do not · nothing at all **0**.
/// The twelve are the six dropped sectors read at both rims, which is the cross-check that the
/// sector census and this one describe the same events.
pub(crate) fn face_spans(
    r: &crate::arrangement::ArcLabel,
    side: SolidSide,
    band_is_above: bool,
) -> Result<bool, BoolError> {
    use crate::arrangement::SegKind;
    let mut answer: Option<bool> = None;
    for (_, kind) in r.marks.iter().filter(|(s, _)| *s == side) {
        let reaches = match *kind {
            SegKind::Seated { .. } => continue, // a planar face's word — see `ArcLabel::marks`
            SegKind::Transversal { .. } => true,
            SegKind::Graze { body_above } => body_above == band_is_above,
        };
        // ★★★ **Two of this solid's faces meeting at one arc and *disagreeing*: refused, and the
        // refusal is a placeholder.** It takes one cylinder class carrying several faces that
        // share a rim — reachable geometry (a split bore whose two bands touch), just not
        // reachable today. The answer that day is almost certainly `.any()`: if any of the
        // solid's faces reaches into the interval, a face is there. It is not written that way
        // now because it cannot be measured, and this ladder has twice shipped an unmeasured rule
        // that turned out wrong. Marks that **agree** decide nothing by themselves, so they are
        // taken: refusing there would be refusing a case with no guess in it.
        if answer.replace(reaches).is_some_and(|prev| prev != reaches) {
            return Err(reject(RejectReason::CylinderFaceUndecided));
        }
    }
    Ok(answer.unwrap_or(false))
}

/// **What the panel road decided, at the layer that decides it** — [`panel_faces`]' instrument.
///
/// ★ The result oracles (volume, `validate`) are the stronger evidence now that the chained
/// operation builds, and they are what pin the *count*. This holds the mechanism instead: that a
/// sector was dropped because the **face was not there**, which no volume can distinguish from a
/// sector dropped by a label. It also survives the day some other refusal moves in front of the
/// assembly again, which is how `ruling_probe` came to exist.
///
/// ★★ **Only [`panel_faces`] records here, never [`face_spans`] itself.** The unit test feeds that
/// function hand-written marks — rows chosen to state the rule, not to describe any geometry — and
/// recording them would make the census describe the test suite instead of the kernel. (That
/// distinction is not academic: this census, once it read production, **refuted** what the cell
/// had written about which rows were unexercised.)
#[cfg(test)]
pub(crate) mod panel_probe {
    use super::SolidSide;
    use crate::arrangement::{ArcLabel, SegKind};
    use std::sync::Mutex;

    /// One entry per sector considered: `(the face reaches this interval, the sector was emitted)`.
    pub(crate) static SECTORS: Mutex<Vec<(bool, bool)>> = Mutex::new(Vec::new());

    /// One entry per rim the panel road read: how many marks of this row's own solid the arc
    /// carried, as `(lateral, seated)`.
    ///
    /// ☑ Measured over the whole binary: 278 reads, `lateral` **always 1**, `seated` **1 once**.
    /// So only the `lateral` half is a zero-population claim and only it is locked; the `Seated`
    /// skip fires for real and is load-bearing. The cell wrote the opposite first, from the four
    /// fixtures it was opening — see the lock.
    pub(crate) static MARKS: Mutex<Vec<(usize, usize)>> = Mutex::new(Vec::new());

    pub(crate) fn sector(spans: bool, kept: bool) {
        SECTORS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .push((spans, kept));
    }

    pub(crate) fn marks(r: &ArcLabel, side: SolidSide) {
        let mine = || r.marks.iter().filter(|(s, _)| *s == side);
        let seated = mine()
            .filter(|(_, k)| matches!(k, SegKind::Seated { .. }))
            .count();
        MARKS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .push((mine().count() - seated, seated));
    }
}

/// The one spelling of "which two bits of a label are this band's chamber": the row's own
/// solid's bit and the counterpart's, on the side of the plane the band occupies. Shared by the
/// whole-disk road ([`chamber`]) and the per-sector panel road, so the two cannot drift.
pub(crate) fn read_bits(
    l: &crate::arrangement::Label,
    side: SolidSide,
    band_is_above: bool,
) -> (bool, bool) {
    let (cyl_bit, other_bit) = match side {
        SolidSide::A => (0usize, 2usize), // [A above, A below, B above, B below]
        SolidSide::B => (2usize, 0usize),
    };
    let i = usize::from(!band_is_above);
    (l[cyl_bit + i], l[other_bit + i])
}

/// The one spelling of the band's keep decision — the wall is a boundary face of its own solid,
/// so that solid's membership flips across it while the counterpart's does not.
pub(crate) fn keep_for(kind: BoolKind, side: SolidSide, in_own: bool, in_other: bool) -> bool {
    match side {
        SolidSide::A => crate::arrangement::keep(kind, in_own, in_other),
        SolidSide::B => crate::arrangement::keep(kind, in_other, in_own),
    }
}

/// **The whole-disk chamber** — both bounding circles uncut, so one question answers the interval.
///
/// ★★★ **This road reads only the label, and that is safe here for a reason worth writing down.**
/// The sector road one screenful up had to ask the *trace* whether the lateral face is even
/// present ([`face_spans`]), because a hole in the face makes two sectors with identical labels.
/// No hole reaches this road: a hole's rims are inside the face (a rim at the face's own end is
/// not an inner loop), they are cut circles, and `bands_of` makes **every** cut rim a band
/// boundary — so an interval a hole passes through has cut ends and goes to the panels. Break
/// that and this function starts answering a question it cannot see.
///
/// ★ So the two roads treat existence differently on purpose, and only one of them has to. They
/// become one question when the lateral gets its own cell complex on its own chart, which is the
/// body of this capability and not this rung.
#[cfg(test)]
pub(crate) fn chamber(
    jd: &Judge<'_, WorkingPlane>,
    row: &CylRow,
    lo: usize,
    hi: usize,
    labels: &crate::arrangement::DiskLabels,
) -> Result<(bool, bool), BoolError> {
    let k = row.class;
    // `above` in the label is the class's **stored** normal side; the band leaves `lo` toward
    // `hi`, i.e. toward increasing axis parameter — which `plus_t_is_above` answers (and where the
    // f64 dot's exactness argument lives).
    let toward_hi = |c: usize| -> bool { crate::planes::plus_t_is_above(&jd.planes[c], &row.def) };
    let read = |c: usize, band_is_above: bool| -> Option<(bool, bool)> {
        let l = labels.get(&(k, c))?;
        Some(read_bits(l, row.side, band_is_above))
    };
    let ends = [(lo, toward_hi(lo)), (hi, !toward_hi(hi))];
    let mut answer: Option<(bool, bool)> = None;
    for (c, band_is_above) in ends {
        let Some(pair) = read(c, band_is_above) else {
            continue; // this end carries no disk cell — the other end speaks
        };
        match answer {
            None => answer = Some(pair),
            // ★ The slab is uniform by the gate's clearance proof — **both** bits must agree at
            // the two ends. Ends disagreeing means that proof failed (something crossed the open
            // slab), and guessing which end to believe is exactly what this kernel does not do.
            // This is also where a flipped side selection surfaces: the two ends read opposite
            // sides, so getting the sign wrong makes them contradict.
            Some(prev) if prev != pair => {
                return Err(reject(RejectReason::CylinderGateUndecided));
            }
            Some(_) => {}
        }
    }
    // Every boundary class of a band carries a disk cell (the emitted circles put them there, and
    // the cylinder's own caps are always arranged) — so this is a wiring failure, not an input.
    answer.ok_or_else(|| reject(RejectReason::CylinderGateUndecided))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planes::{PlaneSetup, plane_index_setup};
    use nacre_math::{Point3, Vector3};
    use nacre_store::Handle;
    use nacre_topo::{Model, Solid};

    /// **The existence rule's whole truth table** — [`face_spans`] against marks written by hand.
    ///
    /// The production lock reaches this through an entire boolean, and a boolean exercises two of
    /// the rows below: a `Transversal` that reaches and a `Graze` pointing away. The rule is one
    /// sentence about every rim, not a description of that fixture, so the rest are stated here —
    /// including the two that can only be reached from geometry the road does not build yet.
    #[test]
    fn face_spans_reads_the_trace_not_the_label() {
        use crate::arrangement::{ArcLabel, SegKind};
        let n = combinatorics::NodeId::ThreePlane([0, 1, 2]);
        let row = |marks: Vec<(SolidSide, SegKind)>| ArcLabel {
            ends: [n, n],
            // Deliberately the label that says "material everywhere": nothing below may read it.
            label: [true; 4],
            marks,
        };
        let (a, b) = (SolidSide::A, SolidSide::B);
        let spans = |marks, above| face_spans(&row(marks), a, above);
        // A face running through the plane is on both sides of it.
        for above in [true, false] {
            assert!(spans(vec![(a, SegKind::Transversal { mat: 1 })], above).unwrap());
        }
        // A face whose boundary stops at the arc is on exactly the side it occupies.
        for body_above in [true, false] {
            for above in [true, false] {
                assert_eq!(
                    spans(vec![(a, SegKind::Graze { body_above })], above).unwrap(),
                    body_above == above,
                    "graze body_above={body_above} against band above={above}"
                );
            }
        }
        // The counterpart's marks are not this row's face, whatever they say.
        assert!(!spans(vec![(b, SegKind::Transversal { mat: 1 })], true).unwrap());
        // A planar face's seated rim says nothing about the lateral — see `ArcLabel::marks`.
        assert!(!spans(vec![(a, SegKind::Seated { body_above: true })], true).unwrap());
        // No mark at all: the face stops short of this arc.
        assert!(!spans(Vec::new(), true).unwrap());
        // Two of this solid's faces on one arc: agreeing decides, disagreeing refuses by name.
        assert!(
            spans(
                vec![
                    (a, SegKind::Graze { body_above: true }),
                    (a, SegKind::Transversal { mat: 1 }),
                ],
                true,
            )
            .unwrap()
        );
        assert!(matches!(
            spans(
                vec![
                    (a, SegKind::Graze { body_above: false }),
                    (a, SegKind::Transversal { mat: 1 }),
                ],
                true,
            ),
            Err(BoolError::Rejected {
                reason: RejectReason::CylinderFaceUndecided,
                ..
            })
        ));
    }

    /// Run the plane arrangement past the C2 stopper and hand the band pass what it needs.
    /// Returns `(band faces, the plane classes' axis parameters keyed by class)`.
    fn bands(
        m: &Model,
        a: Handle<Solid>,
        b: Handle<Solid>,
        kind: BoolKind,
    ) -> (Vec<LocalFace>, Vec<(usize, f64)>) {
        let setup = plane_index_setup(m, a, b).unwrap();
        let PlaneSetup {
            planes: faces_tab,
            geom,
            surf_ix,
            inc_a,
            inc_b,
            plane_ix,
            class_owner,
            n_a,
            standard,
            notes,
            cyls,
            ..
        } = &setup;
        let jd = Judge::new(geom, *standard, notes);
        let trace_in = crate::combinatorics::trace_input(
            m,
            [(a, inc_a), (b, inc_b)],
            surf_ix,
            faces_tab.len(),
            &jd,
            plane_ix,
            cyls,
            Default::default(),
        );
        let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
            m,
            kind,
            a,
            b,
            &jd,
            faces_tab,
            plane_ix,
            cyls,
            *n_a,
            class_owner,
            &trace_in,
        )
        .expect("the drill population traces");
        let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("cylinder rows");
        let out = band_faces(
            kind,
            &plane_faces,
            &rows,
            &jd,
            &curved.disk_labels,
            &curved.cut_rims,
            &curved.arc_labels,
        )
        .expect("bands");
        // The classes' **z**, not their axis parameter: `t` is measured from the cylinder's own
        // origin along its raw `dir`, so a drill starting at z=−1 puts the box's cap at t=1. The
        // assertions read in world z, which is the vocabulary the fixtures are written in.
        let ts = (0..geom.len())
            .filter_map(|c| {
                let t = param_opt(&jd, c, &rows[0].def)?;
                // The class's world z, via the axis point at that parameter.
                let (o, m) = (rows[0].def.origin(), rows[0].def.dir());
                let z = o[2].checked_add(t.checked_mul(m[2])?)?;
                Some((c, z.to_f64()))
            })
            .collect();
        (out, ts)
    }

    /// A band's two ends as world `z` — the assertion vocabulary.
    fn ends(lf: &LocalFace, ts: &[(usize, f64)]) -> (f64, f64) {
        let Bound::Band { lo, hi } = lf.outer else {
            panic!("a band face bounds a band");
        };
        let at = |c: usize| ts.iter().find(|(k, _)| *k == c).expect("a ⊥ class").1;
        (at(lo), at(hi))
    }

    fn box_and_drill(z0: f64, h: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cylinder(
            Point3::from_array([1.0, 1.0, z0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            h,
        );
        m.rebuild_adjacency();
        (m, a, b)
    }

    /// **A through hole.** The drill runs from z=−1 to z=3 through a `[0,2]³` box; `Cut` keeps
    /// exactly the wall between the box's two caps, and the material there is **outside** the
    /// cylinder, so the face is flipped against the stored (outward) normal.
    #[test]
    fn a_through_hole_keeps_the_band_between_the_caps_and_flips_it() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Cut);
        assert_eq!(out.len(), 1, "one surviving band: {out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
        assert!(out[0].flip, "a hole's wall faces its own axis");
        assert!(matches!(out[0].surf, ClassIx::Cyl(0)));
    }

    /// The same drill, fused instead: the material is **inside** the cylinder exactly where the
    /// box is not, so the two protruding stretches survive and the buried one does not — and
    /// neither survivor is flipped.
    #[test]
    fn a_fused_drill_keeps_the_two_protruding_bands_unflipped() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Fuse);
        let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
        spans.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        assert_eq!(spans, vec![(-1.0, 0.0), (2.0, 3.0)], "{out:?}");
        assert!(out.iter().all(|lf| !lf.flip), "a boss's wall faces outward");
    }

    /// `Common` keeps the buried stretch and nothing else — the mirror of the fuse case, and the
    /// one that would pass if the keep rule ignored the operand order.
    #[test]
    fn a_common_keeps_only_the_buried_band() {
        let (m, a, b) = box_and_drill(-1.0, 4.0);
        let (out, ts) = bands(&m, a, b, BoolKind::Common);
        assert_eq!(out.len(), 1, "{out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 2.0));
        assert!(!out[0].flip, "the kept material is inside the wall");
    }

    /// **A blind hole.** The drill stops at z=1 inside the box: the wall survives from the box's
    /// cap down to the drill's own cap, which is where the bottom disk closes it.
    #[test]
    fn a_blind_hole_stops_at_the_drills_own_cap() {
        let (m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
        let (out, ts) = bands(&m, a, b, BoolKind::Cut);
        assert_eq!(out.len(), 1, "{out:?}");
        assert_eq!(ends(&out[0], &ts), (0.0, 1.0));
        assert!(out[0].flip);
    }

    /// ★★ **The first end-to-end cylinder boolean.** A `[0,2]³` box drilled through by a
    /// radius-0.5 bore: seven faces (four walls, two drilled caps, one hole wall), genus 1, and a
    /// volume of `8 − π·0.25·2`. The volume is the **winding lock** — props integrates by the
    /// divergence theorem, so a hole loop wound the wrong way returns `8 + πr²h` instead.
    #[test]
    fn a_through_hole_is_built_and_measures_what_it_should() {
        let (mut m, a, b) = box_and_drill(-1.0, 4.0);
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the drill cuts");
        assert_eq!(out.len(), 1, "one body");
        let s = out[0];
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let faces = crate::planes::solid_shell_handles(&m, s)
            .into_iter()
            .map(|sh| m.shells.get(sh).faces.len())
            .sum::<usize>();
        assert_eq!(faces, 7, "4 walls + 2 drilled caps + the bore's wall");
        // ★ **Genus 1, re-derived rather than quoted.** An earlier plan wrote `V10−E15+F7−L2`
        // from memory; the counts are measured here and what is asserted is the relation
        // (`χ = V − E + F − L = 2(S − G)`, so one shell with one through hole gives `χ = 0`).
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            0,
            "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        let v = nacre_props::mass_props(&m, s)
            .expect("a closed solid has mass props")
            .volume;
        let want = 8.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (v - want).abs() < 1e-6,
            "volume {v} vs {want} — a bore removes material, it does not add it"
        );
        // The same box without the bore, for the sign of the correction rather than its value.
        let mut m2 = Model::new();
        let plain = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        m2.rebuild_adjacency();
        let plain_v = nacre_props::mass_props(&m2, plain)
            .expect("a box has mass props")
            .volume;
        assert!(v < plain_v, "the bore removes material: {v} vs {plain_v}");
    }

    /// ★★ **Two bores in one plate — the case the band pass was rebuilt for.**
    ///
    /// The second cut's counterpart already carries a cylinder face (the first bore's wall), and
    /// the old witness road could not describe such a body, so this used to decline by name. The
    /// band pass now reads the arrangement's own disk labels, and a curved counterpart is no
    /// longer a question anyone has to answer.
    ///
    /// ★ It is also where the wall's **own** membership stopped being assumed: the first bore's
    /// wall bounds the *plate*, and inside that circle the plate has no material — the opposite of
    /// a drill, which fills its own cylinder. Assuming the drill's case (as the pass did while
    /// only drills were tested) puts the second bore's `keep` on the wrong chamber.
    #[test]
    fn a_two_hole_plate_drills_both() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 2.0, 1.0]),
        );
        let d1 = m.add_cylinder(
            Point3::from_array([1.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        let d2 = m.add_cylinder(
            Point3::from_array([3.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let one = crate::boolean(&mut m, BoolKind::Cut, a, d1).expect("first bore");
        let v1 = nacre_props::mass_props(&m, one[0]).expect("props").volume;
        assert!(
            (v1 - (8.0 - std::f64::consts::PI * 0.0625)).abs() < 1e-9,
            "one bore: {v1}"
        );
        let two = crate::boolean(&mut m, BoolKind::Cut, one[0], d2).expect("second bore");
        assert_eq!(two.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v2 = nacre_props::mass_props(&m, two[0]).expect("props").volume;
        let want = 8.0 - 2.0 * std::f64::consts::PI * 0.0625;
        assert!((v2 - want).abs() < 1e-9, "two bores: {v2} vs {want}");
        // Genus 2: two through holes in one shell.
        let reach = m.reachable();
        let l: i64 = reach
            .faces
            .iter()
            .map(|fh| m.faces.get(*fh).inner.len() as i64)
            .sum();
        let chi =
            reach.vertices.len() as i64 - reach.edges.len() as i64 + reach.faces.len() as i64 - l;
        assert_eq!(chi, -2, "two handles: χ = 2(1 − 2)");
    }

    /// **Three bores.** Two was the case the label route was built for; three is the check that
    /// nothing in it counts to two — each cut's counterpart carries one more cylinder face than
    /// the last.
    #[test]
    fn three_bores_in_a_row() {
        let mut m = Model::new();
        let mut solid = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([6.0, 2.0, 1.0]),
        );
        for x in [1.0, 3.0, 5.0] {
            let d = m.add_cylinder(
                Point3::from_array([x, 1.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.25,
                3.0,
            );
            m.rebuild_adjacency();
            solid = crate::boolean(&mut m, BoolKind::Cut, solid, d).expect("a bore")[0];
        }
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, solid).expect("props").volume;
        let want = 12.0 - 3.0 * std::f64::consts::PI * 0.0625;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **A drilled plate cut by a body flush with its faces.** `Common` of a drilled plate with a
    /// box covering half of it: the two operands share the plate's top and bottom planes, and the
    /// bore's rims end on *holed* faces rather than on disks.
    #[test]
    fn a_drilled_plate_can_be_cut_by_a_body_flush_with_its_faces() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 2.0, 1.0]),
        );
        let d = m.add_cylinder(
            Point3::from_array([1.0, 1.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let drilled = crate::boolean(&mut m, BoolKind::Cut, plate, d).expect("bore")[0];
        let half = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Common, drilled, half).expect("flush common");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4.0 - std::f64::consts::PI * 0.0625;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **A cylinder standing far away, whose caps happen to land on the box's own planes.** Two
    /// coplanar planes are one class whichever way they sit, so this cylinder's *disks* share a
    /// class with the box's faces while touching nothing — the shape that made the retired seated
    /// rule mistake a shared class for a contact. `Cut` takes nothing away.
    #[test]
    fn a_distant_cap_on_the_boxs_own_plane_removes_nothing() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        // Far from the box in x/y, but its caps land exactly on z = 0 and z = 2.
        let b = m.add_cylinder(
            Point3::from_array([10.0, 10.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
        assert_eq!(out.len(), 1);
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");
    }

    // ---- Seated caps: a cap flush on the other body's face (M6-2 finish) ----
    //
    // ★★ These five are the population the `SeatedCylinderCap` rule refused. What actually makes a
    // seated circle hard is its boundary meeting the counterpart's — and a boundary is either an
    // edge on a plane (parallel to the axis → `WallMeetsLateral`, oblique → `ObliqueCylinderCut`)
    // or another cylinder's rim (→ `CylinderPairContact`). The two fences at the end of this block
    // are those names still standing, so this file says both what was opened and what was not.
    //
    // ★★★ **Measured against the pre-deletion kernel, all seven of these came back
    // `SeatedCylinderCap` — the two fences included.** The rule stood before the wall and the
    // curved-depth rules and answered for them: an overhanging boss and a pair of coplanar-capped
    // cylinders are refused for reasons that have nothing to do with seating, and the sentence a
    // user got named the seating anyway. Deleting it does not only open the five; it lets the two
    // that stay shut say what is actually in the way.

    /// **(a) A through hole whose two caps are flush with the plate's own faces.** The drill neither
    /// overshoots nor stops short: both cap planes coincide with the plate's, and every ⊥ class in
    /// the operation carries faces of both operands.
    #[test]
    fn a_flush_through_drill_bores_the_plate() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a flush through hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(b) A blind hole seated on the face it is drilled from.** Only the lower cap is flush;
    /// the upper one stops inside the plate and closes the bore itself.
    #[test]
    fn a_blind_drill_seated_on_the_plates_base_bores_it() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a seated blind hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 - std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(c) A boss standing on the plate.** `Fuse` with the cylinder's base cap flush on the
    /// plate's top face — the seating a person draws first, and the one the retired rule refused.
    #[test]
    fn a_boss_seated_on_the_plate_fuses_to_it() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a seated boss");
        assert_eq!(out.len(), 1, "one body, not two");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(d) `Common` where the cylinder's caps sit on the box's own faces.** The intersection is
    /// the whole cylinder — every one of its faces is seated or shared.
    #[test]
    fn a_flush_cylinder_meets_the_box_in_itself() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b = m.add_cylinder(
            Point3::from_array([1.0, 1.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Common, a, b).expect("a flush common");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **(i) Drilling the floor of a pocket.** The seated face here is not an original operand
    /// face at all — it is a *floor the previous boolean made*, so the seating arrives through the
    /// arrangement rather than from the modeller. A blind bore down from that floor.
    #[test]
    fn a_drill_seated_on_a_pocket_floor_bores_it() {
        let mut m = Model::new();
        let block = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([6.0; 3]));
        // Open at the top — a tool that stopped inside would leave a void, not a pocket.
        let pocket_tool = m.add_cuboid(
            Point3::from_array([1.0, 1.0, 2.0]),
            Point3::from_array([5.0, 5.0, 7.0]),
        );
        m.rebuild_adjacency();
        let pocketed =
            crate::boolean(&mut m, BoolKind::Cut, block, pocket_tool).expect("pocket")[0];
        // The bore hangs from the pocket floor (z = 2) down into the material below it.
        let drill = m.add_cylinder(
            Point3::from_array([3.0, 3.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, pocketed, drill).expect("a floor bore");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 216.0 - 64.0 - std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// One rulings-road boolean on the 40×40×20 plate, end to end on the production road:
    /// build, retire the operands, validate clean, answer the exact volume, tessellate
    /// watertight. The volume oracle is what finally *measures* the ladder's sign roster —
    /// a flipped panel winding, disk-side selector, ruling turn or chord sense moves it.
    fn through_boss_builds(kind: BoolKind, base: [f64; 3], want: f64) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            50.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, kind, plate, boss).expect("the rulings road builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the panels tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
    }

    /// ★ **The through-boss builds** — the rulings ladder's milestone: a boss standing through
    /// the plate's wall (axis exactly on the `x = 40` plane), the app-measured refusal that
    /// opened M6-2's remainder. Fuse keeps the outer half (plate + the boss outside it), cut
    /// carves the notch, common keeps the inner half-cylinder — each an exact closed form, and
    /// the three exercise both complementary θ-sectors.
    #[test]
    fn a_through_boss_fuses() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Fuse, [40.0, 20.0, -10.0], 32000.0 + 1000.0 * pi);
    }

    #[test]
    fn a_through_boss_cuts_a_notch() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Cut, [40.0, 20.0, -10.0], 32000.0 - 250.0 * pi);
    }

    #[test]
    fn a_through_boss_common_is_the_inner_half() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Common, [40.0, 20.0, -10.0], 250.0 * pi);
    }

    /// ★ A boss on the plate's **corner** — its axis on *two* wall planes, both recorded, each
    /// circle cut by both walls: the rim pairing (`(wall class, root)`) and the θ-panel
    /// machinery quarter the lateral. Measured to assemble in the lift probe; the volume pins
    /// it (plate + three quarters of the cylinder outside).
    #[test]
    fn a_corner_boss_fuses() {
        let pi = std::f64::consts::PI;
        through_boss_builds(BoolKind::Fuse, [40.0, 40.0, -10.0], 32000.0 + 1125.0 * pi);
    }

    /// ★ The populations the gate's arm deliberately keeps out, each measured lifted before the
    /// arm was shaped: a **tangent** wall (distance exactly `r`) assembles a volume-correct
    /// zero-thickness pinch `validate` cannot see; an **offset** crossing (`0 <` distance `< r`)
    /// walks to `OpenResultShell`, a SuspectedDefect label an honest input must not wear. The
    /// live set survives every refusal.
    ///
    /// ★ A third row lived here until the D2b cutover: the **half-height** boss, whose upper cap
    /// sits inside the plate's material — the band road's `chamber` had no sector answer for that
    /// end and refused it `RulingBoundNotYet`. The chart reads it (a band below the plate, the
    /// outer sector beside it), so it builds now — [`Self::a_half_height_boss_builds`].
    #[test]
    fn the_gate_still_refuses_what_the_road_does_not_serve() {
        for (base, h, want) in [
            ([38.0, 20.0, -10.0], 50.0, RejectReason::WallMeetsLateral),
            ([35.0, 20.0, -10.0], 50.0, RejectReason::WallMeetsLateral),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 20.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(base),
                Vector3::from_array([0.0, 0.0, 1.0]),
                5.0,
                h,
            );
            m.rebuild_adjacency();
            let live = m.live_solids.clone();
            let err =
                crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect_err("outside the road");
            let BoolError::Rejected { reason, .. } = err else {
                panic!("a rejection, not {err:?}");
            };
            assert_eq!(reason, want, "{base:?} h {h}");
            assert_eq!(m.live_solids, live, "the live set survives the refusal");
        }
    }

    /// ★★★★★ **The half-height boss builds — the capability the D2b cutover opened.** A boss
    /// through the plate's wall whose upper cap (z = 10) sits inside the plate: below the plate
    /// the lateral is a whole band, beside it only the outer sector survives, and the cap is a
    /// half-disk. The band road refused this end (`RulingBoundNotYet`); the chart's cells read
    /// it off the same labels. Volumes derived, not copied: the boss outside the plate is the
    /// outer half-cylinder over the full height (`π·25·20/2 = 250π`) plus the inner half below
    /// the plate (`π·25·10/2 = 125π`); the boss inside the plate is the inner half over
    /// `z ∈ [0, 10]` (`125π`).
    ///
    /// The lateral stays **two** faces (band + outer ring): `merge_curved_group` abstains when
    /// a group has a rim of only one sense (the cut rim at z = 0 is refined to arcs, the outer
    /// ring's lower arc erases one, and no `Bound` spells "one rim and a chain of arcs").
    #[test]
    fn a_half_height_boss_builds() {
        let pi = std::f64::consts::PI;
        for (kind, want, lateral) in [
            (BoolKind::Fuse, 32000.0 + 375.0 * pi, vec![2]),
            (BoolKind::Cut, 32000.0 - 125.0 * pi, vec![1]),
            (BoolKind::Common, 125.0 * pi, vec![1]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 20.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array([40.0, 20.0, -10.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                5.0,
                20.0,
            );
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, kind, plate, boss)
                .unwrap_or_else(|e| panic!("{kind:?}: the half-height boss builds: {e:?}"));
            assert_eq!(out.len(), 1, "{kind:?}: one solid");
            m.rebuild_adjacency();
            assert_eq!(m.live_solids, out, "{kind:?}: the operands retired");
            let issues = nacre_validate::validate(&m);
            assert!(issues.is_empty(), "{kind:?}: {issues:?}");
            let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
            assert!((v - want).abs() <= 1e-9 * want, "{kind:?}: {v} vs {want}");
            assert_eq!(
                lateral_face_counts(&m, out[0]),
                lateral,
                "{kind:?}: lateral faces"
            );
            let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
                .expect("the half-height boss tessellates");
            let mut uses: std::collections::HashMap<(u32, u32), usize> =
                std::collections::HashMap::new();
            for (_, tri) in mesh.triangles.iter() {
                for k in 0..3 {
                    let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                    *uses.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            assert_eq!(
                uses.values().filter(|&&n| n != 2).count(),
                0,
                "{kind:?}: the mesh is watertight"
            );
        }
    }

    /// **The straddling boss builds.** The M6-2b milestone fence: a boss hanging over the
    /// plate's edge — its rim circle cut by the plate top's boundary segment — fuses into one
    /// valid solid. The ladder this closes, in the names its rungs wore: `WallMeetsLateral`
    /// (the wall rule read the closed span) → `CircleMeetsSegment` → `ArcBoundNotYet` (the
    /// stopper, walked from the class arrangement to the very end of the assembly) → built.
    ///
    /// The volume and the mixed-loop integrals are pinned by
    /// [`Self::a_cut_rim_boolean_builds_a_complete_solid`]; here the result answers the two
    /// whole-model judges: `validate` (whose winding check reads the arcs as witnesses now) and
    /// the tessellation (watertight, arcs sampled as sub-arcs).
    #[test]
    fn a_boss_overhanging_the_plates_edge_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        assert_eq!(out.len(), 1, "one fused solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the arcs tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "the mesh is watertight");
    }

    /// ★★ **A boss standing over a bore, its whole outline *inside* the rim.** No segment crosses
    /// the circle, so the circle keeps its closed cell — what the shape really asks is that the
    /// **disk cell host a polygon**, and `nest_cells` refused to (a `debug_assert` said a disk
    /// hosts nothing, for a reason that was already wrong: the guard actually keeping the
    /// population out was `circles_meet_no_segment`).
    ///
    /// ★★★ **The population immediately found a live defect in the nesting predicate.** With the
    /// disk allowed to host, the *circle*'s contour came back as a hole of the **square** — the
    /// footprint `[7,9]×[9,11]` straddles the axis `(8,10)`, so "the circle's centre is inside the
    /// ring" is true in the nesting that does not hold. One witness point says *whether* two
    /// disjoint loops nest, never *which way*; the other direction is what settles it
    /// (`cell_in_cell`).
    #[test]
    fn a_segment_inside_the_rim_builds_the_boss_over_the_hole() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // Footprint `[7,9] × [9,11]` sits wholly inside the rim `(8,10)`, `r = 3`; `z ∈ [5,8]`
        // leaves the bore's span `t ∈ [1,6]` clear along the axis, so the wall rule passes it on.
        let boss = m.add_cuboid(
            Point3::from_array([7.0, 9.0, 5.0]),
            Point3::from_array([9.0, 11.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a boss over a bore");
        // ★ The boss's whole footprint is inside the rim, so it stands over the **hole** and
        // touches no material: two bodies, not one. That is forced rather than chosen — a boss
        // that reached material inside the disk would lie within the bore's axial span, and the
        // wall rule refuses that as `WallMeetsLateral`.
        assert_eq!(out.len(), 2, "the boss touches nothing");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [
            2.0 * 2.0 * 3.0,
            40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0,
        ];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "the boss and the bored plate: {vols:?} vs {want:?}"
        );
        // The plate keeps its through hole (χ = 0), the boss is a box (χ = 2).
        let mut chi: Vec<i64> = out
            .iter()
            .map(|&s| {
                let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
                v_n - e_n + f_n - l_n
            })
            .collect();
        chi.sort_unstable();
        assert_eq!(chi, vec![0, 2], "genus 1 and genus 0");
    }

    /// **A tangency builds.** The boss's `x = 11` edge is exactly tangent to the rim `(8,10)`,
    /// `r = 3`, so the locator's quadratic has a double root — the third arm of `CylinderMeet`
    /// that can reach here, and the only one whose point carries no radical (`(11, 10, 5)`,
    /// rational).
    ///
    /// ★★ **And a touch is not a break**: the circle keeps its closed cell, so nothing about the
    /// arrangement had to change for this to work. It was refused only because the guard's
    /// sentence was wider than its proposition — measured before the narrowing landed, with the
    /// whole guard off, this already produced exactly the solid asserted below.
    #[test]
    fn a_segment_tangent_to_the_rim_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let boss = m.add_cuboid(
            Point3::from_array([11.0, 8.0, 5.0]),
            Point3::from_array([15.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("a tangent edge");
        assert_eq!(out.len(), 1, "the boss sits on material");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 40.0 * 20.0 * 5.0 - std::f64::consts::PI * 9.0 * 5.0 + 4.0 * 4.0 * 3.0;
        assert!(
            (v - want).abs() < 1e-9,
            "plate less bore plus boss: {v} vs {want}"
        );
        // The through bore survives the fuse, so the one shell still has one handle.
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            0,
            "genus 1: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        // ★★★★★ **And the solid cannot be meshed — say so here, where the geometry is.**
        // The boss's footprint touches the bore's rim at exactly one point, so the plate's top
        // face has two inner loops meeting there: its interior is pinched, and no triangulation
        // of a pinched region exists in this decomposition. That is a statement about the
        // *tessellator*, not about this solid — `validate` is clean above and the volume is
        // exact — and `TessError`'s own doc draws that line.
        // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
        // census that only agreed with itself would be measuring when it looks rather than what
        // came out. Both spellings say the same thing here.
        m.rebuild_adjacency();
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default());
        assert!(
            matches!(mesh, Err(nacre_tess::TessError::SelfTouchingBoundary)),
            "the tangency's mesh: {:?}",
            mesh.map(|t| t.triangles.len())
        );
    }

    /// ★★ **The fence post: a crossing that lands exactly on a segment's endpoint.** The boss's
    /// corner `(8,13)` sits on the rim `(8,10)`, `r = 3`, so both edges leaving it meet the circle
    /// *at* their own end — the one place the in-segment test's inequality can be open or closed,
    /// and the two spellings give different answers.
    ///
    /// **Measured, both ways:** counting `Zero` as inside names the corner `(8, 13, 5)`; excluding
    /// it drops every root and falls back to the containment witness — a `Segment`, which would
    /// say "this edge lies inside the circle" about an edge that runs *outward* from a single
    /// touching point. Inclusive is the true sentence, and it is the one the arc split will need.
    #[test]
    fn a_crossing_on_a_segments_endpoint_is_inside_it() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let boss = m.add_cuboid(
            Point3::from_array([8.0, 13.0, 5.0]),
            Point3::from_array([12.0, 17.0, 8.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect_err("corner on rim");
        // ★★ **One point, two names** — and the split says so by name. `(8,13,5)` is the boss's
        // corner, so the arrangement already holds it as a three-plane vertex; the rim crossing
        // names the *same* point as a `NodeId::Branch`. The DCEL keys vertices by name, so shipping
        // both would put two vertices where there is one — folding them is its own step, and until
        // then `CoincidentNodes` ("two names for one point") is the true sentence.
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CoincidentNodes,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★★ **The population where the naming rule actually fires — and the first non-`+Z`
    /// cylinder in this repository.**
    ///
    /// Every `add_cylinder*` call in the suite states the axis `[0, 0, 1]` (measured: 87 of them),
    /// and that is not an accident of taste — `add_cuboid` pushes faces `[−Z, +Z, −Y, +Y, −X, +X]`
    /// and `build_prism` pushes "base cap, top cap, then walls", so a `+Z` cylinder's circle always
    /// lives on a class interned **before** the wall that crosses it. Measured over the whole ops
    /// suite: 463 of the 2060 `(class, wall)` pairs the gate examines are descending, but **all
    /// nine that reach the locator are ascending**. So `NodeId::branch`'s canonicalization —
    /// re-sorting the pair and restating the root with it — has never once been exercised by a
    /// production path.
    ///
    /// Turning the boss onto `+X` inverts it: the rim lives on the box's `x = 4` face (class 5,
    /// interned last) and the segment it crosses is on the `z = 2` cap (class 1). The pair swaps.
    ///
    /// ★ **A tangency, so there is exactly one root** — and since 2026-08-21 that also means the
    /// circle is not separated, so this **builds** rather than refusing. The single root is still
    /// what makes the naming rule visible here (it is what the fixture was written for), and the
    /// witness it once asserted is now the *split point that never happens*. That matters: the
    /// straddling fixtures deliberately
    /// accept either of their two crossings, and a fixture copied from that template would be green
    /// whether or not the root rule is right. Here `disc = 0` — the rim (`x = 4`, centre
    /// `(4, 2, 1.5)`, `r = 0.5`) touches `z = 2` at the single point `(4, 2, 2)`, solved from the
    /// fixture's own numbers — and the name it must take is `QuadRoot::Double`, which a swap must
    /// **not** toggle.
    #[test]
    fn a_turned_boss_tangent_to_the_plate_top_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, 1.5]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("a turned boss");
        assert_eq!(
            out.len(),
            1,
            "the boss's base disk sits on the plate's face"
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!(
            (v - want).abs() < 1e-9,
            "plate plus the whole boss: {v} vs {want}"
        );
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
        // ★★★★★ **And the solid cannot be meshed — the same pinch, spelled inner-to-outer.**
        // The boss's base circle is tangent to the plate's `z = 2` edge at `(4, 2, 2)`, so the
        // `x = 4` face's hole touches its own outer ring at one point and the face's interior is
        // pinched there. `validate` is clean above and the volume is exact: this names what the
        // tessellator cannot draw, not a defect in the solid.
        // ★ Rebuilt first, on purpose: `boolean`'s own census meshes *before* the rebuild, and a
        // census that only agreed with itself would be measuring when it looks rather than what
        // came out. Both spellings say the same thing here.
        m.rebuild_adjacency();
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default());
        assert!(
            matches!(mesh, Err(nacre_tess::TessError::SelfTouchingBoundary)),
            "the tangency's mesh: {:?}",
            mesh.map(|t| t.triangles.len())
        );
    }

    /// **The turned boss builds — and its crescent is the winding witness's red switch.**
    ///
    /// The near cap's circle is cut by two different plate edges (top and corner), so this
    /// result carries the arc-dominated crescent face whose chord Newell reads **backwards**
    /// (`cos = −1`, the measured wall) — `validate == []` here is what pins `loop_winding`'s
    /// segment witnesses. ★ The old fence's proposition — a root that fails to follow its pair
    /// through the sort names the wrong crossing — did not retire with the reject: the mint
    /// fence asserts the branch vertices sit **on the derived crossings**, and a wrong root
    /// moves the minted point itself.
    #[test]
    fn a_turned_boss_over_the_plates_corner_builds() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 0.25, 2.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the boss builds");
        assert_eq!(out.len(), 1, "one fused solid");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 32.0 + std::f64::consts::PI * 0.25;
        assert!((v - want).abs() < 1e-12, "{v} vs {want}");
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            2,
            "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }

    /// **The audit and the boolean say the same thing about an arc class — and it says what the
    /// arrangement produced.**
    ///
    /// ★★★★ Three things are pinned here and none had a fence before.
    ///
    /// **One: the two copies of the per-class pipeline agree.** `frame_audit` re-runs pass A and B
    /// inline rather than calling `arrange`, so a split hoisted into one and not the other compiles
    /// perfectly and makes the audit report *no failure* for an input the boolean refuses —
    /// `decline_to_reject`'s doc calls that "the worst possible time to be lying". The existing
    /// guard (`the_audit_does_not_invent_failures`) cannot see it: its fixture is deliberately one
    /// that fails **outside** the class pipeline.
    ///
    /// **Two: what the arrangement produced on an arc class.** The stopper stands after every
    /// arrangement stage and swallows their answer, so without this the whole arc population is
    /// measured by a probe and the probe is deleted before the commit. The numbers are **derived,
    /// not read back**: the
    /// boss's rim crosses one edge of the plate's top face twice, cutting the circle into **2**
    /// arcs and the plane into **4** cells — `rect−disk`, `rect∩disk`, `disk−rect` and the outside.
    /// Only the outside winds `−1`, so there is **1** root group; and it shares nodes with all
    /// three `+1` cells (the plate's corners, the two crossings), so `cell_in_cell` answers "not
    /// comparable" every time and there are **0** holes. Exactly **one** class is cut: the boss
    /// stands *on* the top face, so only its base circle lies in a plane of the plate.
    ///
    /// **Three: the faces it emitted, as coordinates** — and this is the one that sees a wrong
    /// `sense`. Flipping the sense the split carries onto its sub-segments attaches the arcs to
    /// the wrong cells, and **every number in part Two stays put** (the cells still read four with
    /// one `−1`, the roots one, the holes zero, the sorted labels identical) — that blindness was
    /// derived before it was measured, and it held for four rungs. The emitted ring is where it
    /// finally shows, as the **exact reverse**, which is why rotation is normalized here and
    /// reversal is not. See `ClassAudit::outer_rings` for both arguments and for why the red probe
    /// is read on the **turned** boss: the straddling one is blind to the flip *and* its 2-node
    /// ring cannot express a reversal at all.
    ///
    /// ★ **Two earlier claims here were wrong and are recorded rather than quietly dropped**: that
    /// this fence is *green* under that flip (it is red — measured, three times), and that the
    /// flip's first reader is `label_cells`' keep decision (refuted a rung earlier: every
    /// order-independent summary of the labels is identical; it is `emit_faces`).
    #[test]
    fn the_audit_and_the_boolean_agree_about_an_arc_class() {
        // `y = 0` meets the turned boss's circle (centre `(y,z) = (0.25, 2)`, `r = 0.5`) at
        // `(z - 2)² = 0.1875`. The one irrational coordinate in either fixture, and the reason the
        // ring comparison is `near()` rather than `==`.
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        for (origin, axis, n_out, want_rings) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                // The plate's top face carries this class, so the ring is CCW about `+z`.
                [0.0, 0.0, 1.0],
                vec![
                    // `rect - disk`: the plate's top, three corners and the inner arc.
                    vec![
                        [0.0, 0.0, 2.0],
                        [4.0, 0.0, 2.0],
                        [4.0, 1.5, 2.0],
                        [4.0, 2.5, 2.0],
                        [4.0, 4.0, 2.0],
                        [0.0, 4.0, 2.0],
                    ],
                    // `disk - rect`: the overhang's underside — chord and outer arc.
                    vec![[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
                ],
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                // The plate's `x = 4` face carries this one, so CCW is about `+x`.
                [1.0, 0.0, 0.0],
                vec![
                    // `rect - disk`: the plate's right face, the disk biting its top-left corner.
                    vec![
                        [4.0, 0.0, 0.0],
                        [4.0, 4.0, 0.0],
                        [4.0, 4.0, 2.0],
                        [4.0, 0.75, 2.0],
                        [4.0, 0.0, s],
                    ],
                    // `disk - rect`: the boss's base cap outside the plate. The corner `(4,0,2)` is
                    // `0.25` from the circle's centre, so it sits *inside* the disk and is one of
                    // this ring's three nodes.
                    vec![[4.0, 0.0, s], [4.0, 0.0, 2.0], [4.0, 0.75, 2.0]],
                ],
            ),
        ] {
            // The boolean takes `&mut Model` and the audit `&Model`; a rejected boolean leaves
            // arena residue, so each gets its own build of the same fixture.
            let build = || {
                let mut m = Model::new();
                let plate = m.add_cuboid(
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([4.0, 4.0, 2.0]),
                );
                let boss = m.add_cylinder(
                    Point3::from_array(origin),
                    Vector3::from_array(axis),
                    0.5,
                    1.0,
                );
                m.rebuild_adjacency();
                (m, plate, boss)
            };
            let (mut m, plate, boss) = build();
            crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the arc class assembles");
            let (m, plate, boss) = build();
            let audits = crate::arrangement::frame_audit(&m, BoolKind::Fuse, plate, boss).unwrap();
            // The agreement, now that the population is green: the boolean built, and the audit's
            // replay of the same pipeline stops nowhere either.
            let stopped: Vec<_> = audits.iter().filter(|a| a.failed_at.is_some()).collect();
            assert!(stopped.is_empty(), "no class stops: {stopped:?}");
            // ★ The whole audit, not just `produced`: the ring coordinates are a sibling field, and
            // joining two filtered lists by index is the seam where they could come from different
            // classes. ★ By reference: `Produced` carries the labels now, so it is no longer `Copy`.
            let cut: Vec<_> = audits
                .iter()
                .filter(|a| a.produced.as_ref().is_some_and(|p| p.arcs > 0))
                .collect();
            assert_eq!(cut.len(), 1, "exactly one class has its circle cut");

            assert_eq!(
                *cut[0].produced.as_ref().unwrap(),
                crate::arrangement::Produced {
                    cells: 4,
                    arcs: 2,
                    roots: 1,
                    holes: 0,
                    // ★★ **Derived from the geometry, not read back.** `[A_above, A_below,
                    // B_above, B_below]` for the three `+1` regions, sorted: the boss stands *on*
                    // the plate, so plate-minus-disk has the plate below it and nothing above;
                    // disk-minus-plate — the overhang's underside — has the boss above it and
                    // nothing below; and their intersection has both. The outside cell winds `-1`
                    // and is not here.
                    //
                    // ★ **The same three come out for the turned boss, and the story is the same
                    // one in its own frame** — this assertion is inside the two-fixture loop, so
                    // reading it as a sentence about the straddling boss alone would be reading
                    // half of what it checks. There the class is `x = 4`: the plate lies on the
                    // `−x` side, the boss's base cap outside it on `+x`, and the overlap has both.
                    pos_labels: vec![
                        [false, false, true, false],
                        [false, true, false, false],
                        [false, true, true, false],
                    ],
                },
                "the arrangement walked the arcs, nested them, and labelled them"
            );
            // ★★★★★ **And the faces it emitted, as coordinates.** This is what a wrong `sense` on a
            // split segment moves and nothing before it does: the cell counts, the nesting and the
            // sorted labels above are all identical whichever way the rings run.
            //
            // ★★★ **The premise first.** `emit_faces` emits CCW about `n_out(wc)`, and which face
            // is the class root — hence which way `n_out` points — is a plane-table fact, not one
            // of the fixture's numbers. Pinning it here means a root flip fails *this* assertion,
            // with its own sentence, instead of silently reversing every ring below. Both come
            // from the plane table, which no arrangement stage can move.
            let got_n_out: Vec<f64> = cut[0]
                .root_normal
                .iter()
                .map(|c| c * cut[0].orient_sign as f64)
                .collect();
            assert!(
                got_n_out
                    .iter()
                    .zip(&n_out)
                    .all(|(a, b)| (a - b).abs() < 1e-9),
                "the rings below are derived CCW about {n_out:?}, but the class faces {got_n_out:?}"
            );
            let got = cut[0]
                .outer_rings
                .as_ref()
                .expect("the per-class road emits the rings");
            assert_eq!(
                got.len(),
                want_rings.len(),
                "one outer ring per emitted face"
            );
            for (g, w) in got.iter().zip(&want_rings) {
                assert_eq!(g.len(), w.len(), "ring length: got {g:?}, want {w:?}");
                assert!(
                    g.iter()
                        .zip(w)
                        .all(|(p, q)| (0..3).all(|i| (p[i] - q[i]).abs() < 1e-9)),
                    "ring: got {g:?}, want {w:?}"
                );
            }
        }
    }

    /// **The seam realizes a branch vertex and measures it — seen through the door production
    /// uses, because nothing else can see it at all.**
    ///
    /// ★★★ The deferred stopper intercepts the whole seam stretch, so with the branch arm
    /// deleted the seam fails `BranchVertexUnnamed`, the interception swallows it, and **every
    /// boolean-level fence stays green** (measured — that probe is what forced this test). The
    /// arm's only witness is a direct second consumer of `seam_table` on the very faces
    /// production feeds it, so this walks production's stretch step for step: trace, clean,
    /// append the bands, build the seam.
    ///
    /// ★ The tolerance is asserted as a **bound**, never a copied value; the branch coordinates
    /// themselves are pinned by `ClassAudit::outer_rings` through the same `branch_point` road,
    /// so re-asserting them here would be a second copy of an existing lock — and a `Lo`/`Hi`
    /// mix-up cannot hide behind the bound either, because both crossings lie on every defining
    /// surface and `outer_rings` is what tells them apart.
    #[test]
    fn the_seam_realizes_a_branch_vertex_and_measures_it() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
            let PlaneSetup {
                planes: faces_tab,
                geom,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                class_owner,
                n_a,
                standard,
                notes,
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, deferred) =
                crate::arrangement::trace_result_faces_full_for_test(
                    &m,
                    BoolKind::Fuse,
                    plate,
                    boss,
                    &jd,
                    faces_tab,
                    plane_ix,
                    cyls,
                    *n_a,
                    class_owner,
                    &trace_in,
                )
                .expect("the arc population traces");
            // The stopper socket is empty since the population went green — nothing defers.
            assert!(deferred.is_none(), "{deferred:?}");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                band_faces(
                    BoolKind::Fuse,
                    &faces,
                    &rows,
                    &jd,
                    &curved.disk_labels,
                    &curved.cut_rims,
                    &curved.arc_labels,
                )
                .expect("bands"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
                .expect("the seam realizes branch nodes");
            let branch: Vec<_> = seam
                .iter()
                .filter(|sv| crate::combinatorics::branch_name(sv.triple).is_some())
                .collect();
            assert_eq!(branch.len(), 2, "both crossings reach the seam, once each");
            for sv in branch {
                assert!(
                    sv.tol < 1e-12,
                    "a branch realization sits on everything that defines it: tol {}",
                    sv.tol
                );
            }
        }
    }

    /// **Every ring node of the arc population's result faces earns a definition — measured
    /// through the door production uses.**
    ///
    /// ★★★ Same instrument shape as the seam fence, same reason: the deferred stopper stands
    /// behind `name_result_vertices` (at the assembly's very end now) and intercepts
    /// everything, so no reject name can testify the naming completed — a direct second consumer
    /// is the only witness. This walks production's road (trace → clean → bands → seam → naming)
    /// and asserts on its product.
    ///
    /// ★★ Disabling the branch arm reddens both fixtures (and the walls-fallback cannot fake a
    /// branch def past its arc-carrier guard). ★ The fallback itself went **zero-population** when
    /// the split-twin subdivision landed — the bitten corner's twins match now, so its def comes
    /// down the far-plane road and no probe reddens on the fallback alone; the subdivision has
    /// its own fence (`a_subdivided_twin_matches_its_neighbour_edge_for_edge`).
    #[test]
    fn every_result_vertex_of_the_arc_population_is_named() {
        for (origin, axis, bites_corner) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], false),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], true),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
            let PlaneSetup {
                planes: faces_tab,
                geom,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                class_owner,
                n_a,
                standard,
                notes,
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                band_faces(
                    BoolKind::Fuse,
                    &faces,
                    &rows,
                    &jd,
                    &curved.disk_labels,
                    &curved.cut_rims,
                    &curved.arc_labels,
                )
                .expect("bands"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming stages run");
            let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
            // Completeness: every ring node has a definition.
            let mut missing = Vec::new();
            for (fi, lf) in live.iter().enumerate() {
                for ring in lf.poly_rings() {
                    for &n in ring.iter() {
                        if !named.defs.contains_key(&(named.group_of[fi], n)) {
                            missing.push(n);
                        }
                    }
                }
            }
            assert!(missing.is_empty(), "def-less nodes: {missing:?}");
            let branch_defs = named
                .defs
                .values()
                .filter(|d| matches!(d, crate::boolean::Def::Branch { .. }))
                .count();
            assert_eq!(branch_defs, 2, "both crossings are declared, once each");
            // The bitten corner's def names the right point: realize its three planes and land
            // on (4, 0, 2) — the fixture's own number, no class index copied.
            if bites_corner {
                let hit = named.defs.values().any(|d| {
                    let crate::boolean::Def::Three(t) = d else {
                        return false;
                    };
                    nacre_geom::intersect::three_planes(
                        &geom[t[0]].plane,
                        &geom[t[1]].plane,
                        &geom[t[2]].plane,
                    )
                    .is_some_and(|p| {
                        let c = p.as_array();
                        (0..3).all(|i| (c[i] - [4.0, 0.0, 2.0][i]).abs() < 1e-9)
                    })
                });
                assert!(hit, "the bitten corner's def realizes to (4, 0, 2)");
            }
        }
    }

    /// **After the split-twin subdivision, every segment edge has exactly one twin.**
    ///
    /// ★★★ The subdivision cannot be locked through the naming fence — the walls-fallback
    /// rescues the bitten corner whether or not the twins match, so that fence is green either
    /// way. What the subdivision actually changes is the **edge-key census**: a neighbour's whole
    /// edge and the arc class's subdivided pieces share `norm_edge` keys only once the whole edge
    /// is cut at the same branch nodes. So the proposition, with its one exception stated:
    ///
    /// > every segment ring edge (a `Wall::Plane` carrier) has its key used by exactly two
    /// > faces.
    ///
    /// Arc edges are excluded because their far side is the band, which contributes no ring.
    /// ★ Both-ends-branch keys used to be excluded too — the chord and the two arcs between one
    /// branch pair folded into a single `norm_edge` key — but the carrier gave arcs their own
    /// ordered key, so the chord's line key counts exactly its two coplanar faces now (measured:
    /// the exclusion removed, both fixtures stay green — the tightening the carrier cell's plan
    /// predicted).
    ///
    /// red: with the subdivision disabled, the pre-split single-use keys come back.
    #[test]
    fn a_subdivided_twin_matches_its_neighbour_edge_for_edge() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
            let PlaneSetup {
                planes: faces_tab,
                geom,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                class_owner,
                n_a,
                standard,
                notes,
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                band_faces(
                    BoolKind::Fuse,
                    &faces,
                    &rows,
                    &jd,
                    &curved.disk_labels,
                    &curved.cut_rims,
                    &curved.arc_labels,
                )
                .expect("bands"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd).expect("seam");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming stages run");
            let live: &[LocalFace] = named.per_solid.as_deref().unwrap_or(&faces);
            let mut uses: std::collections::HashMap<
                (crate::combinatorics::NodeId, crate::combinatorics::NodeId),
                usize,
            > = std::collections::HashMap::new();
            let mut plain: Vec<(crate::combinatorics::NodeId, crate::combinatorics::NodeId)> =
                Vec::new();
            for lf in live {
                for ring in lf.poly_rings() {
                    let k = ring.nodes.len();
                    for t in 0..k {
                        let (a, b) = (ring.nodes[t], ring.nodes[(t + 1) % k]);
                        if matches!(ring.walls[t], crate::boolean::Wall::Arc { .. }) {
                            continue;
                        }
                        let key = crate::boolean::norm_edge(a, b);
                        *uses.entry(key).or_insert(0) += 1;
                        if !plain.contains(&key) {
                            plain.push(key);
                        }
                    }
                }
            }
            let odd: Vec<_> = plain
                .iter()
                .filter(|k| uses[*k] != 2)
                .map(|k| (*k, uses[k]))
                .collect();
            assert!(odd.is_empty(), "edges without exactly one twin: {odd:?}");
        }
    }

    /// **The branch vertices are minted — canonical, measured, on the derived crossings.**
    ///
    /// Before the boolean no `VertexDef::Branch` exists anywhere in the model, so a whole-store
    /// filter is position-independent; a wrong `QuadRoot` canonicalization moves the minted
    /// point itself, which is what keeps the old toggle-lock alive now that the reject (whose
    /// witness once carried it) is gone.
    ///
    /// ★ The coordinates are the fixtures' own crossing derivations (the same numbers the ring
    /// and seam fences pin) — nothing here is copied from a run. The tolerance is a bound, and
    /// ascending handle order is `VertexDef::Branch`'s own contract, minted through
    /// `QuadRoot::canonical`'s second answer.
    #[test]
    fn a_branch_vertex_is_minted_and_measured() {
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        for (origin, axis, crossings) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
            ),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let branch: Vec<_> = m
                .vertices
                .iter()
                .filter(|(_, v)| matches!(v.def, nacre_topo::VertexDef::Branch { .. }))
                .collect();
            assert_eq!(branch.len(), 2, "both crossings minted, once each");
            for (h, v) in branch {
                let nacre_topo::VertexDef::Branch { planes: [a, b], .. } = v.def else {
                    unreachable!("filtered above");
                };
                assert!(a < b, "planes in ascending handle order: {a:?} vs {b:?}");
                let p = m.vertex_point(h);
                assert!(
                    crossings
                        .iter()
                        .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                    "a minted branch vertex sits on a derived crossing: {p:?}"
                );
                let tol = m
                    .vertex_tol(h)
                    .expect("a discovered vertex carries its tol");
                assert!(tol < 1e-12, "measured against what defines it: {tol}");
            }
        }
    }

    /// **The two complementary arcs are two edges, and the cut rim is none.**
    ///
    /// ★★★ The observable form of the `[A, B]`-CCW convention (`derive_edge_curve`'s circle
    /// arm): between one pair of branch vertices a circle offers two pieces, the endpoints
    /// alone cannot tell them apart, and the *vertex order* is the bit that does — so the store
    /// must hold **two** circle-carrier edges whose vertex pairs are each other's reverse.
    /// Erase the order from the welding key and they fold into one edge (the red probe this
    /// fence was built against).
    ///
    /// ★★ The rim skip's two sides, on one store: the **cut** circle mints no closed `[v, v]`
    /// edge (before the skip, both fixtures minted one that only the reject discarded), while
    /// the **uncut** far rim still mints exactly one — the skip's negative control, pinned to
    /// the far cap's axis coordinate so a skip that turned into "skip every rim" reddens here.
    ///
    /// ★ The two populations differ on the chord, deliberately: the straddling boss's branch
    /// pair is joined by the plate-top chord (welded with the subdivided middle piece into
    /// **one** line edge used by both coplanar faces — the subdivision cell's promise realized
    /// in the store), while the turned boss's pair sits across the plate corner, joined through
    /// it by split boundary edges — no chord at all. Scoped to the edges the boolean minted
    /// (a snapshot, not a whole-store filter: the input cylinder's own rims are `[v, v]` too).
    #[test]
    fn an_arc_and_its_complement_are_minted_as_two_ordered_edges() {
        let s = 2.0 - 3.0f64.sqrt() / 4.0;
        // `seam_split`: where θ = 0 sits. `None` = a branch vertex lies on the seam generator
        // (the straddling boss — the split's own `SeamIncident` case), so no piece splits;
        // `Some(p)` = the seam vertex S is minted at `p` (centre + ref_dir·r, derived) and the
        // wrap arc is cut there into two pieces.
        for (origin, axis, height, crossings, chord_edges, seam_split) in [
            (
                [4.0, 2.0, 2.0],
                [0.0, 0.0, 1.0],
                1.0,
                [[4.0, 1.5, 2.0], [4.0, 2.5, 2.0]],
                1usize,
                None,
            ),
            (
                [4.0, 0.25, 2.0],
                [1.0, 0.0, 0.0],
                1.0,
                [[4.0, 0.75, 2.0], [4.0, 0.0, s]],
                0usize,
                Some([4.0, 0.25, 1.5]),
            ),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                height,
            );
            m.rebuild_adjacency();
            let minted_from = m.edges.len();
            let vertices_from = m.vertices.len();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            let minted: Vec<_> = m.edges.iter().skip(minted_from).collect();
            let is_cyl = |sh| matches!(m.surface(sh), nacre_geom::Surface::Cylinder(_));
            // A circle-carrier edge is a **mixed** pair (the cap plane and the lateral); the
            // band's seam edge is `[lat, lat]` — both carriers the cylinder — and is not a
            // piece of any circle.
            let on_circle = |e: &nacre_topo::Edge| is_cyl(e.surfaces[0]) != is_cyl(e.surfaces[1]);

            // The cut circle's pieces: their directed vertex pairs chain into **one** cycle —
            // the observable of the `[A, B]`-CCW convention (for two pieces the cycle *is* the
            // mutually-reversed pair the first version of this fence asserted), and of the seam
            // split (three pieces when S stands, still one circle).
            let arcs: Vec<_> = minted
                .iter()
                .filter(|(_, e)| on_circle(e) && e.vertices[0] != e.vertices[1])
                .collect();
            let expect_pieces = 2 + usize::from(seam_split.is_some());
            assert_eq!(
                arcs.len(),
                expect_pieces,
                "the cut circle's pieces: {arcs:?}"
            );
            let mut succ: std::collections::HashMap<_, _> = std::collections::HashMap::new();
            for (_, e) in &arcs {
                let prev = succ.insert(e.vertices[0], e.vertices[1]);
                assert!(prev.is_none(), "two pieces leave one vertex CCW: {arcs:?}");
            }
            let start = arcs[0].1.vertices[0];
            let (mut cur, mut steps) = (start, 0usize);
            loop {
                cur = succ[&cur];
                steps += 1;
                if cur == start || steps > expect_pieces {
                    break;
                }
            }
            assert_eq!(steps, expect_pieces, "the pieces close one circle");
            // Endpoints: exactly two branch vertices on the derived crossings, plus — when the
            // seam splits an arc — one OnSeam vertex at the derived seam point, with its
            // tolerance measured.
            let (mut branch, mut on_seam) = (Vec::new(), Vec::new());
            for &v in succ.keys() {
                match m.vertices.get(v).def {
                    nacre_topo::VertexDef::Branch { .. } => branch.push(v),
                    nacre_topo::VertexDef::OnSeam(_) => on_seam.push(v),
                    ref d => panic!("an arc endpoint is neither branch nor seam: {d:?}"),
                }
            }
            assert_eq!(branch.len(), 2, "one branch pair");
            for &v in &branch {
                let p = m.vertex_point(v);
                assert!(
                    crossings
                        .iter()
                        .any(|c| (0..3).all(|i| (p.as_array()[i] - c[i]).abs() < 1e-9)),
                    "an arc endpoint sits on a derived crossing: {p:?}"
                );
            }
            match seam_split {
                None => assert!(on_seam.is_empty(), "the seam is the branch vertex itself"),
                Some(sp) => {
                    assert_eq!(on_seam.len(), 1, "one seam vertex on the cut circle");
                    let p = m.vertex_point(on_seam[0]);
                    assert!(
                        (0..3).all(|i| (p.as_array()[i] - sp[i]).abs() < 1e-12),
                        "S sits on the derived seam point: {p:?}"
                    );
                    let tol = m.vertex_tol(on_seam[0]).expect("S is discovered");
                    assert!(tol < 1e-12, "measured against what defines it: {tol}");
                }
            }
            // The minted OnSeam census: the uncut far rim's vertex, plus S when it stands —
            // and nothing else (a duplicate S at a seam-incident branch vertex would show here).
            let minted_on_seam = m
                .vertices
                .iter()
                .skip(vertices_from)
                .filter(|(_, v)| matches!(v.def, nacre_topo::VertexDef::OnSeam(_)))
                .count();
            assert_eq!(
                minted_on_seam,
                1 + on_seam.len(),
                "far rim + S, nothing else"
            );

            // The rim skip's two sides: no closed edge on the cut circle, exactly one on the
            // uncut far rim (its axis coordinate is the far cap's, derived from the fixture).
            let closed: Vec<_> = minted
                .iter()
                .filter(|(_, e)| e.vertices[0] == e.vertices[1])
                .collect();
            assert_eq!(closed.len(), 1, "one uncut rim, no cut one: {closed:?}");
            let far = (0..3)
                .map(|i| (origin[i] + axis[i] * height) * axis[i])
                .sum::<f64>();
            let p = m.vertex_point(closed[0].1.vertices[0]).as_array();
            let along = (0..3).map(|i| p[i] * axis[i]).sum::<f64>();
            assert!(
                (along - far).abs() < 1e-12,
                "the closed rim is the far cap's: {along} vs {far}"
            );

            // The chord: welded into one line edge on the straddling boss, absent across the
            // turned boss's corner.
            let (ba, bb) = (branch[0], branch[1]);
            let chords = minted
                .iter()
                .filter(|(_, e)| {
                    !on_circle(e) && (e.vertices == [ba, bb] || e.vertices == [bb, ba])
                })
                .count();
            assert_eq!(chords, chord_edges, "the branch pair's line edges");
        }
    }

    /// **The band's loop is one continuous cycle, and the shell closes over it.**
    ///
    /// ★★★ The band assembles its cut rim from the arc chain; this fence counts the closure
    /// directly on the store beside the production guard: every edge the boolean minted is used
    /// exactly twice across its faces, and the band face's outer loop is one vertex-continuous
    /// cycle of the derived length, with the seam edge traversed once in each sense.
    ///
    /// ★ Three fixtures: the straddling boss (lo rim cut, seam ≡ branch), the turned boss
    /// (lo rim cut, seam splits the wrap arc — six half-edges), and the **hung** boss (the
    /// straddling boss mirrored under the plate — the cut circle is the band's **hi** end, so
    /// the chain is walked reversed; measured to pass the gate before this fence was written).
    /// The hi-cut *and* seam-split combination has no fixture yet — the chain logic is shared,
    /// and its population brings one when it arrives.
    #[test]
    fn the_bands_loop_is_one_continuous_cycle() {
        for (origin, axis, band_len) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0], 5usize),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0], 6usize),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0], 5usize),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let faces_from = m.faces.len();
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let garbage: Vec<_> = m.faces.iter().skip(faces_from).collect();
            assert!(!garbage.is_empty(), "the face loop ran to completion");

            // Closure: every edge of the garbage faces is used exactly twice.
            let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
            for (_, f) in &garbage {
                for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                    for he in &lp.half_edges {
                        *uses.entry(he.edge).or_default() += 1;
                    }
                }
            }
            let odd: Vec<_> = uses.iter().filter(|&(_, &n)| n != 2).collect();
            assert!(
                odd.is_empty(),
                "a closed shell uses every edge twice: {odd:?}"
            );

            // The band face: one, on the cylinder, its outer loop a continuous cycle of the
            // derived length, the seam edge once in each sense.
            let bands: Vec<_> = garbage
                .iter()
                .filter(|(_, f)| matches!(m.surface(f.surface), nacre_geom::Surface::Cylinder(_)))
                .collect();
            assert_eq!(bands.len(), 1, "one band face");
            let lp = &bands[0].1.outer;
            assert_eq!(lp.half_edges.len(), band_len, "the derived loop length");
            let ends = |he: &nacre_topo::HalfEdge| {
                let [a, b] = m.edges.get(he.edge).vertices;
                if he.forward { (a, b) } else { (b, a) }
            };
            for w in 0..lp.half_edges.len() {
                let (_, e0) = ends(&lp.half_edges[w]);
                let (s1, _) = ends(&lp.half_edges[(w + 1) % lp.half_edges.len()]);
                assert_eq!(e0, s1, "the loop chains vertex to vertex at step {w}");
            }
            let mut seen: std::collections::HashMap<_, Vec<bool>> =
                std::collections::HashMap::new();
            for he in &lp.half_edges {
                seen.entry(he.edge).or_default().push(he.forward);
            }
            let twice: Vec<_> = seen.values().filter(|v| v.len() == 2).collect();
            assert_eq!(twice.len(), 1, "exactly one edge is walked twice: the seam");
            assert_ne!(twice[0][0], twice[0][1], "once in each sense");
        }
    }

    /// **The grouping joins across a cut rim — through the door production uses.**
    ///
    /// ★★★ A cut circle bounds no whole disk, so the rim-key rule cannot see it, and the
    /// unordered node rule cannot either: on a 2-node circle the chord and both complementary
    /// arcs fold into one `norm_edge` pair (six users where each piece has two). The `JoinKey`'s
    /// ordered arc pairs — "a line is unordered, a circle is ordered", third appearance — give
    /// every piece exactly its cap and the band, so the result comes back as **one** component.
    ///
    /// ★ Asserted on `name_result_vertices`' own product (the subdivision must run first for the
    /// chord's line key to match), across all three fixtures. red: the band's registration
    /// removed → n == 2 and the held grouping is the old `BranchVertexUnnamed`.
    #[test]
    fn the_grouping_joins_across_a_cut_rim() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let setup = plane_index_setup(&m, plate, boss).unwrap();
            let PlaneSetup {
                planes: faces_tab,
                geom,
                surf_ix,
                inc_a,
                inc_b,
                plane_ix,
                class_owner,
                n_a,
                standard,
                notes,
                cyls,
                ..
            } = &setup;
            let jd = Judge::new(geom, *standard, notes);
            let trace_in = crate::combinatorics::trace_input(
                &m,
                [(plate, inc_a), (boss, inc_b)],
                surf_ix,
                faces_tab.len(),
                &jd,
                plane_ix,
                cyls,
                Default::default(),
            );
            let (plane_faces, curved, _) = crate::arrangement::trace_result_faces_full_for_test(
                &m,
                BoolKind::Fuse,
                plate,
                boss,
                &jd,
                faces_tab,
                plane_ix,
                cyls,
                *n_a,
                class_owner,
                &trace_in,
            )
            .expect("the arc population traces");
            let faces =
                crate::boolean::unify_coplanar_faces(plane_faces, &jd, &setup.cyls).expect("unify");
            let rows = cyl_rows(faces_tab, plane_ix, *n_a).expect("rows");
            let mut faces = faces;
            faces.extend(
                band_faces(
                    BoolKind::Fuse,
                    &faces,
                    &rows,
                    &jd,
                    &curved.disk_labels,
                    &curved.cut_rims,
                    &curved.arc_labels,
                )
                .expect("bands"),
            );
            let seam = crate::arrangement::seam_table(&faces, cyls, &jd)
                .expect("the seam realizes branch nodes");
            let named =
                crate::boolean::name_result_vertices(&jd, &seam, &faces, cyls, &curved.cut_rims)
                    .expect("the naming runs");
            let g = named
                .grouping
                .as_ref()
                .expect("the grouping joins across the cut rim");
            assert_eq!(g.n, 1, "one component");
            assert_eq!(g.positives, vec![0], "one material piece, no cavity");
        }
    }

    /// **A cut-rim boolean builds a complete solid, and the integrals pin it exactly.**
    ///
    /// ★★★ The volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap — every fixture
    /// is a contact), so a wrong segment sign or scale moves an exact number; and because all
    /// three fixtures share that one number, the straddling boss's two z = 2 faces additionally
    /// split the segment **signs** (the plate top loses its half-disk bite, the digon *is* the
    /// bite). The structure is pinned beside it: one solid, no cavities, its outer shell exactly
    /// the boolean's minted faces.
    #[test]
    fn a_cut_rim_boolean_builds_a_complete_solid() {
        for (origin, axis) in [
            ([4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ([4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ([4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            let mut m = Model::new();
            let plate = m.add_cuboid(
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let boss = m.add_cylinder(
                Point3::from_array(origin),
                Vector3::from_array(axis),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let (faces_from, solids_from) = (m.faces.len(), m.solids.len());
            let before = m.live_solids.clone();
            let out = crate::boolean(&mut m, BoolKind::Fuse, plate, boss)
                .expect("the cut-rim boolean builds");
            assert_eq!(out.len(), 1, "one fused solid");
            assert_ne!(m.live_solids, before, "the operands retired");
            assert_eq!(m.live_solids, out, "the result lives");
            let pushed: Vec<_> = m.solids.iter().skip(solids_from).collect();
            let [(sh, solid)] = pushed[..] else {
                panic!("one result solid, got {}", pushed.len());
            };
            assert_eq!(sh, out[0], "the pushed solid is the returned one");
            assert!(solid.cavities.is_empty(), "one material piece, no cavity");
            let shell_faces: std::collections::HashSet<_> =
                m.shells.get(solid.outer).faces.iter().copied().collect();
            let minted: std::collections::HashSet<_> =
                m.faces.iter().skip(faces_from).map(|(h, _)| h).collect();
            assert_eq!(
                shell_faces, minted,
                "the outer shell is exactly the boolean's minted faces"
            );
            // ★★ **The arc integrals answer exactly** —
            // the volume is derived (plate 4·4·2 = 32, boss π·0.25·1, zero overlap: every
            // fixture is a contact), so a wrong segment sign or scale moves an exact number.
            let props = nacre_props::mass_props(&m, sh).expect("the mixed loops integrate");
            let expect = 32.0 + std::f64::consts::PI / 4.0;
            assert!(
                (props.volume - expect).abs() < 1e-12,
                "volume {} vs derived {expect}",
                props.volume
            );
            // ★ All three fixtures share that one number, so a global scale error fools them
            // together — the straddling boss's two z = 2 faces split the segment **signs**: the
            // plate top loses its half-disk bite (16 − π/8), the overhang digon *is* the other
            // half-disk (π/8).
            if origin == [4.0, 2.0, 2.0] {
                let (mut top, mut digon) = (None, None);
                for (h, f) in m.faces.iter().skip(faces_from) {
                    let nacre_geom::Surface::Plane(p) = m.surface(f.surface) else {
                        continue;
                    };
                    if (p.normal().as_array()[2] - 1.0).abs() > 1e-9 {
                        continue;
                    }
                    match f.outer.half_edges.len() {
                        2 => digon = Some(h),
                        6 => top = Some(h),
                        1 => {} // the boss's closed top cap
                        n => panic!("an unexpected z-normal face with {n} half-edges"),
                    }
                }
                let area = |h| nacre_props::face_props(&m, h).expect("planar").area;
                let bite = std::f64::consts::PI / 8.0;
                assert!(
                    (area(top.expect("plate top")) - (16.0 - bite)).abs() < 1e-12,
                    "the plate top loses its bite"
                );
                assert!(
                    (area(digon.expect("overhang digon")) - bite).abs() < 1e-12,
                    "the digon is the bite"
                );
            }
        }
    }

    #[test]
    #[ignore = "scratch OBJ dump for eyeballing — run on demand"]
    fn dump_straddling_boss_obj() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = m.add_cylinder(
            Point3::from_array([4.0, 2.0, 2.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Fuse, plate, boss).expect("builds");
        m.rebuild_adjacency();
        let obj = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("tessellates")
            .to_obj();
        // Written only on request — the `--ignored` sweep runs this test too, and a test that
        // writes outside the workspace on every sweep is a side effect nobody asked for.
        match std::env::var("OBJ_OUT") {
            Ok(path) => {
                std::fs::write(&path, obj).expect("write");
                println!("wrote {path}");
            }
            Err(_) => println!("set OBJ_OUT=<path> to write the OBJ ({} bytes)", obj.len()),
        }
    }

    /// The chained bore-then-boss body: `cut` a through-bore at `(12,12)`, then `fuse` a boss
    /// whose rim the plate's edge cuts. The chain is what the milestone fence could not ask:
    /// nesting must place the bore's rim circle inside a **bitten** top ring — two branch
    /// corners and an arc step — where the chart road's parity had no rational corners to read
    /// (`WitnessNotRational`, chaining wall 3) and the mixed parity answers in ℚ(√c).
    fn a_bored_plate_with_a_boss(boss_base: [f64; 3]) -> f64 {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([12.0, 12.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            4.0,
            20.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        let boss = m.add_cylinder(
            Point3::from_array(boss_base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            10.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, bored, boss).expect("the boss fuses");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, out, "the operands retired, the result lives");
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{issues:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default())
            .expect("the arcs tessellate");
        let mut uses: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (a, b) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "the mesh is watertight");
        nacre_props::mass_props(&m, out[0]).expect("props").volume
    }

    /// ★ **A bored plate takes a straddling boss — the first two-boolean chain through the arc
    /// population is green.** Volume `40·40·20 − π·4²·20 + π·5²·10 = 32000 − 70π` (the boss
    /// stands wholly above the plate top, so fuse adds its full cylinder).
    #[test]
    fn a_bored_plate_takes_a_straddling_boss() {
        let v = a_bored_plate_with_a_boss([40.0, 20.0, 20.0]);
        let want = 32000.0 - 70.0 * std::f64::consts::PI;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The mirror: the boss hangs under the plate's bottom edge — same chain, same volume, the
    /// cut circle at the band's hi end instead of lo.
    #[test]
    fn a_bored_plate_hangs_a_straddling_boss() {
        let v = a_bored_plate_with_a_boss([40.0, 20.0, -10.0]);
        let want = 32000.0 - 70.0 * std::f64::consts::PI;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **The chained contact-cut names its refusal.** Cut the same chain instead of fusing:
    /// the boss only *touches* the bored plate's top, so the cut removes nothing — but the
    /// assembly still keeps the top ring's branch vertices, whose definitions name the boss's
    /// cylinder that the result has no face on. That used to be a debug_assert (dev panic,
    /// release shipped the mis-named solid silently); the floor is now an honest reject with
    /// the offending corner as its witness, and the operands stay live.
    #[test]
    fn a_chained_contact_cut_names_its_refusal() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 20.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([12.0, 12.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            4.0,
            20.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        let boss = m.add_cylinder(
            Point3::from_array([40.0, 20.0, 20.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            5.0,
            10.0,
        );
        m.rebuild_adjacency();
        let live_before = m.live_solids.clone();
        let err = crate::boolean(&mut m, BoolKind::Cut, bored, boss).expect_err("contact-cut");
        let BoolError::Rejected { reason, at } = err else {
            panic!("{err:?}");
        };
        assert_eq!(reason, RejectReason::VertexNamesAbsentSurface);
        // The witness is one of the two branch corners where the boss's rim crosses the
        // plate's edge: (40, 20±√(25−0), 20) → y = 15 or 25 (approximate-only assert).
        let Some(crate::RejectWhere::Point(p)) = at else {
            panic!("{at:?}");
        };
        let p = p.as_array();
        assert!(
            (p[0] - 40.0).abs() < 1e-9
                && (p[2] - 20.0).abs() < 1e-9
                && ((p[1] - 15.0).abs() < 1e-9 || (p[1] - 25.0).abs() < 1e-9),
            "{p:?}"
        );
        m.rebuild_adjacency();
        assert_eq!(m.live_solids, live_before, "a reject restores the live set");
    }

    /// **Two cylinders with coplanar caps fuse apart.** They stand `5` apart with their caps in
    /// the same two planes, which is what once made them look like a seating problem; what was
    /// actually hard was classifying **two curved bodies**, and both are now probed from a cap
    /// disk's centre — neither has a vertex to probe from.
    ///
    /// ★ Both volumes are the same number derived the same way (`π·0.5²·2`), so the assertion
    /// cannot pass by matching one body against the other.
    #[test]
    fn two_cylinders_with_coplanar_caps_fuse_apart() {
        let mut m = Model::new();
        let a = m.add_cylinder(
            Point3::from_array([0.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        let b = m.add_cylinder(
            Point3::from_array([5.0, 0.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, a, b).expect("two curved bodies");
        assert_eq!(out.len(), 2, "they do not touch");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let want = std::f64::consts::PI * 0.25 * 2.0;
        for &s in &out {
            let v = nacre_props::mass_props(&m, s).expect("props").volume;
            assert!((v - want).abs() < 1e-9, "an untouched cylinder: {v}");
            let (v_n, e_n, f_n, l_n) = euler_counts(&m, s);
            assert_eq!(
                v_n - e_n + f_n - l_n,
                2,
                "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
            );
        }
    }

    /// `χ = V − E + F − L` over **one solid's** shells — the counts the genus relation
    /// `χ = 2(S − G)` is read from.
    ///
    /// ★ Per solid, not per model: `Model::reachable` spans every live solid, which coincides with
    /// this only while a fixture makes exactly one. The fixtures that make two need the sum split,
    /// and one spelling for all of them is what keeps the relation from being restated.
    fn euler_counts(m: &Model, s: Handle<Solid>) -> (i64, i64, i64, i64) {
        let faces: Vec<_> = crate::planes::solid_shell_handles(m, s)
            .into_iter()
            .flat_map(|sh| m.shells.get(sh).faces.clone())
            .collect();
        let mut verts = std::collections::HashSet::new();
        let mut edges = std::collections::HashSet::new();
        let mut loops = 0i64;
        for &fh in &faces {
            let f = m.faces.get(fh);
            loops += f.inner.len() as i64;
            for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
                for he in &lp.half_edges {
                    edges.insert(he.edge);
                    verts.extend(m.edges.get(he.edge).vertices);
                }
            }
        }
        (
            verts.len() as i64,
            edges.len() as i64,
            faces.len() as i64,
            loops,
        )
    }

    // ---- Wall faces: the gate asks the boundary, not the infinite plane (M6-2b preparation) ----
    //
    // ★★ The uniform-slab theorem's premise is about the other operand's **boundary**. The gate
    // used to test the wall's infinite *plane*, which is a cheaper sufficient condition — and it
    // refused a whole family the engine serves: a body standing well clear of a bore whose wall
    // plane, extended, happens to pass through it. These lock what that opened and what it did not.

    /// A plate with a bore, fused to a boss standing far away in `x` — whose `y = 12` wall **plane**
    /// stands only 2 from the bore's axis (`r = 3`). The boss's face is 20 away; nothing meets.
    fn plate_bore_and_boss(hole_y: f64, boss_z0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, hole_y, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let boss = m.add_cuboid(
            Point3::from_array([28.0, 4.0, boss_z0]),
            Point3::from_array([36.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        (m, holed, boss)
    }

    #[test]
    fn a_boss_whose_wall_plane_crosses_a_distant_bore_still_fuses() {
        let (mut m, holed, boss) = plate_bore_and_boss(10.0, 5.0);
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("the boss fuses");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **The d = 0 member of the family: the wall plane runs exactly through the bore's
    /// axis** (M6-2 rulings ladder). The gate passes it the same way (the boss's face clears the
    /// footprint along the span), and the rulings road stays **silent** — its trigger is the
    /// gate's carried `crossings` record, which is empty for every gate-passed pair. This
    /// geometry had no corpus row (its siblings are all d = 2), and it is the population an
    /// unconditionally-firing contribution arm broke — measured, 14 arc tests red — so it locks
    /// "an empty record changes nothing".
    #[test]
    fn a_boss_whose_wall_plane_holds_the_bores_axis_still_fuses() {
        let (mut m, holed, boss) = plate_bore_and_boss(12.0, 5.0);
        let out = crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("the boss fuses");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The d = 0 member with **caps in play**: a cylinder tool standing clear of the plate, its
    /// axis exactly on the plate's `x = 40` wall plane and its caps strictly inside the plate's
    /// height (no coplanar contact anywhere). Two disjoint bodies is the valid fuse answer, and
    /// the chord arm — like the ruling arm — stays silent behind the empty record.
    #[test]
    fn a_capped_tool_on_the_plates_wall_plane_fuses_as_two_bodies() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let tool = m.add_cylinder(
            Point3::from_array([40.0, 30.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            3.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, plate, tool).expect("a distant fuse");
        assert_eq!(out.len(), 2, "two disjoint bodies");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [std::f64::consts::PI * 9.0 * 3.0, 4000.0];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "{vols:?} vs {want:?}"
        );
    }

    /// The `Cut` twin — the same wall planes reach the gate whichever way the operation runs.
    /// The tool is sunk into the plate (`z` from 3) so it removes material: a body merely *resting*
    /// on the top face is a coplanar contact and meets a different guard entirely, which would make
    /// this fixture measure that instead of the wall rule.
    #[test]
    fn the_same_boss_cuts_the_bored_plate() {
        let (mut m, holed, boss) = plate_bore_and_boss(10.0, 3.0);
        let out = crate::boolean(&mut m, BoolKind::Cut, holed, boss).expect("the boss cuts");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 9.0 * 5.0 - 8.0 * 8.0 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **Two bores, one wall plane crossing both.** The rule is per (cylinder, class), so a second
    /// cylinder is a second set of questions about the same face — and the face answers for each.
    #[test]
    fn one_wall_plane_may_cross_two_bores() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let mut holed = plate;
        for x in [8.0, 20.0] {
            let hole = m.add_cylinder(
                Point3::from_array([x, 10.0, -1.0]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                3.0,
                7.0,
            );
            m.rebuild_adjacency();
            holed = crate::boolean(&mut m, BoolKind::Cut, holed, hole).expect("bore")[0];
        }
        let boss = m.add_cuboid(
            Point3::from_array([28.0, 4.0, 5.0]),
            Point3::from_array([36.0, 12.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out =
            crate::boolean(&mut m, BoolKind::Fuse, holed, boss).expect("two bores and a boss");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - 2.0 * std::f64::consts::PI * 9.0 * 5.0 + 8.0 * 8.0 * 3.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **The fence: a wall face that really does cross the bore.** The plate's own `y = 20` wall
    /// would clear, so the tool here is a slab whose face runs right across the hole — the wall
    /// rule's true population, and M6-2b's.
    #[test]
    fn a_wall_face_that_really_crosses_the_bore_is_refused() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // A slab covering x ∈ [0, 40]: its y = 12 face passes straight through the bore.
        let slab = m.add_cuboid(
            Point3::from_array([0.0, 12.0, 0.0]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        m.rebuild_adjacency();
        let before = m.live_solids.clone();
        let err = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect_err("a real crossing");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
        assert_eq!(m.live_solids, before, "a refused boolean retires nothing");
    }

    /// **The fence: exact tangency.** The slab's face stands exactly `r` from the axis, so it
    /// touches the lateral surface along one ruling — the boundary case, and `Inside` includes it.
    #[test]
    fn a_wall_face_tangent_to_the_bore_is_refused() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let holed = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        let slab = m.add_cuboid(
            Point3::from_array([0.0, 13.0, 0.0]), // y = 13 is exactly r = 3 from y = 10
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, holed, slab).expect_err("tangency");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    /// ★★ **The fence that measures the straight-edge barrier** — the only fixture that does.
    ///
    /// The plate carries a **crosswise** bore, so the plane `x = 30` holds that bore's circular cap
    /// face: an outer loop of one arc and one seam vertex. That plane is also parallel to the
    /// vertical drill's axis and passes within `r` of it, so the face test is reached — and a rule
    /// reading vertices alone would find "every vertex on one side" true of a **single point** and
    /// wave a disk straight through the strip. Asking the edge's curve kind is what stops it.
    #[test]
    fn a_disk_face_on_a_wall_plane_is_undecided_not_waved_through() {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 10.0]),
        );
        // A crosswise bore along +X, ending inside the plate at x = 30: its cap lies on x = 30.
        let cross = m.add_cylinder(
            Point3::from_array([-1.0, 10.0, 5.0]),
            Vector3::from_array([1.0, 0.0, 0.0]),
            2.0,
            31.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, cross).expect("crosswise bore")[0];
        // A vertical drill whose axis stands 1 from the plane x = 30 — inside its radius 3 — and
        // ★ **6 from the crosswise bore's axis**, clear of the radius sum 5, so the pair rule is
        // not what answers here. (At y = 10 the two axes actually meet and this fixture would be
        // measuring `CylinderPairContact` instead — the adjacent proposition.)
        let drill = m.add_cylinder(
            Point3::from_array([31.0, 4.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            12.0,
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, bored, drill)
            .expect_err("a disk face on the wall plane cannot be judged by its vertices");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::CylinderGateUndecided,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    // ---- One surface, several lateral faces: the band belongs to the face ----

    /// The shape at the heart of it: a bored plate whose bore's **middle** is cut away, leaving two
    /// disjoint bands on one lateral surface. Returns the model and that solid.
    fn plate_with_a_split_bore(cut_x0: f64, cut_x1: f64) -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
        let hole = m.add_cylinder(
            Point3::from_array([5.0, 5.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            2.0,
            12.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
        // Walls 3 from the axis (r = 2), so the wall rule is not what this measures.
        let mid = m.add_cuboid(
            Point3::from_array([cut_x0, cut_x0, 4.0]),
            Point3::from_array([cut_x1, cut_x1, 6.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, bored, mid).expect("the middle cut");
        (m, out[0])
    }

    /// The lid tool the locks below cut with.
    fn add_lid(m: &mut Model) -> Handle<Solid> {
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        lid
    }

    /// How many faces of `s` lie on cylinder surfaces, grouped by surface.
    fn lateral_face_counts(m: &Model, s: Handle<Solid>) -> Vec<usize> {
        let mut counts: std::collections::HashMap<Handle<nacre_geom::Surface>, usize> =
            Default::default();
        let solid = m.solids.get(s);
        for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
            for &fh in &m.shells.get(sh).faces {
                let surf = m.faces.get(fh).surface;
                if matches!(
                    m.surface_truth(surf),
                    nacre_topo::SurfaceTruth::Cylinder { .. }
                ) {
                    *counts.entry(surf).or_default() += 1;
                }
            }
        }
        let mut v: Vec<usize> = counts.into_values().collect();
        v.sort_unstable();
        v
    }

    /// **Making it.** The split bore is an ordinary, correct solid — and nothing measured that
    /// until this defect was found, which is why it is its own test: if building it ever breaks,
    /// the test below (which *uses* it) must not be the one that goes red.
    #[test]
    fn a_bore_cut_across_the_middle_leaves_two_bands_on_one_surface() {
        let (m, s) = plate_with_a_split_bore(2.0, 8.0);
        assert_eq!(
            lateral_face_counts(&m, s),
            vec![2],
            "one lateral surface, two faces"
        );
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, s).expect("props").volume;
        // 1000 − bore(π·4·10) − the middle cube outside the bore(6·6·2 − π·4·2)
        let want = 1000.0 - std::f64::consts::PI * 40.0 - (72.0 - std::f64::consts::PI * 8.0);
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **Using it.** The two-banded solid as an operand. `cyl_rows` used to let the first
    /// lateral face speak for the class, so the second band's span was clipped away and never
    /// emitted — the result shell came back open (`OpenResultShell`), naming a symptom of our own
    /// omission rather than anything about the input.
    #[test]
    fn a_solid_with_two_bands_on_one_surface_can_be_cut_again() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        // ★ The property under repair, counted directly — a volume alone does not say "two bands
        // came out", and this defect is exactly one band going missing.
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![2],
            "both bands survive the cut"
        );
        let p = nacre_props::mass_props(&m, out[0]).expect("props");
        let want = 1000.0
            - std::f64::consts::PI * 40.0
            - (72.0 - std::f64::consts::PI * 8.0)
            - (200.0 - std::f64::consts::PI * 8.0);
        assert!((p.volume - want).abs() < 1e-9, "{} vs {want}", p.volume);
        // ★★ **Area is the sharper oracle for this defect.** A missing band does not change the
        // volume — the boolean simply refuses — but it takes its `2πr·h` out of the surface. The
        // two bands here are `z ∈ [0,4]` and `[6,8]`, so `16π + 8π` of lateral area must be in
        // this number: 640 + 8π once every planar face is counted.
        let want_area = 640.0 + std::f64::consts::PI * 8.0;
        assert!(
            (p.area - want_area).abs() < 1e-9,
            "{} vs {want_area}",
            p.area
        );
    }

    /// ★ **The negative control.** Move the middle cut into a corner, away from the bore: the wall
    /// stays **one** face and the same three steps already worked before the fix. Without this, a
    /// change that merely made *any* third boolean succeed would look like a repair.
    #[test]
    fn a_middle_cut_that_misses_the_bore_leaves_one_band() {
        let (mut m, s) = plate_with_a_split_bore(0.0, 2.0);
        assert_eq!(lateral_face_counts(&m, s), vec![1], "the wall is untouched");
        let lid = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 8.0]),
            Point3::from_array([10.0, 10.0, 12.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, lid).expect("the lid cut");
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want =
            1000.0 - std::f64::consts::PI * 40.0 - 8.0 - (200.0 - std::f64::consts::PI * 8.0);
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// **The mechanism, read directly.** Two rows, one per lateral face, with the spans the two
    /// bands actually occupy — and the `t` axis here starts at the drill's origin `z = −1`, so the
    /// bands `z ∈ [0,4]` and `[6,10]` are `t ∈ [1,5]` and `[7,11]`.
    #[test]
    fn cyl_rows_gives_one_row_per_lateral_face() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = add_lid(&mut m);
        let setup = plane_index_setup(&m, s, lid).unwrap();
        let rows = cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).expect("rows");
        let spans: Vec<[f64; 2]> = rows
            .iter()
            .map(|r| [r.span[0].to_f64(), r.span[1].to_f64()])
            .collect();
        assert_eq!(
            spans,
            vec![[1.0, 5.0], [7.0, 11.0]],
            "one row per face, in t order"
        );
        assert!(rows.iter().all(|r| r.class == 0), "both on the one surface");
    }

    /// ★★ **The band pass makes nothing in the gap.** This is where the right fix parts company
    /// with the plausible one: merging the two faces' spans into a single `min..max` would put a
    /// band across `z ∈ [4,6]`, where the solid has no lateral face at all.
    #[test]
    fn no_band_is_invented_where_the_solid_has_no_face() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let lid = add_lid(&mut m);
        let (out, ts) = bands(&m, s, lid, BoolKind::Cut);
        let mut spans: Vec<(f64, f64)> = out.iter().map(|lf| ends(lf, &ts)).collect();
        spans.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for (lo, hi) in &spans {
            // Disjoint from the open gap (4, 6): a band either ends at or below 4, or starts at or
            // above 6. ★ Measured: with the spans merged into one `min..max` this fires with
            // "a band 4..6 crosses the gap", which is the whole reason the clause is here.
            assert!(
                *hi <= 4.0 || *lo >= 6.0,
                "a band {lo}..{hi} crosses the gap z 4..6 where there is no face"
            );
        }
        assert!(
            spans.iter().any(|(lo, _)| (*lo - 0.0).abs() < 1e-9),
            "the lower band is emitted: {spans:?}"
        );
        assert!(
            spans.iter().any(|(lo, _)| (*lo - 6.0).abs() < 1e-9),
            "the upper band is emitted too: {spans:?}"
        );
    }

    // ---- The footprint's other axis: along the cylinder (M6-2b preparation) ----
    //
    // ★★ A wall parallel to the axis meets the cylinder in a **rectangle** of the wall's own
    // plane — the strip across, a lateral face's span along — so "does this face miss it" is one
    // question with two separating axes. The gate used to read only the first, and refused every
    // wall that stood clear of the cylinder along its length.
    //
    // ★ With several lateral faces there are several rectangles: the strip is shared, the spans
    // are not. That is why the gap between two bands is passable at all, and it is not a special
    // case bolted on — `bands_of` clips every cut to its own row's span, so no band is ever built
    // in a gap and there is no premise there for a wall to break.

    /// **The case this opened.** The two-banded bore of `plate_with_a_split_bore` has bands at
    /// `z ∈ [0,4]` and `[6,10]` — `t ∈ [1,5]` and `[7,11]`. A tool whose `y = 6` wall stands only 1
    /// from the axis (`r = 2`) crosses the strip, and its face runs the plate's whole width, so the
    /// axis across cannot clear it. Along the axis it is `t ∈ [5.5, 6.5]` — inside the gap, missing
    /// both rectangles.
    #[test]
    fn a_wall_face_in_the_gap_between_two_bands_clears() {
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let tool = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 4.5]),
            Point3::from_array([10.0, 6.0, 5.5]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the wall clears in t");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![2],
            "the tool passes between the bands and leaves both"
        );
        // 1000 − bore(π·4·10) − the middle cut(6·6·2 − π·4·2) − this tool, which meets the solid
        // over 10×6 minus the 6×4 the middle cut already took, one deep.
        let want = 1000.0
            - std::f64::consts::PI * 40.0
            - (72.0 - std::f64::consts::PI * 8.0)
            - (60.0 - 24.0);
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **The two axes are one judgement, measured as one.** The same tool three ways against
    /// the same two-banded bore, with only the numbers moved: clearing *either* axis is enough,
    /// and clearing neither is the refusal. A rule that had merely gained a second, independent
    /// test would pass (a) and (b) too — what this pins is (c), that the two are OR-ed rather
    /// than each able to wave a face through on its own terms.
    #[test]
    fn either_axis_clears_the_footprint_and_neither_does_not() {
        // (a) Across only: the `y = 6` face sits at `x ∈ [0, 2.5]`, clear of the strip
        //     `x ∈ [3.27, 6.73]`, while its plane still crosses the bore. Along the axis it runs
        //     `t ∈ [4,8]`, straddling both bands. ★ The `x = 2.5` wall must stand clear of the axis
        //     by more than `r`, or *it* becomes the face under test — at `x = 3` it is exactly
        //     tangent and this arm measures that instead (measured: it refuses).
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let across = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 3.0]),
            Point3::from_array([2.5, 6.0, 7.0]),
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, s, across).expect("clear across the strip");

        // (b) Along only: the face spans the full width, so only the gap saves it.
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let along = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 4.5]),
            Point3::from_array([10.0, 6.0, 5.5]),
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, s, along).expect("clear along the axis");

        // (c) Neither: full width *and* straddling both bands.
        let (mut m, s) = plate_with_a_split_bore(2.0, 8.0);
        let neither = m.add_cuboid(
            Point3::from_array([0.0, 0.0, 3.0]),
            Point3::from_array([10.0, 6.0, 7.0]),
        );
        m.rebuild_adjacency();
        let err = crate::boolean(&mut m, BoolKind::Cut, s, neither).expect_err("crosses both");
        assert!(
            matches!(
                err,
                BoolError::Rejected {
                    reason: RejectReason::WallMeetsLateral,
                    ..
                }
            ),
            "{err:?}"
        );
    }

    // ---- A blind bore is usable, not just buildable ----
    //
    // ★★ A cylinder's lateral face touching a ⊥ class at its **rim** contributes a `Graze`, the
    // same as any planar wall meeting a class along an edge. It used to contribute nothing, and
    // `edge_mask`'s `Graze > Seated` precedence exists exactly for the corner where the two
    // disagree — a **blind bore's ceiling**. Missing it, the cap's seated rule flipped the wrong
    // label bit, the band's two ends contradicted each other, and the *next* boolean on that solid
    // came back `CylinderGateUndecided`. A through bore has no ceiling and was always fine.

    /// A plate bored to `hole_h` deep, ready to be operated on again.
    fn plate_with_a_bore(hole_h: f64) -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 20.0, 5.0]),
        );
        let hole = m.add_cylinder(
            Point3::from_array([8.0, 10.0, 0.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            hole_h,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore");
        (m, out[0])
    }

    /// A boss well clear of the bore (walls 4+ from the axis, r = 3), fused onto the plate.
    fn fuse_a_clear_boss(m: &mut Model, s: Handle<Solid>, z0: f64) -> Vec<Handle<Solid>> {
        let boss = m.add_cuboid(
            Point3::from_array([4.0, 4.0, z0]),
            Point3::from_array([12.0, 6.0, 8.0]),
        );
        m.rebuild_adjacency();
        crate::boolean(m, BoolKind::Fuse, s, boss).expect("the boss")
    }

    #[test]
    fn a_blind_bore_can_be_fused_onto_afterwards() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let out = fuse_a_clear_boss(&mut m, s, 5.0);
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![1],
            "the bore's wall survives"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 + 48.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// The `Cut` twin. ★ The tool is sunk into the plate (`z` from 3): a body merely *resting* on
    /// the top face is a coplanar contact that meets `StraightAngle`, a different refusal
    /// entirely, and this fixture would then be measuring that instead of the bore.
    #[test]
    fn a_blind_bore_can_be_cut_afterwards() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let tool = m.add_cuboid(
            Point3::from_array([4.0, 4.0, 3.0]),
            Point3::from_array([12.0, 6.0, 8.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, tool).expect("the cut");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 - 8.0 * 2.0 * 2.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★★ **The script a person actually writes**: a blind hole, then a second hole somewhere
    /// else. It goes through the cylinder-pair rule and two lateral surfaces at once — a road the
    /// boss fixtures above do not take.
    #[test]
    fn a_blind_bore_can_be_drilled_beside() {
        let (mut m, s) = plate_with_a_bore(3.0);
        let second = m.add_cylinder(
            Point3::from_array([20.0, 10.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            3.0,
            7.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, s, second).expect("the second hole");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            lateral_face_counts(&m, out[0]),
            vec![1, 1],
            "two bores, one wall each"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 4000.0 - std::f64::consts::PI * 27.0 - std::f64::consts::PI * 45.0;
        assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    }

    /// ★ **The control**: a through bore and an overshooting one already worked, and must keep
    /// working. Without this, a change that merely made *any* second boolean succeed would look
    /// like the repair.
    #[test]
    fn a_through_bore_is_unaffected() {
        for h in [5.0, 7.0] {
            let (mut m, s) = plate_with_a_bore(h);
            let out = fuse_a_clear_boss(&mut m, s, 5.0);
            assert_eq!(out.len(), 1);
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{:?}",
                nacre_validate::validate(&m)
            );
            let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
            let want = 4000.0 - std::f64::consts::PI * 45.0 + 48.0;
            assert!((v - want).abs() < 1e-9, "h={h}: {v} vs {want}");
        }
    }

    /// **An enclosed cylindrical void is a cavity, not a second body.** A `[0,4]³` box `Cut` by a
    /// radius-`0.5`, height-`2` cylinder wholly inside it: the boolean's two components are the
    /// box's shell and the void's, and the void must come back **depth 1** — odd, so a hole.
    ///
    /// ★★★ **This is the population the coordinate probe exists for, and the only one that makes
    /// it say `true`.** The void's boundary is two disks and a band, so it carries no vertex at
    /// all and `nodes_of` came back empty — which is what `curved_component_depth`, now retired,
    /// used to refuse.
    /// Two *disjoint* bodies exercise the same road, but their answer is "outside" either way, so
    /// a road that always said `false` would pass them; here it would make the void a **second
    /// solid** and the box's volume whole.
    #[test]
    fn an_enclosed_cylindrical_void_is_a_cavity() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let void = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
        assert_eq!(out.len(), 1, "one body, hollow — not two");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            m.solids.get(out[0]).cavities.len(),
            1,
            "the void is a cavity"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (v - want).abs() < 1e-9,
            "the box less the void: {v} vs {want}"
        );
        // Two shells, genus 0 each: `χ = V − E + F − L = 2(S − G) = 4`.
        let (v_n, e_n, f_n, l_n) = euler_counts(&m, out[0]);
        assert_eq!(
            v_n - e_n + f_n - l_n,
            4,
            "two genus-0 shells: V{v_n} E{e_n} F{f_n} L{l_n}"
        );
    }

    /// **A bore and a sealed void in one body.** `[0,8]×[0,4]×[0,4]` with a through bore at
    /// `(2,2)` and, `4` away, a wholly enclosed cylinder at `(6,2)` — `128 − π − π/2`.
    ///
    /// ★ It is the only fixture where the **coordinate** probe meets a face with a hole: the
    /// void has no vertex, so it is probed from a cap centre, and the body it asks about has
    /// annular caps. Two cylinders in one boolean also need `cylinders_clear`, which the `4`
    /// between the axes supplies.
    #[test]
    fn a_bore_and_a_sealed_void_live_in_one_body() {
        let mut m = Model::new();
        let a = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([8.0, 4.0, 4.0]),
        );
        let bore = m.add_cylinder(
            Point3::from_array([2.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            6.0,
        );
        m.rebuild_adjacency();
        let bored = crate::boolean(&mut m, BoolKind::Cut, a, bore).expect("the bore");
        let void = m.add_cylinder(
            Point3::from_array([6.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Cut, bored[0], void).expect("the sealed void");
        assert_eq!(out.len(), 1, "one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        assert_eq!(
            m.solids.get(out[0]).cavities.len(),
            1,
            "the void is a cavity"
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 128.0 - std::f64::consts::PI * 0.25 * 4.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!((v - want).abs() < 1e-9, "the box less both: {v} vs {want}");
    }

    /// **A cavity finds its owner when there is more than one material.** The void fixture above
    /// takes a shortcut the code makes explicit — one material owns every cavity, no search — so
    /// the search itself is only reached when two materials stand apart. Here the hollow box is
    /// fused with a distant second box, and the void must still come back as the **hollow box's**
    /// cavity rather than a third body.
    ///
    /// ★ It is the only fixture that runs `first_deciding` over a **coordinate** probe in the
    /// cavity-owner search: the void has no vertex, and every other multi-material case in the
    /// suite has cavities with corners.
    #[test]
    fn a_cavity_with_no_vertex_still_finds_its_owner() {
        let mut m = Model::new();
        let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([4.0; 3]));
        let void = m.add_cylinder(
            Point3::from_array([2.0, 2.0, 1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            2.0,
        );
        m.rebuild_adjacency();
        let hollow = crate::boolean(&mut m, BoolKind::Cut, a, void).expect("a void inside a box");
        let b = m.add_cuboid(
            Point3::from_array([10.0, 10.0, 10.0]),
            Point3::from_array([14.0, 14.0, 14.0]),
        );
        m.rebuild_adjacency();
        let out = crate::boolean(&mut m, BoolKind::Fuse, hollow[0], b).expect("two bodies apart");
        assert_eq!(out.len(), 2, "the void is a cavity, not a third body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let mut seen: Vec<(f64, usize)> = out
            .iter()
            .map(|&s| {
                (
                    nacre_props::mass_props(&m, s).expect("props").volume,
                    m.solids.get(s).cavities.len(),
                )
            })
            .collect();
        seen.sort_by(|x, y| x.0.partial_cmp(&y.0).expect("finite"));
        let hollowed = 64.0 - std::f64::consts::PI * 0.25 * 2.0;
        assert!(
            (seen[0].0 - hollowed).abs() < 1e-9 && seen[0].1 == 1,
            "the hollow box keeps its one cavity: {seen:?}"
        );
        assert!(
            (seen[1].0 - 64.0).abs() < 1e-9 && seen[1].1 == 0,
            "the second box is solid: {seen:?}"
        );
    }

    /// **A drill that misses.** `Cut` returns the box untouched and `Fuse` returns two bodies —
    /// the population where a cylinder is present but no circle is ever emitted, so the whole
    /// curved path has to stay out of the way.
    #[test]
    fn a_cylinder_that_misses_changes_nothing_and_fuses_apart() {
        let (mut m, a, b) = {
            let mut m = Model::new();
            let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
            let b = m.add_cylinder(
                Point3::from_array([10.0, 10.0, 0.5]),
                Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            (m, a, b)
        };
        let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that misses");
        assert_eq!(cut.len(), 1);
        let v = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
        assert!((v - 8.0).abs() < 1e-12, "the box is untouched: {v}");

        // ★★ **The fuse this test's name promises — and the case that used to abort the kernel.**
        // Two disjoint bodies means `n = 2`, so each is classified by a ray probe. That probe
        // first panicked (it asked a band face for its plane class), then refused by name
        // (`curved_component_depth`, since retired), and now answers: the ray counts a
        // cylinder's crossings, and
        // the bare cylinder — whose boundary carries **no vertex at all** — is probed from a cap
        // disk's centre instead of from a corner.
        //
        // ★ The values are derived from the inputs, not read back from the result: a `[0,2]³` box
        // is `8`, and a radius-`0.5`, height-`1` cylinder is `π/4`. Two separate bodies, each a
        // sphere topologically (`χ = 2`), and `validate` clean.
        let mut m2 = Model::new();
        let a2 = m2.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
        let b2 = m2.add_cylinder(
            Point3::from_array([10.0, 10.0, 0.5]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            1.0,
        );
        m2.rebuild_adjacency();
        let out = crate::boolean(&mut m2, BoolKind::Fuse, a2, b2).expect("two bodies apart");
        assert_eq!(out.len(), 2, "a fuse of two disjoint bodies is two bodies");
        assert!(
            nacre_validate::validate(&m2).is_empty(),
            "{:?}",
            nacre_validate::validate(&m2)
        );
        let mut vols: Vec<f64> = out
            .iter()
            .map(|&s| nacre_props::mass_props(&m2, s).expect("props").volume)
            .collect();
        vols.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        let want = [std::f64::consts::PI * 0.25, 8.0];
        assert!(
            (0..2).all(|i| (vols[i] - want[i]).abs() < 1e-9),
            "the cylinder and the box, each untouched: {vols:?}"
        );
        for &s in &out {
            let (v_n, e_n, f_n, l_n) = euler_counts(&m2, s);
            assert_eq!(
                v_n - e_n + f_n - l_n,
                2,
                "genus 0: V{v_n} E{e_n} F{f_n} L{l_n}"
            );
        }
    }

    /// **A blind hole, end to end.** The bore stops inside the box, so its own cap closes the
    /// bottom — a disk face that only the seated-circle path can emit.
    #[test]
    fn a_blind_hole_is_built_and_measures_what_it_should() {
        let (mut m, a, b) = box_and_drill(-1.0, 2.0); // z ∈ [−1, 1]
        let out = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("the blind bore cuts");
        assert_eq!(out.len(), 1);
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "{:?}",
            nacre_validate::validate(&m)
        );
        let v = nacre_props::mass_props(&m, out[0]).expect("props").volume;
        let want = 8.0 - std::f64::consts::PI * 0.25 * 1.0; // the bore reaches z = 1
        assert!((v - want).abs() < 1e-9, "volume {v} vs {want}");
    }

    /// ★ **The uniform-slab theorem's counterexample, as a band case.** A cylinder standing in an
    /// L-prism's notch is clear of every wall by more than `r` — the refuted "z-range" rule would
    /// keep its wall — but it is outside the material, so `Cut` changes nothing and no band
    /// survives.
    #[test]
    fn a_cylinder_in_the_notch_keeps_no_band() {
        let profile = crate::ops::Profile2d::polygon(vec![
            nacre_math::Point2::from_array([0.0, 0.0]),
            nacre_math::Point2::from_array([2.0, 0.0]),
            nacre_math::Point2::from_array([2.0, 1.0]),
            nacre_math::Point2::from_array([1.0, 1.0]),
            nacre_math::Point2::from_array([1.0, 2.0]),
            nacre_math::Point2::from_array([0.0, 2.0]),
        ])
        .unwrap();
        let mut m = crate::ops::replay(&[crate::tests::extrude_log_op(profile, 1.0)]).unwrap();
        let a = m.live_solids[0];
        // In the notch (x,y ∈ [1,2]²), a slim drill standing clear of both notch walls.
        let b = m.add_cylinder(
            Point3::from_array([1.5, 1.5, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.25,
            3.0,
        );
        m.rebuild_adjacency();
        let (out, _) = bands(&m, a, b, BoolKind::Cut);
        assert!(
            out.is_empty(),
            "the notch is void — nothing to cut: {out:?}"
        );

        // ★ …and end to end: `Cut` returns the prism untouched. That is the counterexample's real
        // assertion — a band pass answering by z-range would run a wall through the notch, and
        // this volume would move.
        let before = nacre_props::mass_props(&m, a).expect("props").volume;
        let cut = crate::boolean(&mut m, BoolKind::Cut, a, b).expect("a cut that changes nothing");
        assert_eq!(cut.len(), 1);
        let after = nacre_props::mass_props(&m, cut[0]).expect("props").volume;
        assert!(
            (after - before).abs() < 1e-12,
            "the notch drill removes nothing: {before} → {after}"
        );
    }
}
