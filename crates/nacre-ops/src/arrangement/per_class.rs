use super::*;
/// **The per-class stages, in one place** — the walk and the nesting.
///
/// **It exists because the sequence would otherwise be written twice.** `trace_result_faces`'
/// `arrange` and `frame_audit`'s replay run the same stages against the same edges, and the second
/// is what makes a reject debuggable — so a drift between them would be, in `decline_to_reject`'s
/// words, *"the worst possible time to be lying"*. One copy, and the audit is auditing what the
/// boolean ran.
///
/// **And it puts every stage behind one `Result`.**
///
/// **`emit_faces` is the last stage**, and it fails only by the name a whole circle's marks are
/// refused with (`MergedCircle::whole_marks`) — the name `label_cells` meets first, reading the
/// same marks for every circle's mask; what follows it only reads its faces (the along-the-axis
/// check).
pub(super) fn per_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    kind: BoolKind,
    wc: usize,
    edges: &ClassEdges<'_>,
) -> Result<Staged, BoolError> {
    let (cells, face_of) = timed!(C_EXTRACT, walk_cells(jd, cyls, wc, edges))?;
    let nesting = timed!(C_NEST, nest_cells(jd, cyls, wc, &cells, edges))?;
    // **The seed is `[false; 4]`, and the argument is why it stays an argument.** The
    // arrangement covers all of space, so its unbounded cells reach infinity, where neither solid
    // is. That is a fact about arranging the *whole* model — restrict the input to a region of
    // space and the unbounded cells become an artifact of the restriction, which is what the
    // parameter records.
    let labels = timed!(
        C_LABEL,
        label_cells(&cells, &face_of, edges, &nesting, [false; 4])
    )?;
    // The cells beside a knife carry one label: its mask is zero.
    let knife_cuts = (edges.knives_from..edges.pieces.len())
        .map(|i| {
            face_of
                .get(&(2 * i))
                .is_some_and(|&c| knife_parts(kind, &edges.pieces[i].merged, labels[c]))
        })
        .collect();
    let (faces, disk_labels, arc_labels) = timed!(
        C_EMIT,
        emit_faces(kind, &labels, &cells, edges, jd, wc, &nesting.holes)
    )?;
    // **A face bounded along a ruling of a section the class draws along the axis**
    // (`Silence::Sectioned`, a class with no world statement): the gate proved every face on the
    // class apart from the lateral boxes, so such a face would have no lateral partner — the stages
    // disagree. Checked here, where the faces are. (A section oblique to the axis is a line the
    // class draws like any other; a face bounded along it carries an elliptical edge.)
    use crate::combinatorics::Wall;
    let along = |cyl: usize| cyls.get(cyl).is_some_and(|cy| cy.sections_along(wc));
    if faces.iter().any(|f| {
        std::iter::once(&f.outer).chain(&f.inner).any(|b| {
            b.ring().is_some_and(|r| {
                r.walls
                    .iter()
                    .any(|w| matches!(*w, Wall::Ruling { cyl, .. } if along(cyl)))
            })
        })
    }) {
        return Err(reject(RejectReason::CylinderStagesDisagree));
    }
    // The cylinder chart's **vertical** lines, read off the pieces this class already made. Where
    // a ruling ends on the chart is the chart's to read, off the ends' names (`cyl_chart::chart_of`).
    let ruling_extents: Vec<(usize, RulingExtent)> = edges
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let OnCurve::Ruling { cyl, side, .. } = p.on else {
                return None;
            };
            // A sectioned class's ruling is not a line of the lateral's chart: no face of the class
            // reaches the lateral, so it splits no cell there (see the check above).
            if along(cyl) {
                return None;
            }
            // The chart's vertical answer — see [`RulingExtent::label`] and
            // [`RulingExtent::outer`], the two cells beside the ruling. Each is an `Option` rather
            // than a `?` deliberately: dropping the piece here would make the recorded population
            // differ between a test build and a release one, and then the census would be measuring
            // a corpus production never sees.
            let (label, outer) = {
                let cell_at = |even: bool| -> Option<Label> {
                    Some(labels[*face_of.get(&(2 * i + usize::from(!even)))?])
                };
                // A **tangent** ruling (`side == 0`) has no interior side on this class:
                // the plane touches the cylinder, and both cells beside the line lie outside its
                // disk. Its extent carries no label — the chart reads it as a station where the
                // face ends, not as a wall with a chamber behind it.
                let inside = if side == 0 {
                    None
                } else {
                    ruling_interior_is_even(jd, wc, side)
                };
                let out = inside.and_then(cell_at);
                let outer = inside.and_then(|inside| cell_at(!inside));
                #[cfg(test)]
                {
                    // The postcondition, checked rather than assumed — see `ruling_probe::SIDE_CHECK`.
                    // It asks whether a **cell** carries the lateral's solid, not whether the two
                    // cells *differ*: crossing a ruling on this wall crosses the **lateral**, so they
                    // differ always (measured: every piece, without
                    // exception) and that says nothing about the side.
                    // Recorded only where a label was formed, so `None` means one thing — the content
                    // did not distinguish — and never "there was nothing to check".
                    if let Some(s) = p.merged.first().map(|m| m.0) {
                        let has = |l: Label| {
                            let b = match s {
                                crate::planes::SolidSide::A => 0,
                                crate::planes::SolidSide::B => 2,
                            };
                            l[b] || l[b + 1]
                        };
                        if let (Some(i_l), Some(e_l)) = (out, outer) {
                            // The lateral's material is inside the cylinder for a boss and
                            // **outside** for a bore or a notch (`orient` — re-operated Cut
                            // results are such a population on a through-axis class), so
                            // the content agrees with the derivation when the cell on the material's
                            // side has the solid and the other does not.
                            let material_inside =
                                matches!(p.on, OnCurve::Ruling { orient, .. } if orient > 0);
                            let verdict = match (has(i_l), has(e_l)) {
                                (a, b) if a != b => Some(a == material_inside),
                                _ => None,
                            };
                            // **Asserted where the fact is made**, so the coverage is total and
                            // the panic names the offending test — a test reading the ledger afterwards
                            // sees only what ran before it.
                            assert_ne!(
                                verdict,
                                Some(false),
                                "{}: class {wc}, side {}",
                                ruling_probe::WRONG_SIDE,
                                side
                            );
                            ruling_probe::SIDE_CHECK.push(verdict);
                        }
                    }
                    // A tangent piece (`side == 0`) carries no label by design; the
                    // ledger's proposition is about the rulings a chamber lies behind.
                    if side != 0 {
                        ruling_probe::LABELLED.push(out.is_some());
                    }
                }
                (out, outer)
            };
            Some((
                cyl,
                RulingExtent {
                    wall: wc,
                    side,
                    end: p.end,
                    label,
                    outer,
                    #[cfg(test)]
                    marks: p.merged.clone(),
                },
            ))
        })
        .collect();
    Ok(Staged {
        #[cfg(test)]
        cells,
        #[cfg(test)]
        nesting,
        #[cfg(test)]
        labels,
        faces,
        disk_labels,
        arc_labels,
        ruling_extents,
        knife_cuts,
    })
}

/// **A knife parts the kept material beside it** — the rule that makes a knife an edge of the
/// class's faces. The result's faces on a plane end where the material they bound parts, and a
/// label change says so for every edge but a knife: the knife's solid `X` has both faces leaving
/// the plane to one side `S`, so the label is one on both sides of it, and what differs is the
/// wedge between the two faces — the chamber on `S` with `X`'s bit there turned over (two faces of
/// one solid, crossed in turn). The cells beside the knife are a face where `keep` differs between
/// their two chambers; where the kept one is `S`'s and the wedge is not kept, the kept material on
/// `S` reaches the plane on both sides of the knife and meets itself only along it, so the face
/// ends there — two wings resting on one floor, the floor two faces. A knife whose wedge is the
/// kept one, or whose kept chamber is the other side's, leaves the face whole: the material it
/// bounds passes under the knife.
///
/// `label` is the cells' beside the knife — one label on both its sides; `merged` the knife's
/// contributions, which are grazes only ([`is_knife`]). A knife two solids graze at once parts
/// nothing ([`knife_grazer`]): the premise below holds for each of them, so all four bits are one
/// on both sides of it, `keep` is too, and the class has no face there to part.
///
/// **`X`'s two bits in `label` are equal**: `X` has no face on the plane where its knife lies — a
/// face there would hold `X`'s own edge inside it, which no valid solid has.
fn knife_parts(kind: BoolKind, merged: &[(SolidSide, SegKind)], label: Label) -> bool {
    let Some((x, side)) = knife_grazer(merged) else {
        debug_assert!(
            label[0] == label[1] && label[2] == label[3],
            "a knife both solids graze, where one has a face on the plane: {label:?}"
        );
        return false;
    };
    debug_assert_eq!(
        label[label_bit(x, true)],
        label[label_bit(x, false)],
        "a knife of {x:?} where {x:?} has a face on the plane: {label:?}"
    );
    let kept = |l: Label, above: bool| {
        crate::draft::keep(
            kind,
            l[label_bit(SolidSide::A, above)],
            l[label_bit(SolidSide::B, above)],
        )
    };
    let mut wedge = label;
    wedge[label_bit(x, side)] ^= true;
    // Counter-clockwise from the plane's half on one side of the knife, turning to `side`: the
    // chamber on `side`, the knife's two faces round the wedge between them, the chamber again,
    // the plane's other half, and the far side.
    let (plane, (first, second)) = (0, (1, 2));
    let star = crate::draft::LineStar::new(
        vec![plane, first, second, plane],
        vec![
            kept(label, side),
            kept(wedge, side),
            kept(label, side),
            kept(label, !side),
        ],
    );
    star.continues(plane) == Some(false)
}

/// The solid whose knife this is and the side its faces leave to — `None` unless every
/// contribution is a graze of one solid.
fn knife_grazer(merged: &[(SolidSide, SegKind)]) -> Option<(SolidSide, bool)> {
    let mut grazer: Option<(SolidSide, bool)> = None;
    for &(solid, k) in merged {
        let SegKind::Graze { body_above } = k else {
            return None;
        };
        match grazer {
            None => grazer = Some((solid, body_above)),
            Some((x, _)) if x != solid => return None,
            Some(_) => {}
        }
    }
    grazer
}

/// The label's bit for one solid's chamber on one side — [`edge_mask`]'s layout.
fn label_bit(solid: SolidSide, above: bool) -> usize {
    let base = match solid {
        SolidSide::A => 0,
        SolidSide::B => 2,
    };
    base + usize::from(!above)
}

/// **Can a knife part the material under this `kind` at all** — [`knife_parts`] over every label
/// a knife can lie under: its solid's two bits equal (the premise [`knife_parts`] states), the
/// other solid's free. Under `Fuse` none can — the one kept chamber beside a face keeps the wedge
/// too — and neither can a knife two solids graze at once.
fn knife_can_part(kind: BoolKind, merged: &[(SolidSide, SegKind)]) -> bool {
    let Some((x, _)) = knife_grazer(merged) else {
        return false;
    };
    (0u8..16)
        .map(|bits| -> Label { std::array::from_fn(|k| bits >> k & 1 == 1) })
        .filter(|l| l[label_bit(x, true)] == l[label_bit(x, false)])
        .any(|l| knife_parts(kind, merged, l))
}

/// **One class's faces, with the knives that part its material** — [`per_class`](fn@per_class)
/// over the split's pieces and the knives that are edges of its faces ([`knife_parts`]).
///
/// Which knives part the material is read off the labels of the cells beside them, so the walk
/// that carries a knife is the one that answers it: the knives that can be edges at all are walked
/// with the pieces, and the class is walked again with only those that part the material when some
/// do not ([`rewalk`]). A knife that cannot be an edge — one that ends inside a face, or one
/// joining two pieces of the skeleton no path of edges joins ([`knives_as_edges`]) — is not walked:
/// a knife that parts the material there is the result touching itself, the material passing round
/// its end or between the two boundaries, which the self-touch verdict names. A class with no knife
/// to carry is walked once.
///
/// **The second walk's knives all part.** A cell's label is the seed turned over by the masks on a
/// path to it, and a knife's mask is zero — so the label beside a knife is the same whichever
/// other knives are walked, and what parted in the first walk parts in the second.
pub(super) fn arrange_class(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    kind: BoolKind,
    wc: usize,
    split: &ClassSplit,
) -> Result<Staged, BoolError> {
    let carried = knives_as_edges(
        split,
        (0..split.knives.len())
            .filter(|&i| knife_can_part(kind, &split.knives[i].merged))
            .collect(),
    );
    let plain = || per_class(jd, cyls, kind, wc, &ClassEdges::of(cyls, wc, split));
    if carried.is_empty() {
        return plain();
    }
    let walk = |knives: &[usize]| {
        let pieces: Vec<Piece> = split
            .pieces
            .iter()
            .cloned()
            .chain(knives.iter().map(|&i| split.knives[i].clone()))
            .collect();
        per_class(
            jd,
            cyls,
            kind,
            wc,
            &ClassEdges::with_knives(cyls, wc, split, &pieces),
        )
    };
    let all = walk(&carried)?;
    match rewalk(split, &carried, &all.knife_cuts) {
        Rewalk::Kept => Ok(all),
        Rewalk::Plain => plain(),
        Rewalk::Again(kept) => {
            let again = walk(&kept)?;
            debug_assert!(
                again.knife_cuts.iter().all(|&c| c),
                "a knife that parted in the first walk does not in the second"
            );
            Ok(again)
        }
    }
}

/// What [`arrange_class`] does after the walk that carried every knife that can be an edge.
#[derive(Debug, PartialEq, Eq)]
enum Rewalk {
    /// Every knife walked parts the material — that walk is the answer.
    Kept,
    /// No knife both parts the material and can be an edge — the class is walked without knives.
    Plain,
    /// Walk again with these knives: those that part the material and can still be edges.
    Again(Vec<usize>),
}

/// **The knives the class is walked with: those that part the material and can be edges** —
/// `carried` the knives the first walk carried, `knife_cuts` its answer per knife. The knives that
/// part go through [`knives_as_edges`] again: one of them can be left with a free end when a knife
/// beside it does not part.
fn rewalk(split: &ClassSplit, carried: &[usize], knife_cuts: &[bool]) -> Rewalk {
    let parting: Vec<usize> = carried
        .iter()
        .zip(knife_cuts)
        .filter_map(|(&i, &cuts)| cuts.then_some(i))
        .collect();
    let kept = knives_as_edges(split, parting);
    if kept == carried {
        Rewalk::Kept
    } else if kept.is_empty() {
        Rewalk::Plain
    } else {
        Rewalk::Again(kept)
    }
}

/// **The knives among `candidates` that can be edges of the class's faces** — indices into
/// `split.knives`: those whose every end is a vertex of the skeleton or of another such knife, and
/// whose knives join no two pieces of the skeleton that are not already one. A knife with a free
/// end ends inside a face and cannot cut it; a chain of knives between two separate pieces of the
/// skeleton (a face's outer boundary and a hole's) joins them rather than cutting a face in two.
/// The answer is its own fixed point: a chain dropped as a bridge shares no vertex with the chains
/// kept.
fn knives_as_edges(split: &ClassSplit, candidates: Vec<usize>) -> Vec<usize> {
    let mut carried = candidates;
    if carried.is_empty() {
        return carried;
    }
    // Vertex ids over the skeleton's ends and the knives'.
    let mut id: HashMap<NodeId, usize> = HashMap::new();
    let mut of = |n: NodeId| {
        let next = id.len();
        *id.entry(n).or_insert(next)
    };
    let skel: Vec<[usize; 2]> = split.pieces.iter().map(|p| p.end.map(&mut of)).collect();
    let knife: Vec<[usize; 2]> = split.knives.iter().map(|p| p.end.map(&mut of)).collect();
    let n = id.len();
    let mut deg_skel = vec![0usize; n];
    for e in &skel {
        for &v in e {
            deg_skel[v] += 1;
        }
    }
    // A free end, again and again: removing a knife can free its neighbour's end.
    loop {
        let mut deg = deg_skel.clone();
        for &k in &carried {
            for &v in &knife[k] {
                deg[v] += 1;
            }
        }
        let before = carried.len();
        carried.retain(|&k| knife[k].iter().all(|&v| deg[v] >= 2));
        if carried.len() == before {
            break;
        }
    }
    // Which piece of the skeleton each vertex is on, and the knives' own chains.
    let mut skel_part: Vec<usize> = (0..n).collect();
    for e in &skel {
        let (a, b) = (
            crate::planes::uf_find(&mut skel_part, e[0]),
            crate::planes::uf_find(&mut skel_part, e[1]),
        );
        skel_part[a] = b;
    }
    let mut chain: Vec<usize> = (0..n).collect();
    for &k in &carried {
        let (a, b) = (
            crate::planes::uf_find(&mut chain, knife[k][0]),
            crate::planes::uf_find(&mut chain, knife[k][1]),
        );
        chain[a] = b;
    }
    // Per chain, the skeleton pieces its ends rest on; a chain resting on two joins them.
    let mut rests: HashMap<usize, Vec<usize>> = HashMap::new();
    for &k in &carried {
        for &v in &knife[k] {
            let c = crate::planes::uf_find(&mut chain, v);
            let on = rests.entry(c).or_default();
            if deg_skel[v] > 0 {
                on.push(crate::planes::uf_find(&mut skel_part, v));
            }
        }
    }
    carried.retain(|&k| {
        let c = crate::planes::uf_find(&mut chain, knife[k][0]);
        rests[&c].windows(2).all(|w| w[0] == w[1])
    });
    carried
}

/// **Does this circle actually meet that ruling?** — asked of the class's own tables: `Ok` is the
/// measured answer ([`crossing_measured`]), `Err` the name a pair that could not be measured is
/// refused by. On a class sectioning the circle's cylinder obliquely
/// ([`crate::planes::CylRecord::Sections`]) that is the oblique coverage's end
/// (`ObliqueCylinderCut` — an ellipse against a ruling is not measured yet), and anywhere else the
/// net's own backstop (`CircleCrossesRuling`). Named here, once, so every way the measure can fail
/// carries the same name.
///
/// This is the check `RejectReason::CircleCrossesRuling` is raised from, so it ships.
fn circle_crosses_ruling(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    circle: &MergedCircle,
    ruling: &MergedRuling,
) -> Result<bool, RejectReason> {
    let unmeasured = if cyls[circle.cyl].sections_oblique(wc) {
        RejectReason::ObliqueCylinderCut
    } else {
        RejectReason::CircleCrossesRuling
    };
    crossing_measured(jd, cyls, wc, circle, ruling).ok_or(unmeasured)
}

/// **[`circle_crosses_ruling`]'s measure.**
///
/// **First, the pair rule's own question, of the two edges as pieces**
/// ([`crate::planes::pieces_clear`]): an arc is the piece of the axis band its curve sweeps × its
/// arc on its cylinder — `[t, t]` for a circle, `t ± Δ` for an oblique section's ellipse — and a
/// ruling its axis span × the whole circle on its own (one ruling's angle is seldom rational, and
/// the whole circle is wider, never narrower). So the net and the gate ask one predicate, and
/// what it is asked about is the bounded piece — the ruling's segment, not the line it lies on.
/// Of its three witnesses only two can speak here: the class runs along the ruling's axis and
/// meets the circle's, so the two axes are never parallel and the parallel cross-section never
/// runs; on a class across the circle's axis the common perpendicular is the strip direction
/// below, where the ruling's whole-circle reach says less than the strip does. What this guards, then, is that the edges the
/// arrangement drew lie in the faces the gate proved apart (a whole circle drawn where a face used
/// a quarter of it is wider than the face).
///
/// Then the strip, for a circle: [`nacre_exact::cylinder_ruling_reached_extent`] — the arc's reach
/// across the strip, against the one named ruling's line — which reads a tangency as no crossing,
/// as this net must (a disjointness proof is about closed sets and cannot). The strip's arc ends
/// are a circle's ([`crate::planes::arc_ends_along`]), so an ellipse the pieces do not part is not
/// measured. `None` stays "could not be measured", never "no"; a piece that cannot be stated only
/// skips the first question.
fn crossing_measured(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    circle: &MergedCircle,
    ruling: &MergedRuling,
) -> Option<bool> {
    let coeffs = jd.world_rat(wc)?;
    let (circle_def, ruling_def) = (&cyls[circle.cyl].def, &cyls[ruling.cyl].def);
    let centre = combinatorics::circle_centre_rat(jd, wc, circle_def)?;
    let (o, m, r2) = (ruling_def.origin(), ruling_def.dir(), ruling_def.r2());
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // The predicate's own precondition, asked rather than assumed: a class that is not parallel
    // to the axis carries no ruling of it, and the decomposition would be about other geometry.
    if nacre_exact::dot_sign_rat(&n, &m) != nacre_exact::Orient::Zero {
        return None;
    }
    // The ruling as a piece: how far it reaches along the axis — each end stands on its section
    // ([`combinatorics::section_of`]), within the band that section sweeps, so the hull of the two
    // bands — the whole circle round. The arc's band is the class's own.
    let ruling_piece = (|| {
        let band = |n: NodeId| {
            crate::planes::trace_band(&combinatorics::section_of(jd, ruling_def, n)?, ruling_def)
        };
        let ([a0, a1], [b0, b1]) = (band(ruling.end[0])?, band(ruling.end[1])?);
        Some(crate::planes::Footprint {
            span: Some([a0.min(b0), a1.max(b1)]),
            theta: crate::planes::Theta::Whole,
        })
    })();
    let oblique = !crate::planes::class_carries_circle(&n, &circle_def.dir());
    let circle_band = crate::planes::trace_band(&coeffs, circle_def);
    let apart = |theta: crate::planes::Theta| -> bool {
        match (&ruling_piece, circle_band) {
            (Some(rp), Some(band)) => crate::planes::pieces_clear(
                circle_def,
                &crate::planes::Footprint {
                    span: Some(band),
                    theta,
                },
                ruling_def,
                rp,
            ),
            _ => false,
        }
    };
    // The strip runs along `e = n × m`, so that is the direction a piece states its reach in.
    let e = nacre_exact::cross3_rat(&n, &m)?;
    let ask = |theta: &crate::planes::Theta| -> Option<bool> {
        let (lo, hi) =
            crate::planes::arc_ends_along(&centre, circle_def.r2(), &circle_def.dir(), theta, &e)?;
        Some(nacre_exact::cylinder_ruling_reached_extent(
            &coeffs,
            &nacre_exact::StripReach {
                lo: (&nacre_exact::MeetPoint::Narrow(lo.0), &lo.1),
                hi: Some((&nacre_exact::MeetPoint::Narrow(hi.0), &hi.1)),
            },
            &o,
            &m,
            r2,
            // **The two vocabularies are opposite.** `ruling.side` is measured against `m × n`
            // (`ruling_side_signed`), and the scalar family writes its sides about `n × m`. The
            // restatement is one negation, and it belongs here — at the boundary between the two
            // index spaces — rather than inside a predicate that would then have to guess.
            -ruling.side,
            // **A touch is not what this asks about**. The node no road mints is a
            // **crossing**; an arc tangent to a ruling divides nothing, which is the same sentence
            // as for a tangency at a vertex. Measured: a slab wall at `x = ±12` puts a
            // `d 30` cylinder's rulings exactly on a `d 18` tool's circle, and the boolean that
            // follows validates clean.
            false,
        ))
    };
    // **The question is about the arc, not the circle**. A contribution states the
    // angular extent its face covers, and *that* is the edge — the rest of the circle is a
    // continuation no face uses. Asking the whole circle refuses a plate whose corner
    // fillets never come near the drill's rulings.
    //
    // An extent that cannot be realized, and a circle with nothing partial to read, fall back to
    // the whole circle: a superset of every arc on it, so the answer stays sound and is the one
    // this predicate always gave.
    if circle.merged.is_empty() {
        let whole = crate::planes::Theta::Whole;
        return if apart(whole) {
            Some(false)
        } else if oblique {
            None
        } else {
            ask(&whole)
        };
    }
    for (_, _, extent) in &circle.merged {
        let theta = extent.map_or(crate::planes::Theta::Whole, |ends| {
            combinatorics::rim_arc_of(jd, cyls, &centre, &circle_def.dir(), ends)
        });
        if apart(theta) {
            continue;
        }
        if oblique {
            return None;
        }
        if ask(&theta)? {
            return Some(true);
        }
    }
    Some(false)
}

/// How many `(circle, ruling)` pairs the net above looked at, and how many it **could not
/// measure** — the denominator of "zero crossings", without which the zero could mean
/// "the instrument saw nothing".
#[cfg(debug_assertions)]
pub(crate) struct MixedClassAudit {
    pub(crate) pairs: std::sync::atomic::AtomicUsize,
    pub(crate) unmeasured: std::sync::atomic::AtomicUsize,
}

#[cfg(debug_assertions)]
pub(crate) static MIXED_CLASS_AUDIT: MixedClassAudit = MixedClassAudit {
    pairs: std::sync::atomic::AtomicUsize::new(0),
    unmeasured: std::sync::atomic::AtomicUsize::new(0),
};

/// Run the arc split and, when the class carries rulings, the ruling split, and hold the result
/// as one list of pieces and the knives apart ([`ClassSplit`]): **the one place a class's edges are
/// split** — a step of the alias rounds
/// ([`class_split`]), so a tie between a crossing it names and a name already there is recorded
/// in `aliases` ([`Aliases::record_same_point`]) for the next round rather than refused.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_class_edges(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    wc: usize,
    segs: &[MergedSeg],
    circles: &[MergedCircle],
    rulings: &[MergedRuling],
    aliases: &mut Aliases,
) -> Result<ClassSplit, BoolError> {
    use std::borrow::Cow;
    // **A class may carry a circle (⊥ one cylinder's axis) and rulings (∥ another's)
    // at once, and that is not a refusal**. What would be unarrangeable is the two
    // **meeting**: neither split cuts the other's kind — [`split_circles`] cuts circles
    // against segments, [`split_rulings`] cuts rulings against segments — so a crossing they
    // both walked past would be a node no road ever mints, and the point itself is
    // `plane ∩ cylinder ∩ cylinder`, a nested radical with no name.
    //
    // **That population cannot reach here, by two layers.** A circle only becomes an edge
    // *whole* (a partial trace forces its own cut, and `PartialCircleUncut` names the
    // alternative) and a ruling is a segment of a **face**, so two such edges meeting means
    // two lateral **faces** share a point — exactly what [`crate::planes::lateral_faces_clear`]
    // is written to deny ("the two classes share no face — the proposition the arrangement
    // needs"). A pair not shown to clear is `CylinderPairContact`, raised before any
    // arrangement runs. That covers pairs from **different** solids; the pair loop skips
    // same-solid pairs on purpose. But a same-solid pair cannot be built either: two cylinder
    // faces whose axes are *not* parallel need two extrude directions, hence a boolean, and
    // at that boolean they are on opposite sides and meet the very gate above. (Two arcs of one
    // sketch do make one body whose cylinders meet — with **parallel** axes, which never put a
    // circle and a ruling on one class, and such an operand is refused at the boolean's entry
    // anyway, `CylinderPairEdge`.) A circle and rulings that do not meet make an ordinary
    // class: the cap of a pin nested crosswise in a larger cylinder holds the pin's circle
    // between the larger one's rulings, and the cap a crosswise stud stands proud of holds the
    // stud's rulings inside its circle.
    //
    // **So the check below is a net over an argument, not a filter over a population**,
    // and it ships rather than sitting under `debug_assertions`: the day the cylinder-pair
    // refusal opens, this population becomes real, and the obligation to carry a shipped check
    // then is written at that refusal, where it would be read.
    //
    // **The question is asked of the two edges as pieces** — the arc a contribution states,
    // not the whole circle (which is not an edge of any face when the face uses a quarter of
    // it, and reading it that way refuses a plate whose corner fillets never come near a
    // crosswise drill), and the ruling's segment, not its line (whose line runs past the end
    // of a short crosswise cylinder and through a cap circle its segment never reaches). Asked
    // so, it is the pair rule's own question ([`crate::planes::pieces_clear`]), and the
    // argument holds edge for edge: a pair the gate proved apart has its edges apart.
    //
    // **One sentence covers all three uncut pairs.** Circle×circle, circle×ruling and
    // ruling×ruling are cut by nothing, and all three are safe for the same reason: **a
    // class's edge is some face's boundary** (or, for an oblique section, a curve on a
    // lateral face), so two of them from different cylinders meeting puts two lateral
    // faces on one point. Only this pair carries a check, because only this
    // pair has a cheap exact predicate; building detectors for populations that cannot arrive
    // measures nothing.
    //
    // A pair that could not be measured refuses too — by the net's name, or on an oblique
    // section by the coverage's (`circle_crosses_ruling`): honest-reject over silent-wrong, and
    // measured to fire on nothing today (`unmeasured` is 0 across the corpus, the census and
    // the kit).
    if !rulings.is_empty() && !circles.is_empty() {
        for c in circles {
            for ru in rulings {
                let crosses = circle_crosses_ruling(jd, cyls, wc, c, ru);
                #[cfg(debug_assertions)]
                {
                    use std::sync::atomic::Ordering;
                    MIXED_CLASS_AUDIT.pairs.fetch_add(1, Ordering::Relaxed);
                    if crosses.is_err() {
                        MIXED_CLASS_AUDIT.unmeasured.fetch_add(1, Ordering::Relaxed);
                    }
                }
                match crosses {
                    Ok(false) => {}
                    Ok(true) => return Err(reject(RejectReason::CircleCrossesRuling)),
                    Err(reason) => return Err(reject(reason)),
                }
            }
        }
    }
    type AfterCircles<'a> = (
        Cow<'a, [MergedSeg]>,
        Vec<Piece>,
        Cow<'a, [MergedCircle]>,
        Vec<(usize, CutRim)>,
    );
    // **Circles first, then rulings — one step per curve kind, in this order.** The step that
    // refuses first names the boolean, so which name a refused input carries rests on it.
    let (segs2, arcs, circles2, cut_rims): AfterCircles<'_> =
        match split_circles(jd, cyls, wc, segs, circles, aliases)? {
            Some((s, c, a, r)) => {
                // The names below rest on this: `split_circles` answers `Some` only when some
                // circle collected a node — an arc's end or a crossing — and either yields at
                // least one arc (a tangency is skipped before it is collected). So "was split"
                // and "has arcs" are one fact.
                debug_assert!(!a.is_empty(), "a split that cut no circle into arcs");
                (Cow::Owned(s), a, Cow::Owned(c), r)
            }
            None => (
                Cow::Borrowed(segs),
                Vec::new(),
                Cow::Borrowed(circles),
                Vec::new(),
            ),
        };

    let (segs3, rulings2) = match split_rulings(jd, cyls, faces, wc, &segs2, rulings, aliases)? {
        Some((s, r)) => (Cow::Owned(s), Cow::Owned(r)),
        None => (segs2, Cow::Borrowed(rulings)),
    };
    // The knives cut every curve above and are set apart here, last: whether one is an edge of the
    // class's faces is asked of the material beside it (`arrange_class`).
    let mut knives = Vec::new();
    let mut segs4 = Vec::with_capacity(segs3.len());
    for sg in segs3.into_owned() {
        if is_knife(&sg.merged)? {
            knives.push(Piece::from(sg));
        } else {
            segs4.push(sg);
        }
    }
    // Segments, then arcs, then ruling pieces — the order [`ClassEdges`] numbers them in.
    let pieces = segs4
        .into_iter()
        .map(Piece::from)
        .chain(arcs)
        .chain(rulings2.into_owned().into_iter().map(Piece::from))
        .collect();
    Ok(ClassSplit {
        pieces,
        circles: circles2.into_owned(),
        cut_rims,
        knives,
    })
}

impl<'a> ClassEdges<'a> {
    /// The class's split edges as the cells read them — one numbering over pass A's product.
    pub(super) fn of(
        cyls: &'a [crate::planes::WorkingCyl],
        wc: usize,
        split: &'a ClassSplit,
    ) -> Self {
        Self::with_knives(cyls, wc, split, &split.pieces)
    }

    /// The split's edges with knives after them: `pieces` is `split.pieces` followed by the knives
    /// the walk is to carry ([`arrange_class`]).
    pub(super) fn with_knives(
        cyls: &'a [crate::planes::WorkingCyl],
        wc: usize,
        split: &'a ClassSplit,
        pieces: &'a [Piece],
    ) -> Self {
        debug_assert!(pieces.len() >= split.pieces.len());
        ClassEdges {
            pieces,
            cyls,
            circles: &split.circles,
            wc,
            knives_from: split.pieces.len(),
        }
    }

    /// Where the walk's half-edges end and the circles' pseudo-half-edges begin.
    pub(super) fn he_count(&self) -> usize {
        2 * self.pieces.len()
    }

    /// **The one classifier.** Every reader of a half-edge's kind asks this.
    pub(super) fn kind(&self, he: usize) -> HalfEdgeKind {
        if he < self.he_count() {
            HalfEdgeKind::Piece(he / 2)
        } else {
            HalfEdgeKind::Circle((he - self.he_count()) / 2)
        }
    }

    /// The vertex a half-edge leaves. `he % 2 == 0` takes `end[0]` — one rule, all pieces.
    pub(super) fn origin(&self, he: usize) -> NodeId {
        match self.kind(he) {
            HalfEdgeKind::Piece(i) => self.pieces[i].end[he % 2],
            HalfEdgeKind::Circle(_) => {
                unreachable!("a circle's pseudo-half-edge has no vertex to leave")
            }
        }
    }

    /// **The one place a `RingEdge` is made from a half-edge.** The carrier, sense included, is the
    /// piece's ([`Piece::carrier`]) — a second spelling that passed `Carrier::plane(wall)`, that is
    /// `sense: None`, would drop the sense a split carries onto its sub-segments, the one fact the
    /// endpoints cannot supply.
    pub(super) fn edge_at(&self, he: usize) -> combinatorics::RingEdge {
        let HalfEdgeKind::Piece(i) = self.kind(he) else {
            unreachable!("a circle's pseudo-half-edge is not a ring edge")
        };
        let p = &self.pieces[i];
        // The endpoints as handles on this edge's line, carried by the producer. Not recovered from
        // the names: a canonical name need not mention `wc` or the wall (see
        // `combinatorics::RingEdge`).
        let pins = p.pins();
        combinatorics::RingEdge {
            node: self.origin(he),
            to: self.origin(he ^ 1),
            carrier: p.carrier(self.cyls, he % 2 == 0, self.wc),
            from_h: pins[he % 2],
            to_h: pins[1 - he % 2],
        }
    }
}

#[cfg(test)]
#[path = "../tests/per_class.rs"]
mod tests;
