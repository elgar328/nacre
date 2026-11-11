use super::*;
/// **The population gate** — decides, exactly, whether this operand pair stays inside the
/// population the cylinder arrangement serves — plane classes across the axis, along it, or
/// oblique to it in whole sections — and names the refusal otherwise. Signs over world-stated
/// descriptions are answered totally by `nacre-exact` (denominators cleared, `BigInt`); a class
/// with no world statement is sided by the judge's filter (`WorkingPlane::misses`), which proves
/// or abstains. The values the gate forms are checked `Rat` on world-stated descriptions, and a
/// refusal says which of two things stopped it: a
/// description that cannot be put in the world in one frame is geometry,
/// [`RejectReason::CylinderGateUndecided`]; one the world holds and `Rat` does not is width,
/// [`RejectReason::WitnessNotRational`].
///
/// Per (plane class, cylinder) pair, [`class_meet`] answers what the class receives
/// ([`ClassMeet`]), and every pair is answered before anything reads the answers. The gate is
/// then what can and cannot be drawn: it settles the sectioned candidates against the cylinder's
/// other records, refuses an oblique section the class cannot draw and a class carrying the
/// traces of axes that are not parallel ([`RejectReason::ObliqueCylinderCut`]), and hands the
/// arrangement the records it reads — the cuts, oblique sections, touchings and silences
/// ([`WorkingCyl::records`]), the [`Tangency`] rows and the lines two classes share on a lateral
/// face ([`SharedRuling`]) — refusing a tangency row it could not state.
///
/// Per cylinder pair **of different owners** (two classes of one valid solid keep their faces
/// apart by construction and are not asked): axes clear of each other (`dist > r₁+r₂`, whatever
/// their orientation) pass outright; one surface stated under two handles ([`same_surface`]) is
/// refused as the coincident pair it is; for parallel axes one cylinder strictly inside the other
/// clears outright ([`nacre_exact::cylinders_nested`]); otherwise the pair passes when every
/// **face** of one class provably misses every face of the other ([`lateral_faces_clear`] — each
/// face pair apart along a direction, or one face inside the other cylinder, or for parallel axes
/// apart in the cross-section: [`pieces_clear`]) — and a pair whose faces cannot be shown to miss
/// is recorded and refused as [`RejectReason::CylinderPairContact`] (the quartic intersection
/// curve is not built).
#[allow(clippy::type_complexity)]
pub(crate) fn cylinder_gate(
    model: &Model,
    cyl_surfs: &[Handle<Surface>],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
    jd: &Judge<'_, WorkingPlane>,
) -> Result<(Vec<WorkingCyl>, Vec<Tangency>), BoolError> {
    let geom = jd.planes;
    // Every question below is a **sign**, and the `nacre-exact` layer answers signs totally: local
    // checked-`Rat` closures would decline on overflow, putting a width limit inside
    // `CylinderGateUndecided` and making that name say less than it claims.
    use nacre_exact::Orient;
    let undecided = || reject(RejectReason::CylinderGateUndecided);

    let mut tangencies: Vec<Tangency> = Vec::new();
    let mut cyls = Vec::with_capacity(cyl_surfs.len());
    for &surf in cyl_surfs {
        // **The world statement or nothing.** A moved cylinder's def is written before its
        // motion, and every test below compares it against world planes. A chain that folds
        // carries it out exactly (`Model::world_cylinder_def` — the same door the face rows
        // take); a frame node or a turn off the quarters has no world description here and is
        // refused rather than measured across two frames, and an overflow carrying it is refused
        // for its width ([`world_cylinder`] tells the two apart).
        let def = world_cylinder(model, surf).map_err(|e| {
            reject(
                e.stated_reason()
                    .expect("a cylinder statement fails for its frame or its width"),
            )
        })?;
        let Some(cache) = model.surface_cache(surf).cylinder() else {
            unreachable!("push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache")
        };
        // **One surface, two solids.** Cylinders intern by their exact statement, so two operands
        // whose laterals coincide arrive as one class carrying rows of both — the coaxial pair of
        // equal radius, by another spelling. The chart reads a class as one solid's surface
        // (`cyl_chart::chart_of`), so this is refused here, by the name the pair rule gives that
        // pair, instead of asserting inside the chart: two identical circle prisms fused would
        // reach that assertion. Otherwise the side that owns the class is written on it — the pair
        // loop reads it.
        let owned = |rows: &[FaceRow]| {
            rows.iter()
                .any(|r| matches!(r, FaceRow::Cylinder(cf) if cf.surf == surf))
        };
        let owner = match (owned(&faces[..n_a]), owned(&faces[n_a..])) {
            (true, true) => return Err(reject(RejectReason::CylinderPairContact)),
            (true, false) => SolidSide::A,
            (false, true) => SolidSide::B,
            // A class is named by some face row; a table with none is the class table and the
            // face rows disagreeing.
            (false, false) => return Err(reject(RejectReason::CylinderStagesDisagree)),
        };
        cyls.push(WorkingCyl {
            surf,
            def,
            realized: *cache,
            owner,
            shared: Vec::new(),
            grazed: Vec::new(),
            records: Vec::new(),
            boxes: Vec::new(),
        });
    }

    // **One row per `cyl_surfs` entry, in order** — the loop above either pushes or returns, so
    // this table is indexed by the very `ClassIx::Cyl` number that named the surface. Consumers
    // index it directly (`merge_circles`), which is only sound while that holds.
    debug_assert_eq!(
        cyls.len(),
        cyl_surfs.len(),
        "the cylinder class table is index-aligned with the class numbering"
    );

    // **A class with no world statement is asked first whether the whole cylinder stands on one
    // side of it** (`nameless_stance`) — the question every other arm of `class_meet` asks
    // needs the class's rational coefficients, and this one needs none. Then the cylinder has
    // nothing to say on that class, and the class nothing to say about the cylinder; the roads that
    // would ask read the record (`WorkingCyl::records`). The boxes are built only when such a
    // class is present — the common boolean pays nothing.
    let unstated = geom.iter().any(|wp| wp.world.is_none());
    let mut whole_rims = vec![false; cyls.len()];
    if unstated {
        for (cyl, whole) in cyls.iter_mut().zip(&mut whole_rims) {
            let fps = lateral_footprints(faces, cyl.surf);
            cyl.boxes = fps
                .iter()
                .map(|fp| lateral_box(&cyl.def, fp))
                .collect::<Option<Vec<_>>>()
                .unwrap_or_default();
            *whole = only_whole_rims(model, faces, cyl.surf);
        }
    }
    let cx = PairCx {
        model,
        faces,
        plane_ix,
        n_a,
        jd,
        whole_rims: &whole_rims,
    };
    let mut meets: Vec<Vec<ClassMeet>> = Vec::with_capacity(cyls.len());
    for (ci, cyl) in cyls.iter().enumerate() {
        // The footprint's second axis, gathered **lazily and at most once** per cylinder: it does
        // not depend on the plane class the loop below walks, but almost no boolean ever asks for
        // it — the plane-level test decides first (`wall_faces_clear`'s doc has the measurement).
        let mut footprints: Option<Vec<Footprint>> = None;
        let mut row = Vec::with_capacity(geom.len());
        for c in 0..geom.len() {
            let (meet, rows) = class_meet(&cx, c, ci, cyl, &mut footprints)?;
            tangencies.extend(rows);
            row.push(meet);
        }
        meets.push(row);
    }
    // The classes whose answer about cylinder `ci` is `want`, and the tangent ones carrying a face
    // of its other operand (`ClassMeet::Tangent`'s `other_owner`), whether or not its faces
    // wrote rows.
    let meeting = |ci: usize, want: ClassMeet| {
        meets[ci]
            .iter()
            .enumerate()
            .filter(move |(_, meet)| **meet == want)
            .map(|(c, _)| c)
    };
    let tangent_walls_of = |ci: usize| {
        meets[ci]
            .iter()
            .enumerate()
            .filter(|(_, meet)| {
                matches!(
                    meet,
                    ClassMeet::Tangent {
                        other_owner: true,
                        ..
                    }
                )
            })
            .map(|(c, _)| c)
    };

    // **A sectioned class shares a line of a lateral face with another record of that cylinder only
    // where the gate states it.** `shared_rulings` writes the pairs whose two planes have world
    // coefficients; this class has none, so a line it meets another recorded class in, lying on a
    // lateral face, would be a ruling no record states. Every such pair is shown off the faces
    // ([`line_off_lateral`]), or — for a cutting class — shown to be this class's hinge
    // ([`hinge_ruling`]) and recorded here, or refused. The hinge is the one world-rational line a
    // class with no world statement can be shown to hold exactly: a turn off the quarters fixes
    // only its own axis line, and an interval proves two points apart, never one.
    // Two walls whose corner line lies on a fillet's quarter are refused here although their faces
    // stand apart; past the gate the split refuses the same shape by the same name.
    //
    // An apart class's line with another record is not asked off the lateral: its arrangement
    // draws none of the cylinder's rulings ([`WorkingCyl::draws_rulings`]), and that includes a
    // sectioned class of the same cylinder. An apart plane that does not run along the axis meets
    // any other record's plane along the axis in a line across it, which is no ruling. One that
    // runs along the axis is apart only where neither the intervals nor the exact reading before
    // its turn part it from the radius (`class_meet`); such a plane has population 0 in the
    // suite, the census, the ignored sweep and the kit, so that its line with another record is
    // safe on a lateral face is a risk, not a proof.
    let mut hinged: Vec<(usize, SharedRuling)> = Vec::new();
    let mut grazed: Vec<(usize, GrazeLine)> = Vec::new();
    for (ci, cyl) in cyls.iter().enumerate() {
        let mut pieces: Option<Vec<Footprint>> = None;
        for c in meeting(ci, ClassMeet::Silent(Silence::Sectioned)) {
            let fps: &[Footprint] = pieces.get_or_insert_with(|| lateral_pieces(faces, cyl.surf));
            let others = meeting(ci, ClassMeet::Cuts)
                .chain(tangencies.iter().filter(|t| t.cyl == ci).map(|t| t.wall))
                .chain(tangent_walls_of(ci))
                .chain(meeting(ci, ClassMeet::Silent(Silence::Sectioned)).filter(move |&x| x != c));
            for x in others {
                if line_off_lateral(jd, c, x, &cyl.def, fps) {
                    continue;
                }
                // A cutting class states the line as a chart station; the sectioned class draws
                // it on its own arrangement only, so it is not one of `stated_by`.
                let at = (meets[ci][x] == ClassMeet::Cuts)
                    .then(|| hinge_ruling(model, jd, c, x, &cyl.def))
                    .flatten()
                    .filter(|at| ruling_on_a_face(&cyl.def, fps, at) == Some(true))
                    .ok_or_else(undecided)?;
                hinged.push((
                    ci,
                    SharedRuling {
                        line: (c.min(x), c.max(x)),
                        stated_by: vec![x],
                        at,
                        tangent_halves: None,
                    },
                ));
            }
        }
        // A grazing class answers the same way: its tangent line is a ruling of the cylinder, so
        // another record's line with it on a lateral face is shown off the faces, or — for a
        // sectioned class, which draws that ruling — shown to be it before their common turn
        // ([`graze_ruling`]) and recorded, or refused. A named plane can hold the tangent line too
        // (a turn of 45° puts its points on `x − y = c`), so every record is asked.
        // A named partner has no record to take: turned back by an angle off the quarters, a named
        // wall is not rational, so no proof before the turn stands, and the pair is refused. That
        // refusal is wider than the defect it guards — the defect is a wall face reaching the line,
        // which makes one point with two names at a ruling's end, while the question is whether the
        // line lies on a lateral face. Measured: an `x = y` prism across a half cylinder, with the
        // wall of a box turned 45° tangent at that ruling, puts no face on the line, and with the
        // refusal off it builds in both orders by derived volume. The narrowed question cannot be
        // measured yet: a shape whose face does run through the line is refused first, at a corner
        // where a named face calls a nameless class (`corners.rs`).
        for g in meeting(ci, ClassMeet::Silent(Silence::Grazes)) {
            let fps: &[Footprint] = pieces.get_or_insert_with(|| lateral_pieces(faces, cyl.surf));
            let others = meeting(ci, ClassMeet::Cuts)
                .chain(tangencies.iter().filter(|t| t.cyl == ci).map(|t| t.wall))
                .chain(tangent_walls_of(ci))
                .chain(meeting(ci, ClassMeet::Silent(Silence::Sectioned)))
                .chain(meeting(ci, ClassMeet::Silent(Silence::Grazes)).filter(move |&x| x != g));
            for x in others {
                if line_off_lateral(jd, g, x, &cyl.def, fps) {
                    continue;
                }
                if meets[ci][x] != ClassMeet::Silent(Silence::Sectioned)
                    || !graze_ruling(model, jd, x, g, &cyl.def)
                {
                    return Err(undecided());
                }
                // The fold reads the faces' ends on the line: a face whose ends are not heights (a
                // slanted rim) leaves the line's extent unstated, which refuses rather than folding
                // the names past where the face stops.
                let spans = lateral_stretches(faces, cyl.surf).ok_or_else(undecided)?;
                grazed.push((
                    ci,
                    GrazeLine {
                        section: x,
                        graze: g,
                        spans,
                    },
                ));
            }
        }
    }
    for (cyl, row) in cyls.iter_mut().zip(&meets) {
        cyl.records = row
            .iter()
            .enumerate()
            .filter_map(|(c, meet)| match meet {
                ClassMeet::Silent(s) => Some((c, CylRecord::Silent(*s))),
                ClassMeet::Cuts => Some((c, CylRecord::Cuts)),
                ClassMeet::Sections => Some((c, CylRecord::Sections)),
                ClassMeet::Tangent { .. } => Some((c, CylRecord::Touches)),
                ClassMeet::Across | ClassMeet::Clear | ClassMeet::Undrawable => None,
            })
            .collect();
    }
    debug_assert!(
        cyls.iter().all(|cyl| (0..geom.len())
            .filter(|&c| geom[c].world.is_none())
            .all(|c| cyl.is_silent(c))),
        "every nameless class is silent about every cylinder, or the gate refused"
    );

    // An oblique section no road draws (`ClassMeet::Undrawable`).
    if meets
        .iter()
        .flatten()
        .any(|meet| *meet == ClassMeet::Undrawable)
    {
        return Err(reject(RejectReason::ObliqueCylinderCut));
    }
    // **A class's sections of two cylinders nest only across parallel axes** ([`carries_unnested`]).
    let mixed = (0..geom.len()).any(|c| {
        let row: Vec<([nacre_exact::Rat; 3], &ClassMeet)> = cyls
            .iter()
            .zip(&meets)
            .map(|(cy, meets)| (cy.def.dir(), &meets[c]))
            .collect();
        carries_unnested(&row)
    });
    if mixed {
        return Err(reject(RejectReason::ObliqueCylinderCut));
    }

    // **The pair rule asks classes of different owners whether their faces share a point, and
    // writes what it cannot prove.** The proposition is "the two classes share no face", decided in
    // order of cost. Two classes of one solid share none by construction — a valid solid's faces
    // meet only along their edges — so those pairs are not asked (the fillets of one plate, the two
    // half cylinders of one slot). The axis distance exceeding the radius **sum** decides most
    // inputs (a drill crossing a bore with room to spare). It is not the whole rule — that would
    // make a fact about surfaces read as a fact about faces: a stud fused through a cube and then a
    // second stud across it would be refused because their *axes* cross, though the first stud's
    // remaining faces sit past `|z| = 0.5` and the second's whole surface within `|z| = 0.2`. So a
    // pair the distance cannot clear asks the faces themselves — the plane–cylinder arm's per-face
    // question (`face_clears_footprint`), spelled for a lateral face's reach against the other's
    // along every direction the pair can state (`separating_dirs`: for skew axes the common
    // perpendicular is where a fillet beside a crosswise drill is seen apart), or one face lying
    // inside the other's cylinder, which is where a pin nested crosswise in a larger cylinder is
    // seen apart (`pieces_clear`). For parallel axes the nested test (the arm below) comes first,
    // and the reach along a parallel axis has no radial term, so what is left is the axis spans,
    // the common cross-section and containment.
    let mut cyl_pairs: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for (i, a) in cyls.iter().enumerate() {
        for (j, b) in cyls.iter().enumerate().skip(i + 1) {
            if matches!(
                (a.owner, b.owner),
                (SolidSide::A, SolidSide::A) | (SolidSide::B, SolidSide::B)
            ) {
                continue;
            }
            let clear = match nacre_exact::cylinders_clear(
                &a.def.origin(),
                &a.def.dir(),
                a.def.r2(),
                &b.def.origin(),
                &b.def.dir(),
                b.def.r2(),
            ) {
                Orient::Positive => true,
                _ if same_surface(&a.def, &b.def) => false,
                // Parallel axes with one infinite cylinder strictly inside the other — a pin in a
                // bore, a boss under a smaller pin — never meet as surfaces either: the second
                // sufficient condition the distance rule has for parallel pairs.
                _ if nacre_exact::parallel_rat(&a.def.dir(), &b.def.dir())
                    && nacre_exact::cylinders_nested(
                        &a.def.origin(),
                        &a.def.dir(),
                        a.def.r2(),
                        &b.def.origin(),
                        b.def.r2(),
                    ) == Orient::Negative =>
                {
                    true
                }
                _ => lateral_faces_clear(faces, a, b),
            };
            if !clear {
                cyl_pairs.insert((i, j));
            }
        }
    }
    // A recorded pair is refused (`CylinderPairContact`).
    //
    // **An obligation inherited by whatever opens this refusal.** The arrangement leans on this
    // line for a proposition of its own: a class may carry a circle (⊥ one cylinder) and rulings (∥
    // another) at once, and that is safe **only while the two never meet** — neither split cuts the
    // other's kind, so a crossing they both walked past would be a node no road mints, and the
    // point is `plane ∩ cylinder ∩ cylinder`, a nested radical with no name. Two such edges meeting
    // means two lateral faces share a point, which is exactly what this refusal denies.
    // `split_class_edges` carries the derivation and a shipped check for the circle × ruling pair, the
    // one of the three that has a cheap exact predicate; when this population is admitted, the
    // other two pairs need theirs.
    if !cyl_pairs.is_empty() {
        return Err(reject(RejectReason::CylinderPairContact));
    }
    for (ci, sr) in shared_rulings(faces, jd, &cyls, &meets, &tangencies)
        .ok_or_else(|| reject(RejectReason::WitnessNotRational))?
        .into_iter()
        .chain(hinged)
    {
        // Both classes of a shared line hold it as a station of their own arrangements — the
        // readers that fold a class's segments onto the ruling ask only whether the line is on
        // the class (`split_rulings`), so the record is where that is true.
        debug_assert!(
            cyls[ci].holds_ruling(sr.line.0) && cyls[ci].holds_ruling(sr.line.1),
            "a shared line pairs two classes that hold a ruling of the cylinder"
        );
        cyls[ci].shared.push(sr);
    }
    for (ci, gl) in grazed {
        cyls[ci].grazed.push(gl);
    }
    // **A row that could not be stated does not pass** — refused here, by the row's own name,
    // rather than arranged: the verdict would have nothing to read, and where both solids run
    // through the line it is the only reader. So is a row whose stretches through the line could
    // not be read: read as running through, a solid that ends on the line would be judged by the
    // three regions it does not have; read as ending there, a pinch nothing else sees would pass.
    for t in &tangencies {
        if let Some(e) = t.undecided {
            return Err(e);
        }
        if t.through.is_none() {
            return Err(crate::reject_at(
                RejectReason::CylinderGateUndecided,
                crate::RejectWhere::Point(t.witness),
            ));
        }
    }
    Ok((cyls, tangencies))
}

/// **What a plane class receives from a cylinder** — [`class_meet`]'s answer: the curve the
/// class's arrangement meets on the cylinder, or what proves it meets none. Read only by the gate:
/// its refusals (an ellipse, the sectioned settlement) and the records it hands on.
#[derive(Debug, PartialEq, Eq)]
enum ClassMeet {
    /// The plane is perpendicular to the axis: it meets the lateral in circles, which the faces'
    /// rims carry.
    Across,
    /// No face of the class reaches a lateral face: the plane clears the radius or misses every
    /// lateral face, or the class's faces clear the strip with no stated footprint to prove a cut,
    /// or (oblique) every lateral face misses the plane.
    Clear,
    /// The plane runs within the radius and cuts a lateral face: two rulings, drawn on the class.
    Cuts,
    /// The plane stands exactly `r` from the axis: one grazing line. `holders` are the classes
    /// holding that line ([`classes_holding_the_line`]); `other_owner` says the class carries a
    /// face of the cylinder's other operand — such a tangent class is half of [`SharedRuling`]'s
    /// population, and a record the sectioned settlement weighs.
    Tangent {
        holders: Vec<usize>,
        other_owner: bool,
    },
    /// A class with no world statement, silent about the cylinder. `Sectioned` is a candidate
    /// until [`cylinder_gate`] settles it.
    Silent(Silence),
    /// An oblique plane touching none of the circles the lateral faces' boundaries ride: the class
    /// draws the section — a whole ellipse or arcs of one — on its own arrangement
    /// ([`CylRecord::Sections`]).
    Sections,
    /// A section the class cannot draw: an oblique plane not shown to miss the lateral faces that
    /// touches one of their boundary sections at a point, or whose faces' boundaries cannot be
    /// read; or a plane across the axis touching a slanted boundary section at a point.
    Undrawable,
}

/// What [`class_meet`] reads besides the pair: the operands' face table and the judge, and which
/// cylinders leave only whole rims.
#[derive(Clone, Copy)]
struct PairCx<'a, 'j> {
    model: &'a Model,
    faces: &'a [FaceRow],
    plane_ix: &'a [ClassIx],
    n_a: usize,
    jd: &'a Judge<'j, WorkingPlane>,
    whole_rims: &'a [bool],
}

/// **What plane class `c` receives from cylinder `cyl`** — the gate's one answer per pair.
///
/// With `n` the class's rational normal and `m`/`o`/`r` the cylinder's raw axis/origin/radius:
/// - a class with no world statement is silent about the cylinder or refused here
///   ([`Silence`], asked beside, grazing, sectioned, apart in that order; `Sectioned` is a
///   candidate [`cylinder_gate`] settles once every pair is asked); a `Wide` name is refused for
///   its width;
/// - `n × m = 0` — a perpendicular cut, [`ClassMeet::Across`], a cap seated flush on the other
///   body included;
/// - `n · m = 0` — a wall parallel to the axis. It must provably miss the **rectangle** the
///   cylinder occupies in that plane. The infinite plane clearing the axis by more than `r`
///   (`(n·o + d)² > r²·|n|²`) settles it outright and decides most inputs; so does the plane
///   missing every lateral **face**, read on the faces' arcs ([`laterals_clear_plane`] — a plane
///   that cuts the infinite cylinder only off a fillet's quarter). Otherwise each face on the
///   class answers for itself, across the strip or along a lateral face's span
///   ([`face_clears_footprint`]). A face not shown to miss either rides the rulings road — the
///   wall's plane within the radius (`0 ≤ d < r`, through the axis or offset from it),
///   [`ClassMeet::Cuts`] — or, the plane exactly `r` from the axis, is [`ClassMeet::Tangent`]
///   with its [`Tangency`] rows.
///
///   **One clearance call, three answers.** Past the lateral faces' question, `Positive` clears;
///   `Negative` cuts (two rulings); `Zero` is tangent (one grazing line) — not `Cuts`, for the
///   reason the record-and-pass arm gives. A row that could not be stated is refused at the gate's
///   end, by its own name, before anything is arranged.
/// - anything else — an oblique plane, [`ClassMeet::Clear`] when every lateral face of the
///   cylinder provably misses the plane — the face's reach along `n` against the plane's station
///   ([`lateral_reach`], [`laterals_clear_plane`]) — then [`ClassMeet::Sections`] when the plane
///   touches none of the circles the lateral faces' boundaries ride ([`oblique_section_meet`]),
///   and otherwise [`ClassMeet::Undrawable`] — questions about the lateral faces, their reach and
///   their circles, so a slanted plane is not refused for cutting the infinite cylinder alone
///   ([`RejectReason::ObliqueCylinderCut`] lists what passes).
///
/// The order of the questions is the order of the refusals they can raise, and so part of what
/// this function answers: a `Zero` pair asks the lateral faces before the classes holding its
/// line, and those before its faces.
fn class_meet(
    cx: &PairCx<'_, '_>,
    c: usize,
    ci: usize,
    cyl: &WorkingCyl,
    footprints: &mut Option<Vec<Footprint>>,
) -> Result<(ClassMeet, Vec<Tangency>), BoolError> {
    let PairCx {
        model,
        faces,
        plane_ix,
        n_a,
        jd,
        whole_rims,
    } = *cx;
    use nacre_exact::Orient;
    let undecided = || reject(RejectReason::CylinderGateUndecided);
    let (o, m, r2) = (cyl.def.origin(), cyl.def.dir(), cyl.def.r2());
    let wp = &jd.planes[c];
    // **A world description, or the proof that none is needed** — the question this function
    // asks is geometric (does this class's plane clear that cylinder), and both sides have to
    // speak about the world. `rotated` is not that question: a plane whose truth carries a
    // translation is *judged* through its chain and still has exact world coefficients, and that
    // is the population the rulings road serves.
    let Some(coeffs) = jd.world_rat(c) else {
        if wp.world.is_some() {
            return Err(reject(RejectReason::WitnessNotRational));
        }
        let stance = nameless_stance(model, jd, c, cyl);
        if let Stance::Beside(side) = stance {
            return Ok((ClassMeet::Silent(Silence::Beside(side)), Vec::new()));
        }
        let apart = faces_apart(model, faces, plane_ix, c, &cyl.boxes);
        // A tangent plane divides nothing, so the rims are not asked — before the filled reading,
        // which asks them.
        if stance == Stance::Grazes && apart {
            return Ok((ClassMeet::Silent(Silence::Grazes), Vec::new()));
        }
        // With every face of the class apart from the lateral boxes, a plane cutting the
        // lateral in rulings draws them on its own arrangement only ([`Silence::Sectioned`]),
        // and a disk of the cylinder on the class its chord — the section drawn truly, whatever
        // the rims. Asked before the filled reading: a class that read one cylinder filled and
        // another, inside it, truly would meet itself across the two readings (`LabelConflict`).
        if stance == Stance::Rulings && apart {
            return Ok((ClassMeet::Silent(Silence::Sectioned), Vec::new()));
        }
        if whole_rims[ci] && apart {
            return Ok((ClassMeet::Silent(Silence::Apart), Vec::new()));
        }
        return Err(undecided());
    };
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // A perpendicular cut is the circle population, and it passes — **including a cap seated
    // flush on the other body's face**. Refusing that seating would be wider than anything the
    // refusal could name: what makes a seated circle hard is its boundary meeting the
    // counterpart's boundary, and a boundary is either an edge on a plane (that plane is parallel
    // to the axis → the wall rule below, or oblique → the oblique rule; both are conservative: a
    // wall plane that cuts a lateral face is recorded whatever its faces answer, and an oblique
    // plane touching a rim is refused — `oblique_section_meet` asks the rim's whole circle) or another
    // cylinder's rim (two circles can only overlap when the axes stand closer than r₁+r₂ → the
    // pair rule below). So what a seated rule would turn away is exactly the population whose
    // circle lies wholly inside the counterpart's face — measured across the family (through hole,
    // blind hole, boss fuse, common, drilling a pocket floor): every one exact and validating
    // clean. Those rules **name** a tangency or a crossing before the arrangement sees it, so no
    // degenerate seating reaches the arrangement unnamed: a touch passes with a `Tangency` row,
    // refused here only where its verdict could not judge the row.
    // The one relation rule ([`nacre_exact::axis_relation`]) — the same one that decides an
    // edge's curve (`Model::derive_edge_curve`), so the gate and the edges cannot part.
    let relation = nacre_exact::axis_relation(&n, &m);
    if relation != nacre_exact::AxisRelation::Across {
        if relation == nacre_exact::AxisRelation::Oblique {
            // **The oblique arm asks the faces**: if every lateral face provably misses the
            // infinite plane (`laterals_clear_plane`), no face of that plane's class can meet
            // the lateral, and the caps are the plane–plane arrangement's business. A pair not
            // shown to miss meets the cylinder in an ellipse: the class's arrangement sections the
            // cylinder's own solid by the whole plane, and that section runs along the ellipse —
            // whole between a band's rims, or arcs of it ending where a rim, a ruling or a
            // segment of the class crosses it — which the class draws on its own arrangement and
            // the lateral's chart reads as a line, wherever the plane touches none of the lateral's
            // circles ([`ClassMeet::Sections`], [`oblique_section_meet`]).
            let fps = footprints.get_or_insert_with(|| lateral_footprints(faces, cyl.surf));
            let meet = if laterals_clear_plane(&cyl.def, fps, &coeffs) {
                ClassMeet::Clear
            } else {
                oblique_section_meet(model, faces, cyl, fps, &coeffs)
                    .unwrap_or(ClassMeet::Undrawable)
            };
            return Ok((meet, Vec::new()));
        }
        // A parallel wall must provably miss the lateral surface. **The question is about the
        // wall's *faces*, not its plane** — the uniform-slab theorem this feeds says so in its own
        // words ("the other operand's **boundary** does not meet the open cylinder slab"). The two
        // cheaper *sufficient* conditions (this function's doc) pass the pair with nothing
        // recorded: the class's arrangement gets nothing from this cylinder, whose rim arcs bound
        // its faces and so miss the plane too. Refusing on the plane alone would turn away a whole
        // family the engine serves — a boss standing far away whose wall plane, extended, happens
        // to pass through a hole.
        //
        // A face is asked whether it misses the **rectangle** this cylinder occupies in that
        // plane (`face_clears_footprint`). Its strip is the whole circle's, and past the second
        // question that decides no record wrongly: a plane that cuts a lateral face is recorded
        // below whatever its faces answer, a tangent plane that reaches a face's arc has that
        // tangent line for its strip, and a cylinder with no usable footprint has no arc to read.
        let clearance = nacre_exact::point_plane_clearance_rat(&coeffs, &o, r2);
        if clearance != Orient::Positive
            && laterals_clear_plane(
                &cyl.def,
                footprints.get_or_insert_with(|| lateral_footprints(faces, cyl.surf)),
                &coeffs,
            )
        {
            return Ok((ClassMeet::Clear, Vec::new()));
        }
        // **The classes holding a tangent line, asked once for the pair** — `SharedRuling` reads
        // them. It is asked here, before any wall face is: the rows are written only for faces that fail to
        // clear the footprint and touch the line within a lateral's span, so a face standing on a
        // cap, or one the line passes through a hole of, writes none — yet its line is the same
        // line.
        let mut rows = Vec::new();
        let holders = if clearance == Orient::Zero {
            let Some(h) = classes_holding_the_line(jd, c, &coeffs, &cyl.def) else {
                return Err(reject(RejectReason::WitnessNotRational));
            };
            let other_owner = faces.iter().enumerate().any(|(i, row)| {
                matches!(row, FaceRow::Plane(_))
                    && plane_ix[i] == ClassIx::Plane(c)
                    && (i < n_a) != (cyl.owner == SolidSide::A)
            });
            Some((h, other_owner))
        } else {
            None
        };
        if clearance != Orient::Positive {
            // The face-level test reads each face's own vertices, which are realized
            // world coordinates — so it needs no frame guard of its own; `coeffs` above is
            // already the world description (a class without one never reaches here). A
            // frame guard here would refuse every moved class outright and keep a
            // translated body out of the cylinder roads.
            let spans: Vec<[nacre_exact::Rat; 2]> = footprints
                .get_or_insert_with(|| lateral_footprints(faces, cyl.surf))
                .iter()
                .map(|f| f.span.expect("a listed footprint has a span"))
                .collect();
            if !wall_faces_clear(model, faces, plane_ix, c, &coeffs, &o, &m, r2, &spans)? {
                // **The record-and-pass arm**: a wall whose plane runs **within** the radius —
                // any `0 ≤ d < r`, the through-axis wall included — and whose faces did not clear
                // the strip is recorded and passed; the tracer's ruling/chord arms fire only on
                // recorded pairs and read the *faces* (a ruling piece comes from a face's own
                // cycles, a cap's chord is the true meet of that class's plane with the disk), so
                // a recorded face that misses the lateral contributes nothing. The record is «may
                // meet; the tracer decides».
                //
                // **No offset is refused** (`0 < d < r`): bosses and bores, rational and
                // irrational rulings, axis inside or outside the plate, a split bore — every one
                // assembles, validates clean and answers the exact volume oracle. The arithmetic
                // (`plane_plane_cylinder`'s roots, `ruling_side`, the chart's circular order)
                // assumes no diameter. And a contact is a question a fence cannot answer:
                // `validate` cannot see it (`Cut` returns `Ok`, no net sees the line), which is
                // why the answer is a judge.
                //
                // **A tangency is a graze, not a crossing** — it divides nothing, so it passes,
                // and it is **not** a cut. A cut's proposition is "the plane runs *within* the
                // radius", which a tangent plane does not; the three roads gated on that record
                // assert it in a `debug_assert` (`arrangement`' circular-hole arm,
                // `chord_stations`, `rulings_on_class`), and all three stay true because this pair
                // never reaches them — a tangency is a double root where those roads spell two
                // distinct ones, so a cut's record would put the pair on the ruling road and
                // assemble nothing. Kept apart, the arrangement needs no tangency arm: of 21
                // measured cells, the 15 outside the third-plane population (a third class's face
                // reaching the line) assemble, `validate` clean and
                // volumes exact. Its record is its own (`CylRecord::Touches`): the line divides
                // no cell of this plane, nor the chart where the wall's face runs through it, and is
                // the plane's one station on the cylinder. A face of the wall **ending** on the line
                // is an edge there, which the chart's station on a secant wall's ruling reads off the
                // pair's record (`SharedRuling::tangent_halves`).
                //
                // What the pair can still do is pinch the *result*, and that is a question about
                // the operation — so the geometry is stated in a row and
                // `assembly::tangency_reject` asks `keep` — beside `self_touch_reject`, where the
                // grouping can also say whether the pieces share a solid. A row that could not be
                // stated, or whose stretches through the line could not be read, is refused at this
                // gate's end.
                if clearance == Orient::Zero {
                    rows = tangency_rows(
                        model, faces, plane_ix, n_a, c, ci, &coeffs, &cyl.def, cyl.surf, cyl.owner,
                        jd,
                    );
                } else {
                    return Ok((ClassMeet::Cuts, Vec::new()));
                }
            } else if clearance == Orient::Negative
                // Only a footprint that can be **stated** proves a cut — the plane did not
                // clear the faces above, so a stated one is cut; an unstatable one (a
                // rotated class, no span) proves nothing either way and records nothing —
                // recording on it draws rulings no face has (measured, the rigid-motion
                // oracle).
                && !footprints
                    .get_or_insert_with(|| lateral_footprints(faces, cyl.surf))
                    .is_empty()
            {
                // **The class cuts this lateral face though no face of the other solid touches
                // it**: the refusal question is face against face and it is clear, but the
                // arrangement on this class still needs the ruling — every face of *this* solid
                // that crosses the class ends on it (a cap's section ends where its rim arc meets
                // the plane), and without the ruling that end dangles and the walk doubles back
                // (measured: a gusset beside a filleted plate, its side plane within the fillet's
                // radius). So the record is «the plane cuts the face», not «the faces touch».
                return Ok((ClassMeet::Cuts, Vec::new()));
            }
        }
        return Ok(match holders {
            Some((holders, other_owner)) => (
                ClassMeet::Tangent {
                    holders,
                    other_owner,
                },
                rows,
            ),
            None => (ClassMeet::Clear, Vec::new()),
        });
    }
    // A plane across the axis touching a slanted rim at a point is the touch the oblique arm
    // refuses: one point of the two sections, no crossing to name it by.
    if touches_a_boundary(model, faces, cyl, &coeffs) {
        return Ok((ClassMeet::Undrawable, Vec::new()));
    }
    Ok((ClassMeet::Across, Vec::new()))
}

/// **Two plane classes along a cylinder's axis whose meet line lies on the cylinder** — a ruling
/// of one solid that is also the line of the other's edge (or of the same solid's flat through
/// its own axis). Recorded once, here, and read by every road that must say so: the alias seed
/// (one point, several names), the rulings split (one line, two vocabularies), the chart (one
/// station, its secant walls in chord order), the result's edge subdivision (the pierce points on
/// the line), the straight-angle cleaning (a ruling piece and a plane segment, one line) and the
/// self-touch verdict (an edge along the line inside the solid's own lateral).
///
/// The population is the pairs [`class_meet`] already answered against the cylinder: a class that
/// cuts it ([`ClassMeet::Cuts`] — two rulings) with another that cuts it, or with a **tangent**
/// class that carries a face of the cylinder's other operand ([`ClassMeet::Tangent`]'s
/// `other_owner`). A tangent class of the cylinder's own solid is its smooth edge — a fillet — and
/// stays out: a valid solid's faces do not meet. And a class with no world statement that
/// sections the cylinder ([`Silence::Sectioned`]) with a cutting class, where the line is the
/// sectioned class's hinge ([`hinge_ruling`]) — written by the sectioned settlement in
/// [`cylinder_gate`], since [`shared_rulings`] reads world coefficients.
///
/// **And the line must lie on a lateral *face*, not only on the surface** — like the gate's
/// questions, this record speaks about faces. A fillet is a quarter of its cylinder: a slab
/// standing on the plate beside it, with a face through the fillet's axis and another tangent to
/// the cylinder at the quarter turn the fillet does not have, puts its corner line on the infinite
/// cylinder and nowhere near the lateral face
/// (`a_plane_through_the_fillet_axis_shares_the_tangent_ruling`: measured, without this condition
/// the pair would be recorded). The line is asked of each lateral face ([`line_apart`] — on the arc
/// or on an end is on the face: a line on a face's end ruling is the face's edge, the
/// half-cylinder family) — which a [`Tangency`] row asks too, where the face's own stretches of the
/// line cannot be read.
#[derive(Clone, Debug)]
pub(crate) struct SharedRuling {
    /// The two plane classes, ascending.
    pub(crate) line: (usize, usize),
    /// Which of the two state the ruling as a station of the cylinder's chart — the ones that cut
    /// the cylinder, so one (a tangent or a sectioned partner) or both. Both classes of the pair fold
    /// their own segments on the line onto it (`split_rulings`); the chart gathers every secant
    /// wall through the line into one station (`cyl_chart::stations_of_lines`).
    pub(crate) stated_by: Vec<usize>,
    /// A point of the line, exact: the line is the ruling through it, so a reader needs neither
    /// class's coefficients to say where it is — a sectioned partner has none.
    pub(crate) at: [nacre_exact::Rat; 3],
    /// **For a tangent partner, where its faces lie on the line's `+τ` and `−τ` sides** — the
    /// [`Tangency`] rows' halves of that class and cylinder, gathered (closed pieces in the axis
    /// parameter, joined). One side alone is a face ending on the line, which the line's station
    /// reads as a spoke; both sides is the face running through it — a contact, which divides
    /// nothing. `None` for a pair of secant or sectioned partners, and for a tangent one whose rows
    /// could not be read, which the gate refuses by the rows' own name before anything reads this.
    pub(crate) tangent_halves: Option<[Vec<[nacre_exact::QuadVal; 2]>; 2]>,
}

impl SharedRuling {
    /// The vector from the axis of `def` out to the line, square to the axis — the line's
    /// direction across the cylinder, read off [`SharedRuling::at`]. `None` is overflow.
    pub(crate) fn radial(&self, def: &nacre_topo::CylinderDef) -> Option<[nacre_exact::Rat; 3]> {
        let o = def.origin();
        let mut w = [nacre_exact::Rat::from_int(0); 3];
        for k in 0..3 {
            w[k] = self.at[k].checked_sub(o[k])?;
        }
        nacre_exact::off_axis_rat(&w, &def.dir())
    }

    /// **Whether a point on plane `wx` lies on the line** — the point `line + s·dir` against the
    /// plane through the line with normal `n_x × m`. A plane holding a ruling runs along the axis
    /// `m`, so that normal is not zero, and the two planes meet in the line alone: a point of `wx`
    /// is on it exactly when it is on this plane — and that is so of a point of `wx`'s other
    /// ruling, or of another cylinder, alike. `None` is overflow.
    pub(crate) fn holds(
        &self,
        wx: &[nacre_exact::Rat; 4],
        m: &[nacre_exact::Rat; 3],
        line: &nacre_exact::quad::MeetLine,
        s: &nacre_exact::quad::QuadVal,
    ) -> Option<bool> {
        let n = nacre_exact::cross3_rat(&[wx[0], wx[1], wx[2]], m)?;
        let across = nacre_exact::plane_from_point_normal(n, self.at)?;
        Some(nacre_exact::quad::plane_side(&across, line, s) == nacre_exact::Orient::Zero)
    }
}

/// **A sectioned class and a grazing class whose meet line is a ruling of the cylinder** — the
/// ruling the sectioned class draws is the line the grazing class touches along. Every class across
/// the axis meets that line in one point with two names — the three-plane one and a root of its
/// pair with the sectioned class — and intervals never prove two names one point; the gate proves
/// it before the two classes' common turn ([`graze_ruling`]) and records it here for the name table
/// (`Aliases::record_graze_line`). [`SharedRuling`] cannot carry it: its point of the line is world
/// rational, and a ruling turned off the quarters about the axis has none.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct GrazeLine {
    /// The class that sections the cylinder ([`Silence::Sectioned`]).
    pub(crate) section: usize,
    /// The class touching the cylinder along the line ([`Silence::Grazes`]).
    pub(crate) graze: usize,
    /// The lateral faces' axis spans: a class across the axis meets the line on a face only where
    /// its station falls within one.
    pub(crate) spans: Vec<[nacre_exact::Rat; 2]>,
}

/// The plane classes other than `wall` that hold `wall`'s tangent line to the cylinder `def` —
/// `plane_plane_cylinder` answering `OnRuling` (or the two classes being one plane). The one
/// spelling of «this line lies on the cylinder» the record below uses for every pair. `None` is
/// arithmetic that could not answer.
fn classes_holding_the_line(
    jd: &Judge<'_, WorkingPlane>,
    wall: usize,
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
) -> Option<Vec<usize>> {
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let mut out = Vec::new();
    for k in 0..jd.planes.len() {
        if k == wall {
            continue;
        }
        let Some(w) = jd.world_rat(k) else { continue };
        if matches!(
            nacre_exact::quad::plane_plane_cylinder(coeffs, &w, &o, &m, r2)?,
            nacre_exact::quad::CylinderMeet::OnRuling(_)
                | nacre_exact::quad::CylinderMeet::CoincidentPlanes
        ) {
            out.push(k);
        }
    }
    Some(out)
}

/// Build the [`SharedRuling`] record from the gate's own records — see its doc for the
/// population. `None` is arithmetic that could not answer.
fn shared_rulings(
    faces: &[FaceRow],
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    meets: &[Vec<ClassMeet>],
    tangencies: &[Tangency],
) -> Option<Vec<(usize, SharedRuling)>> {
    let mut out = Vec::new();
    for (ci, cyl) in cyls.iter().enumerate() {
        // In class order, so sorted.
        let secants: Vec<usize> = meets[ci]
            .iter()
            .enumerate()
            .filter(|(_, meet)| **meet == ClassMeet::Cuts)
            .map(|(c, _)| c)
            .collect();
        let (o, m, r2) = (cyl.def.origin(), cyl.def.dir(), cyl.def.r2());
        let fps = lateral_pieces(faces, cyl.surf);
        let on_a_face = |p: &[nacre_exact::Rat; 3]| ruling_on_a_face(&cyl.def, &fps, p);
        for (i, &s) in secants.iter().enumerate() {
            let ws = jd.world_rat(s)?;
            for &t in &secants[i + 1..] {
                let wt = jd.world_rat(t)?;
                if let nacre_exact::quad::CylinderMeet::OnRuling(line) =
                    nacre_exact::quad::plane_plane_cylinder(&ws, &wt, &o, &m, r2)?
                    && on_a_face(&line.base())?
                {
                    out.push((
                        ci,
                        SharedRuling {
                            line: (s, t),
                            stated_by: vec![s, t],
                            at: line.base(),
                            tangent_halves: None,
                        },
                    ));
                }
            }
        }
        for (t, meet) in meets[ci].iter().enumerate() {
            let ClassMeet::Tangent {
                holders,
                other_owner: true,
            } = meet
            else {
                continue;
            };
            let foot = tangency_foot(&jd.world_rat(t)?, &o)?;
            if !on_a_face(&foot)? {
                continue;
            }
            // Every row of this class on this cylinder — one per wall face and lateral face that
            // touch — reads the same faces of the wall's solid, over its lateral's stretches.
            let halves = tangencies
                .iter()
                .filter(|row| row.wall == t && row.cyl == ci)
                .try_fold([Vec::new(), Vec::new()], |mut acc, row| {
                    let h = row.halves.as_ref()?;
                    for (to, from) in acc.iter_mut().zip(h) {
                        to.extend(from.iter().copied());
                    }
                    Some(acc)
                })
                .and_then(|[plus, minus]| Some([joined(plus)?, joined(minus)?]));
            for &s in holders {
                if secants.binary_search(&s).is_ok() {
                    out.push((
                        ci,
                        SharedRuling {
                            line: (s.min(t), s.max(t)),
                            stated_by: vec![s],
                            at: foot,
                            tangent_halves: halves.clone(),
                        },
                    ));
                }
            }
        }
    }
    Some(out)
}

/// Whether the line through `p` along the axis of `def` lies on one of the lateral faces `fps`
/// ([`line_apart`]); `None` is arithmetic that could not say.
fn ruling_on_a_face(
    def: &nacre_topo::CylinderDef,
    fps: &[Footprint],
    p: &[nacre_exact::Rat; 3],
) -> Option<bool> {
    let line = Footprint {
        span: None,
        theta: Theta::at_line(p, &def.origin(), &def.dir())?,
    };
    for fp in fps {
        if !line_apart(def, &line, fp)? {
            return Some(true);
        }
    }
    Some(false)
}

/// **The ruling a sectioned class `c` shares with the cutting class `x`, as a proof** — `c`'s
/// plane holds the line its chain turns it about ([`Model::plane_hinge`]), and that line lies on
/// `x`'s plane and on the cylinder `def`: then `c ∩ x` is that ruling. `None` where any of it is
/// not shown.
fn hinge_ruling(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    x: usize,
    def: &nacre_topo::CylinderDef,
) -> Option<[nacre_exact::Rat; 3]> {
    use nacre_exact::Orient::Zero;
    let (p, d) = model.plane_hinge(jd.planes[c].surf)?;
    let w = jd.world_rat(x)?;
    let (o, m) = (def.origin(), def.dir());
    let along = nacre_exact::cross3_rat(&d, &m)?
        .iter()
        .all(|v| *v == nacre_exact::Rat::from_int(0));
    let on_x = nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], &d) == Zero
        && nacre_exact::plane_residual_sign(&nacre_exact::PlaneName::Narrow(w), p) == 0;
    let on_cyl = nacre_exact::quad::cylinder_radial_side(&p, &o, &m, def.r2()) == Zero;
    (along && on_x && on_cyl).then_some(p)
}

/// **Does every face on plane class `c` provably miss this cylinder?** — the boundary question
/// [`cylinder_gate`]'s wall rule asks once the cheaper plane-level one has failed.
///
/// `Err` is the honest "could not decide exactly"; `Ok(false)` means some face was not shown to
/// clear, which is the wall refusal's whole content.
///
/// **No table is built for this.** The scan runs only on the class that failed the plane test —
/// rare — so walking the face rows there costs nothing on the common path and allocates nothing.
/// A per-class face table computed for every boolean would be built for everyone and read by
/// almost no one. The `spans` this takes is held to the same rule: the caller builds it on first
/// use, not per cylinder — measured with it eager, an ordinary cut over a 16-bore plate builds it
/// 16 times and reads it **zero** times.
#[allow(clippy::too_many_arguments)]
fn wall_faces_clear(
    model: &Model,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    c: usize,
    coeffs: &[nacre_exact::Rat; 4],
    o: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    r2: &nacre_exact::BigRat,
    spans: &[[nacre_exact::Rat; 2]],
) -> Result<bool, BoolError> {
    let mut seen = 0usize;
    for (i, row) in faces.iter().enumerate() {
        let (ClassIx::Plane(k), FaceRow::Plane(fi)) = (plane_ix[i], row) else {
            continue;
        };
        if k != c {
            continue;
        }
        let fh = fi.face.expect("collect_planes yields real faces");
        seen += 1;
        if !face_clears_footprint(model, model.face(fh), coeffs, o, m, r2, spans)? {
            return Ok(false);
        }
    }
    // A class exists because faces made it, so finding none is the class table and the face rows
    // disagreeing rather than an input — and "every one of no faces clears" is a pass this must
    // not hand out by default.
    if seen == 0 {
        return Err(reject(RejectReason::CylinderStagesDisagree));
    }
    Ok(true)
}

/// **A wall plane exactly `r` from a cylinder's axis, and everything a verdict needs about it.**
///
/// **A tangency is a graze, not a crossing**: it earns no cut ([`CylRecord::Cuts`]), for the
/// reason [`class_meet`]'s record-and-pass arm gives. [`WorkingCyl::draws_rulings`] also takes a
/// class that sections the cylinder, `Silence::Sectioned`, whose premise is the same «within the
/// radius» inequality — read exactly before the class's turn where it turns about the axis, on the
/// interval road otherwise. A class with no world statement that touches the cylinder is no cut
/// either: it grazes (`Silence::Grazes`).
///
/// What such a plane can still do is **pinch the result**: near the tangent line the material
/// splits into three regions — the lens inside the cylinder, the **two** wedges between the
/// parabola and the plane, and the far half-space — and whether the kept ones hang together is a
/// question only [`crate::draft::keep`] can answer. So this is a *record*, not a verdict: the
/// gate states the geometry, the operation decides, and the gate refuses a row the verdict cannot
/// read ([`cylinder_gate`]'s last check).
///
/// **A row is about a lateral *face*, not its surface.** A half cylinder's surface is tangent to
/// a wall on the side the face does not have, where nothing touches: a row there made the verdict
/// convict a result that does not touch itself (census `arcprofile slot fuse`, a slot fused with a
/// box tangent to its half-circle end's circle on the missing side). So a row is written only for
/// a lateral face that lies on the line, and only for a wall face that meets the line **within
/// that lateral face's own stretches of it** ([`line_runs_through_face`] over each stretch
/// `lateral_cover_on_ruling` reads off the face's loops; where the face cannot be read there, its
/// angular extent stands in ([`line_apart`], the question [`SharedRuling`] asks), the wider
/// answer). A wall face that reaches the line only through a notch or a hole of its own, or only
/// where a window in the lateral leaves the line bare, or only above or below this lateral,
/// touches nothing here and writes no row, as a face that clears the footprint writes none: the
/// line is still recorded, class by class, for the shared rulings. So on a surface with two
/// lateral faces — the lower holding a half turn, the upper the other half — a wall touching the
/// lower's height on the upper's side writes no row for the upper (no built shape has two such
/// faces yet).
#[derive(Clone, Debug)]
pub(crate) struct Tangency {
    /// The plane class the wall face lies on.
    pub(crate) wall: usize,
    /// The cylinder class the lateral face lies on.
    pub(crate) cyl: usize,
    pub(crate) wall_solid: SolidSide,
    pub(crate) cyl_solid: SolidSide,
    /// Is the cylinder's side of the wall plane the wall **face**'s material side? Read from that
    /// face's own outward statement, not from the class frame — one class can carry faces of both
    /// operands with opposite material sides.
    pub(crate) lens_in_wall_solid: bool,
    /// `+1` material inside the cylinder (a boss), `-1` outside (a bore) — the lateral face's own
    /// [`CylFaceInfo::orient_sign`].
    pub(crate) cyl_orient: i8,
    /// **Something above could not be stated, and the refusal that says what**: the wall face's
    /// own, where the footprint reader could not read it (the name [`wall_faces_clear`] raises for
    /// that face), else [`RejectReason::CylinderGateUndecided`] (a face with no world description
    /// of its own, an overflow). The gate raises it; no row the verdict sees has one.
    undecided: Option<BoolError>,
    /// **Where both solids run through the tangent line, on this row's ground** — closed pieces
    /// `[lo, hi]` in the axis parameter ([`axis_param_of_plane`]'s scale), within those where the
    /// closed wall face meets the line inside the lateral face's open stretches in positive length
    /// ([`line_runs_through_face`]), and where the wall's solid lies on both sides of the line in the wall's plane
    /// and the cylinder's solid on both sides of it on the cylinder (a seam is no boundary). This
    /// is the verdict's ground — and the self-touch check's tangent hand-off's: a point of the lateral
    /// it finds on the line is this verdict's question exactly there, and anywhere else on the line
    /// the contact is the check's own. A valid solid meets a line with one face running through it or two
    /// faces ending on it, so here no third face reaches the line and nothing in the arrangement
    /// draws it — the three regions are the whole picture — and elsewhere on the contact one of
    /// the solids has an edge on the line, which the arrangement draws and the structure judges.
    ///
    /// **Interior, holes and span included** ([`line_runs_through_face`]). Corners on both
    /// sides of the line say this only of a convex face without holes: a U-shaped face or one with
    /// a hole has them while the line runs through the gap alone, and a face can cross the line
    /// below the lateral and not beside it. Read as that, a row convicts `block − boss` of a pinch
    /// where the boss stands in a wall's slot and cuts the material beside it (measured). The span is
    /// open, and a face flush with the cap planes still has its interior inside it (the frozen
    /// `bore-slab` shape). A face meeting the line at a corner meets it in no length and so in no
    /// piece — a point tangency is a valid solid, and [`nacre_exact::cylinder_strip_side`]'s
    /// `StripSide::Inside` for `U = 0` (*"merely conservative"* on a tangent plane's empty strip)
    /// says nothing on its own here.
    ///
    /// `None` where a face could not be read — the gate refuses such a row: no default is safe,
    /// since reading a solid that ends on the line as running through convicts it of three regions
    /// it does not have, and reading one that runs through as ending passes a pinch nothing else
    /// sees.
    pub(crate) through: Option<Vec<[nacre_exact::QuadVal; 2]>>,
    /// **Where the wall's solid has a face on the line's `+τ` side and on its `−τ` side** — `τ = m × r`
    /// the lateral's `+θ` tangent at the line, `r` from the axis out to it; closed pieces in the axis
    /// parameter, over this row's lateral's stretches, from every face of the wall's solid on the
    /// class. [`Self::through`] is the stretch where both sides hold one, met with the lateral's own;
    /// the record of the line ([`SharedRuling::tangent_halves`]) gathers these for the chart, where
    /// one side alone is a face that **ends** on the line — its edge, a spoke of the line's star.
    /// `None` where a face of the wall's solid could not be read (so `through` is `None` too) or the
    /// lens' side could not be stated (`undecided`) — rows the gate refuses either way.
    halves: Option<[Vec<[nacre_exact::QuadVal; 2]>; 2]>,
    /// The tangency point — `RejectWhere`'s doc asks for a witness that is actually there: the
    /// middle of the first stretch where the line runs through the wall face, else of the first it
    /// runs along, when that stretch has rational ends; otherwise the middle of the lateral face's
    /// own span.
    pub(crate) witness: Point3,
}

/// The foot of the perpendicular from the axis origin to the plane — **rational**, because it is
/// `o - ((n·o + d)/(n·n))·n` and nothing there leaves the field. This is the tangency line's base
/// point; the line itself is that point plus `t·m`.
fn tangency_foot(
    coeffs: &[nacre_exact::Rat; 4],
    o: &[nacre_exact::Rat; 3],
) -> Option<[nacre_exact::Rat; 3]> {
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let mut num = coeffs[3];
    for k in 0..3 {
        num = num.checked_add(n[k].checked_mul(o[k])?)?;
    }
    let den = nacre_exact::dot3_rat(&n, &n)?;
    // `Rat` has no division: the reciprocal is the exact inverse and `Rat::new` refuses `0`,
    // which is precisely the degenerate normal this must not divide by.
    let q = num.checked_mul(nacre_exact::Rat::new(den.denom(), den.numer())?)?;
    let mut out = [nacre_exact::Rat::from_int(0); 3];
    for k in 0..3 {
        out[k] = o[k].checked_sub(q.checked_mul(n[k])?)?;
    }
    Some(out)
}

/// A lateral face's stretches of one ruling, in the axis parameter (`lateral_cover_on_ruling`).
type Stretches = Vec<[nacre_exact::Rat; 2]>;

/// One lateral face as [`tangency_rows`] reads it against the tangent line.
struct LateralOnLine<'a> {
    /// Its row in the face table.
    ix: usize,
    cf: &'a CylFaceInfo,
    /// Its own stretches of the line; `None` where they could not be read.
    cover: Option<Stretches>,
    /// The stretches a row reads: `cover`, else the face's axial span.
    stretches: Option<Stretches>,
}

/// **The rows a tangent `(class, cylinder)` pair writes** — one per (wall face, lateral face) pair
/// whose lateral face lies on the line ([`Tangency`]'s doc says how that is read) and whose wall
/// face did not clear the footprint. A face that *did* clear cannot reach the tangent line, so it
/// contributes nothing; that is why collecting here neither widens nor narrows the rule.
///
/// **This walk cannot fail.** [`wall_faces_clear`] short-circuits on its first non-clearing
/// face and raises on shapes it cannot read; a raise propagated from the middle of this walk
/// would outrank the gate's other refusals by the order the faces happen to be in. Every failure
/// here becomes `undecided` on the row instead, carrying its refusal, and the gate raises the
/// first such row's after its other questions — the same name for an unreadable face wherever it
/// stands.
#[allow(clippy::too_many_arguments)]
fn tangency_rows(
    model: &Model,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
    c: usize,
    ci: usize,
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
    surf: Handle<Surface>,
    owner: SolidSide,
    jd: &Judge<'_, WorkingPlane>,
) -> Vec<Tangency> {
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let side_of = |i: usize| {
        if i < n_a { SolidSide::A } else { SolidSide::B }
    };
    let base = tangency_foot(coeffs, &o);
    // Which side of the wall plane the cylinder is on — the lens' side. A tangency puts the
    // whole cylinder on one side, and `0` would mean the axis lies *in* the plane, which needs
    // `r = 0`; `CylinderDef` refuses that at construction (`NonPositiveRadius`, "positive by
    // construction"). Filtered anyway rather than assumed: a `Some(0)` here would compare unequal
    // to every material side and hand the verdict a silent `false`, and this file's own rule is
    // that an answer it cannot state becomes `undecided`, never a default. The sign is the
    // class's world name's own, the frame `mat_side` below is read in.
    let lens_side = jd.planes[c]
        .world
        .as_ref()
        .map(|w| nacre_exact::plane_residual_sign(&w.name, o))
        .filter(|&s| s != 0);
    // Built **once**, not once per face (the rule `wall_faces_clear`'s doc gives).
    let spans = lateral_spans(faces, surf);
    // The lateral faces of this cylinder (`LateralOnLine`), read once per pair: where each lies
    // on the tangent line depends on the line and the lateral, not on the wall face. A cover that could not be read keeps the face's bounding
    // rectangle.
    let laterals: Vec<LateralOnLine<'_>> = faces
        .iter()
        .enumerate()
        .filter_map(|(i, row)| match row {
            FaceRow::Cylinder(cf) if cf.surf == surf => {
                let cover = base
                    .as_ref()
                    .zip(cf.face)
                    .zip(cf.def.as_ref())
                    .and_then(|((b, fh), d)| lateral_cover_on_ruling(model, model.face(fh), d, b));
                let stretches = match &cover {
                    Some(c) => Some(c.clone()),
                    None => cf.stretch.map(|sp| vec![sp]),
                };
                Some(LateralOnLine {
                    ix: i,
                    cf,
                    cover,
                    stretches,
                })
            }
            _ => None,
        })
        .collect();
    // **Where the cylinder's solid runs through the line** — its lateral faces on this cylinder
    // holding the generators on both sides of the line (a seam is no boundary); `None` where a face
    // could not be read. Read once, and only for a pair that writes a row: most tangent pairs'
    // wall faces clear the footprint, and reading every lateral here cost the census 1.1% of its
    // instructions.
    let read_lateral_through = || {
        let b = base.as_ref()?;
        let mut sides: [Vec<[nacre_exact::QuadVal; 2]>; 2] = [Vec::new(), Vec::new()];
        for row in faces {
            let FaceRow::Cylinder(cf) = row else { continue };
            if cf.surf != surf {
                continue;
            }
            let read = lateral_sides_on_ruling(model, model.face(cf.face?), cf.def.as_ref()?, b)?;
            for (to, from) in sides.iter_mut().zip(read) {
                to.extend(
                    from.into_iter()
                        .map(|sp| sp.map(nacre_exact::QuadVal::from_rat)),
                );
            }
        }
        meet_unions(&sides[0], &sides[1])
    };
    let mut lateral_through: Option<Option<Vec<[nacre_exact::QuadVal; 2]>>> = None;
    let mut out = Vec::new();
    for (fi_ix, row) in faces.iter().enumerate() {
        let (ClassIx::Plane(k), FaceRow::Plane(fi)) = (plane_ix[fi_ix], row) else {
            continue;
        };
        if k != c {
            continue;
        }
        let fh = fi.face.expect("collect_planes yields real faces");
        // **Same-solid pairs write no row** (a valid solid's faces do not meet — the pair loop's
        // premise). A wall tangent to its own solid's cylinder is that solid's smooth edge — a
        // fillet — not a contact between operands. The judge would skip such a row anyway (the
        // face ends on the line, so no piece runs through), but a row it could not decide would
        // refuse the boolean by a name that is about nothing.
        if matches!(
            (side_of(fi_ix), owner),
            (SolidSide::A, SolidSide::A) | (SolidSide::B, SolidSide::B)
        ) {
            continue;
        }
        let read = {
            let got = face_clears_footprint(model, model.face(fh), coeffs, &o, &m, r2, &spans);
            #[cfg(feature = "tangency-trace")]
            if got.is_err() {
                let f = model.face(fh);
                let arcs = f
                    .outer
                    .half_edges
                    .iter()
                    .filter(|he| !matches!(model.edge_curve(he.edge), nacre_geom::Curve::Line(_)))
                    .count();
                #[allow(clippy::print_stderr)]
                {
                    eprintln!(
                        "TFACE unread edges={} arcs={} holes={}",
                        f.outer.half_edges.len(),
                        arcs,
                        f.inner.len()
                    );
                }
            }
            got
        };
        // **A face this could not read is «unknown», not «innocent».** Folded to «did not clear»,
        // it would write a row whose stretches through the line were **never computed**, and an
        // empty answer acquits. So the row carries the face's own refusal — the name
        // `wall_faces_clear` raises for such a face (this function's doc says why the order of
        // the class's faces does not change it).
        let (cleared, unread) = match read {
            Ok(cleared) => (cleared, None),
            Err(e) => (false, Some(e)),
        };
        if cleared {
            continue; // a face that clears the footprint cannot reach the tangent line
        }
        // The wall face's material side, w.r.t. the class's coefficients. The outward direction
        // is `orient_sign` times the way the face's plane faces, material lies opposite it, so the
        // answer is `−orient_sign · κ` where `κ` relates the plane's facing to the class's rational
        // name. `world_rat` is **not** that relation — the seed plane `x = 0` carries
        // `world_rat = (1,0,0,0)` while it faces `−x`. The face's own world sense is
        // ([`WorldName`]), where its world name *is* the class's; a face that joined the class by
        // its definitions, with no world name of its own, is left `undecided`.
        let mat_side = fi
            .world
            .as_ref()
            .filter(|w| w.name.narrow() == Some(coeffs))
            .map(|w| -fi.orient_sign * w.sense.sign());
        for LateralOnLine {
            ix: cy_ix,
            cf,
            cover,
            stretches,
        } in &laterals
        {
            // The line must lie on this lateral **face** (`Tangency`'s doc): read off the
            // face's stretches of the line when it could be, else off its angular extent. `None`
            // (arithmetic) keeps the row, undecided.
            let on_face = match cover {
                Some(c) => Some(!c.is_empty()),
                None => base.as_ref().and_then(|b| {
                    let line = Footprint {
                        span: None,
                        theta: Theta::at_line(b, &o, &m)?,
                    };
                    Some(!line_apart(def, &line, &cf.footprint)?)
                }),
            };
            if on_face == Some(false) {
                #[cfg(feature = "tangency-trace")]
                #[allow(clippy::print_stderr)]
                {
                    eprintln!("TSKIP off the lateral face");
                }
                continue;
            }
            // **Per lateral, over that lateral's own stretches of the line** (`Tangency`'s
            // doc). A face the clearance road could not read is not asked again: it keeps the row,
            // carrying its refusal, and reads as a contact through the interior — dropped here,
            // that refusal would vanish with it (`line_runs_through_face`); and so does a row any
            // stretch of which could not be read.
            let per_stretch = base
                .as_ref()
                .zip(stretches.clone())
                .filter(|_| unread.is_none())
                .and_then(|(b, stretches)| {
                    stretches
                        .into_iter()
                        .map(|sp| {
                            line_runs_through_face(model, model.face(fh), coeffs, b, &m, sp)
                                .map(|c| (sp, c))
                        })
                        .collect::<Option<Vec<_>>>()
                });
            let chosen = per_stretch.as_ref().and_then(|v| {
                v.iter()
                    .find(|(_, c)| c.runs_through)
                    .or_else(|| v.iter().find(|(_, c)| c.touches))
                    .map(|(sp, c)| (*sp, c.contact_at))
            });
            let contact = match &per_stretch {
                Some(v) => LineContact {
                    touches: v.iter().any(|(_, c)| c.touches),
                    runs_through: v.iter().any(|(_, c)| c.runs_through),
                    contact_at: chosen.and_then(|(_, at)| at),
                    contact: v
                        .iter()
                        .flat_map(|(_, c)| c.contact.iter().copied())
                        .collect(),
                    sides: v
                        .iter()
                        .flat_map(|(_, c)| c.sides.iter().copied())
                        .collect(),
                },
                // Read as a contact through the interior — over every stretch of this lateral, the
                // ground that reading claims.
                None => LineContact {
                    touches: true,
                    runs_through: true,
                    contact_at: None,
                    contact: stretches
                        .iter()
                        .flatten()
                        .map(|sp| sp.map(nacre_exact::QuadVal::from_rat))
                        .collect(),
                    sides: stretches.iter().flatten().map(|_| [true; 2]).collect(),
                },
            };
            if !contact.touches {
                #[cfg(feature = "tangency-trace")]
                #[allow(clippy::print_stderr)]
                {
                    eprintln!("TSKIP the wall face does not reach the line on the lateral");
                }
                continue;
            }
            // **Where both solids run through the line, on this row's ground** — this wall face's
            // contact, where the wall's solid lies on both sides of the line in the class's plane
            // (over all its faces there) and the cylinder's solid on both sides of it on the
            // cylinder. `None` where a face could not be read.
            let lateral_through = lateral_through.get_or_insert_with(read_lateral_through);
            // The wall's solid's faces on the class, by the side of the line they lie on (`±w`).
            let wall_sides = (|| {
                let b = base.as_ref()?;
                let mut wall: [Vec<[nacre_exact::QuadVal; 2]>; 2] = [Vec::new(), Vec::new()];
                for (j, row) in faces.iter().enumerate() {
                    let (ClassIx::Plane(kj), FaceRow::Plane(fj)) = (plane_ix[j], row) else {
                        continue;
                    };
                    if kj != c || side_of(j) != side_of(fi_ix) {
                        continue;
                    }
                    let face = model.face(fj.face?);
                    for &sp in stretches.as_ref()? {
                        let read = line_runs_through_face(model, face, coeffs, b, &m, sp)?;
                        for (piece, on) in read.contact.iter().zip(&read.sides) {
                            for (to, &here) in wall.iter_mut().zip(on) {
                                if here {
                                    to.push(*piece);
                                }
                            }
                        }
                    }
                }
                let [plus, minus] = wall;
                Some([joined(plus)?, joined(minus)?])
            })();
            let through = (|| {
                let wall = wall_sides.as_ref()?;
                let both =
                    meet_unions(&meet_unions(&wall[0], &wall[1])?, lateral_through.as_ref()?)?;
                per_stretch.as_ref().and_then(|v| {
                    let own: Vec<[nacre_exact::QuadVal; 2]> = v
                        .iter()
                        .flat_map(|(_, c)| c.contact.iter().copied())
                        .collect();
                    meet_unions(&own, &both)
                })
            })();
            // `w = n × m` and `τ = m × r`; the line's foot is `o − q·n` with `q = (n·o + d)/|n|²`, so
            // `r = −q·n` and `τ = q·w`: `+w` is `+τ` exactly where the axis lies on the plane's
            // positive side — the lens' side.
            let halves = wall_sides.zip(lens_side).map(|([plus_w, minus_w], lens)| {
                if lens > 0 {
                    [plus_w, minus_w]
                } else {
                    [minus_w, plus_w]
                }
            });
            let witness_span = chosen.map(|(sp, _)| sp).or(cf.stretch);
            let witness = base.as_ref().zip(witness_span).and_then(|(b, span)| {
                let mid = match contact.contact_at {
                    Some(v) => v,
                    None => span[0]
                        .checked_add(span[1])?
                        .checked_mul(nacre_exact::Rat::new(1, 2)?)?,
                };
                let mut p = [0.0f64; 3];
                for j in 0..3 {
                    p[j] = b[j].checked_add(mid.checked_mul(m[j])?)?.to_f64();
                }
                Some(Point3::from_array(p))
            });
            let undecided = unread.or_else(|| {
                (base.is_none()
                    || on_face.is_none()
                    || lens_side.is_none()
                    || mat_side.is_none()
                    || witness.is_none()
                    || cf.def.is_none())
                .then(|| reject(RejectReason::CylinderGateUndecided))
            });
            #[cfg(feature = "tangency-trace")]
            #[allow(clippy::print_stderr)]
            {
                eprintln!(
                    "TROW runs_through={} through={:?} undecided={}",
                    contact.runs_through,
                    through.as_ref().map(Vec::len),
                    undecided.is_some()
                );
            }
            out.push(Tangency {
                wall: c,
                cyl: ci,
                wall_solid: side_of(fi_ix),
                cyl_solid: side_of(*cy_ix),
                lens_in_wall_solid: lens_side.is_some() && lens_side == mat_side,
                cyl_orient: cf.orient_sign,
                undecided,
                through,
                halves,
                witness: witness.unwrap_or(Point3::from_array([f64::NAN; 3])),
            });
        }
    }
    // A tangency of a solid with **itself** is not this operation's business: for one solid to
    // own both the wall face and a lateral tangent to it along that line, its own material would
    // have to be the two wedges alone or the lens and the far side alone — and both of those *are*
    // the pinch. Such an operand is invalid before any boolean runs, so dropping the row discards
    // no information about a valid input.
    out.retain(|t| t.wall_solid != t.cyl_solid);
    out
}

/// **What the tangent line meets of one wall face, over one lateral face's axis span** — the
/// answer [`line_runs_through_face`] gives when it could read the face.
#[derive(Clone, Debug)]
struct LineContact {
    /// The closed face meets the open stretch of the line in positive length — through its
    /// interior or along its boundary. Without this the face touches nothing a row could speak of.
    touches: bool,
    /// The line runs through the face's **interior** over positive length.
    runs_through: bool,
    /// The middle of the first stretch where the face meets the line — through the interior when
    /// there is one, else along the boundary — when that stretch has rational ends: the row's
    /// witness, which then sits on the contact.
    contact_at: Option<nacre_exact::Rat>,
    /// **Where** the closed face meets the open stretch in positive length — the closed pieces
    /// `[lo, hi]` of the line, in the axis parameter, that make `touches` true, through the
    /// interior or along the boundary alike.
    contact: Vec<[nacre_exact::QuadVal; 2]>,
    /// Beside each piece of `contact`, whether the face lies on the line's `+w` side and on its
    /// `−w` side there (`w = n × m`, the direction across the line in the face's plane): both is
    /// the line running through the face's interior, one is the line along its boundary.
    sides: Vec<[bool; 2]>,
}

/// **`(∪a) ∩ (∪b)`** for closed pieces of the line, as sorted disjoint pieces of positive length —
/// `None` where two ends could not be compared.
fn meet_unions(
    a: &[[nacre_exact::QuadVal; 2]],
    b: &[[nacre_exact::QuadVal; 2]],
) -> Option<Vec<[nacre_exact::QuadVal; 2]>> {
    use core::cmp::Ordering;
    let mut met: Vec<[nacre_exact::QuadVal; 2]> = Vec::new();
    for x in a {
        for y in b {
            let lo = if x[0].cmp_value(&y[0])? == Ordering::Less {
                y[0]
            } else {
                x[0]
            };
            let hi = if x[1].cmp_value(&y[1])? == Ordering::Less {
                x[1]
            } else {
                y[1]
            };
            if lo.cmp_value(&hi)? == Ordering::Less {
                met.push([lo, hi]);
            }
        }
    }
    joined(met)
}

/// **Closed pieces of the line, sorted by the lower end with overlapping and touching pieces
/// joined** — `None` where two ends could not be compared.
fn joined(met: Vec<[nacre_exact::QuadVal; 2]>) -> Option<Vec<[nacre_exact::QuadVal; 2]>> {
    use core::cmp::Ordering;
    let mut sorted: Vec<[nacre_exact::QuadVal; 2]> = Vec::with_capacity(met.len());
    for p in met {
        let mut at = sorted.len();
        for (i, q) in sorted.iter().enumerate() {
            if p[0].cmp_value(&q[0])? == Ordering::Less {
                at = i;
                break;
            }
        }
        sorted.insert(at, p);
    }
    let mut out: Vec<[nacre_exact::QuadVal; 2]> = Vec::with_capacity(sorted.len());
    for p in sorted {
        match out.last_mut() {
            Some(last) if p[0].cmp_value(&last[1])? != Ordering::Greater => {
                if last[1].cmp_value(&p[1])? == Ordering::Less {
                    last[1] = p[1];
                }
            }
            _ => out.push(p),
        }
    }
    Some(out)
}

/// How one loop of the wall face meets the tangent line.
enum LoopOnLine {
    /// The loop's pieces all lie strictly on one side: its region meets the line nowhere.
    Clear,
    /// A disk the line crosses: its region on the line is the open chord between the two values.
    Chord(nacre_exact::QuadVal, nacre_exact::QuadVal),
    /// A polygon: the values where its steps cross the line, under each of the two half-open
    /// readings — `[0]` counts a piece on the line as below it, `[1]` as above.
    Crossings([Vec<nacre_exact::QuadVal>; 2]),
}

/// **Where does the tangent line run through a face whose plane holds it, over one stretch of the
/// line a lateral face covers?** — of the wall face, where a row's contact is; of every face of the
/// wall's solid on its class, where that solid lies on each side of the line
/// ([`Tangency::through`]).
///
/// The line is `foot + v·m` (`foot` is the perpendicular foot from the axis origin onto the wall,
/// so `v` is the axis parameter itself — the scale the span is in), and the question is about the
/// **open** stretch `span[0] < v < span[1]`: the lateral face is there, and a row reads each of the
/// face's stretches in turn. `coeffs` are the face's own plane's, and `foot` lies on it.
///
/// **The line is a strip of width zero.** A piece's side is asked as the strip test asks it
/// ([`nacre_exact::cylinder_strip_side`]), about `foot` at radius `0`: the strip a radius-`0`
/// cylinder on the line cuts from any plane holding the line is the line itself. The wall's own
/// strip about the axis is that line already — its plane is `r` from the axis, so the half-width is
/// `0` — so the wall reads as it does about the axis; a secant plane holding the line cuts a strip
/// of positive width about the axis, of which the line is one edge, and only the radius-`0` strip
/// reads it.
///
/// **It is point-in-polygon with the line itself as the ray.** Every loop — the outer and
/// every hole — is read through [`corner_of`], and a step whose ends lie on opposite sides of the
/// line crosses it once; the parity of the crossings below a point says whether the point is in
/// that loop's region, and the face's region is inside the outer and outside every hole. The
/// crossing rule is [`nacre_geom::intersect::ray_straddle`]'s, the one half-open rule the ring
/// walks read: a piece *on* the line is «not above». Asked twice — once as is, once with the
/// sides swapped — it reads the face along the line pushed off to either side, and a stretch is
/// **interior** exactly when both readings put it inside: a face that ends on the line along an
/// edge, or meets it at one corner, is inside on one side only. The arrangement's own sweep
/// (`trace_transversal_face`) answers the same kind of question, but `planes` sits below the
/// engine and runs before it.
///
/// **The crossing is placed from one end of the step and the step's direction** — the edge's
/// two ends when both are rational, else its carrier pair (`Edge::surfaces`): a plane neighbour's
/// meet with the wall, a cylinder neighbour's ruling. That is linear in the end's coordinates, so
/// a pierce corner's value keeps its own one radical.
///
/// `None` is «could not read the face»: a step whose crossing could not be placed, an arc in a loop
/// the line reaches, an outer piece [`corner_of`] cannot spell or that is off the class's plane.
/// The caller reads that as a contact through the interior, never a silent acquittal. A **hole**
/// that cannot be read is left out instead: the face only grows, so the answer can only lean
/// towards a contact.
#[allow(clippy::too_many_arguments)]
fn line_runs_through_face(
    model: &Model,
    face: &Face,
    coeffs: &[nacre_exact::Rat; 4],
    foot: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    span: [nacre_exact::Rat; 2],
) -> Option<LineContact> {
    use super::corners::Corner;
    use core::cmp::Ordering;
    use nacre_exact::{BigRat, Orient, QuadVal, Rat, StripSide, cross3_rat, dot3_rat};
    debug_assert!(
        nacre_exact::point_on_plane_exact(coeffs, &nacre_exact::MeetPoint::Narrow(*foot)),
        "the line lies on the face's plane"
    );
    let line_r2 = BigRat::zero();
    let zero = Rat::from_int(0);
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    // `w` is the direction across the line in the wall's plane; `StripSide::Plus` is its side.
    let w = cross3_rat(&n, m)?;
    let mm = dot3_rat(m, m)?;
    let inv_mm = Rat::new(mm.denom(), mm.numer())?;
    // A point's `(u, v)`: its offset across the line along `w`, and its axis parameter.
    let uv = |base: &[Rat; 3], dir: &[Rat; 3], s: QuadVal| -> Option<(QuadVal, QuadVal)> {
        let mut rel = [zero; 3];
        for k in 0..3 {
            rel[k] = base[k].checked_sub(foot[k])?;
        }
        let along = |e: &[Rat; 3], scale: Rat| -> Option<QuadVal> {
            QuadVal::from_rat(dot3_rat(&rel, e)?.checked_mul(scale)?)
                .checked_add(&s.checked_mul_rat(dot3_rat(dir, e)?.checked_mul(scale)?)?)
        };
        Some((along(&w, Rat::from_int(1))?, along(m, inv_mm)?))
    };
    let point_uv = |c: &Corner| -> Option<(QuadVal, QuadVal)> {
        match c {
            Corner::Rational(p) => uv(p.narrow()?, &[zero; 3], QuadVal::from_rat(zero)),
            Corner::Pierce(line, s) => uv(&line.base(), &line.dir(), *s),
            Corner::Round { .. } => None,
        }
    };
    // The direction of the straight edge `he`, `a → b` being its two ends' pieces.
    let step_dir = |he: &nacre_topo::HalfEdge, a: &Corner, b: &Corner| -> Option<[Rat; 3]> {
        if let (Corner::Rational(pa), Corner::Rational(pb)) = (a, b)
            && let (Some(pa), Some(pb)) = (pa.narrow(), pb.narrow())
        {
            let mut d = [zero; 3];
            for k in 0..3 {
                d[k] = pb[k].checked_sub(pa[k])?;
            }
            return Some(d);
        }
        let [s0, s1] = model.edge(he.edge).surfaces;
        let other = match (s0 == face.surface, s1 == face.surface) {
            (true, false) => s1,
            (false, true) => s0,
            _ => return None,
        };
        match model.surface(other) {
            nacre_topo::Surface::Plane { .. } => {
                let q = model.world_plane_name(other)?.narrow().copied()?;
                cross3_rat(&n, &[q[0], q[1], q[2]])
            }
            nacre_topo::Surface::Cylinder { .. } => Some(model.world_cylinder_def(other)?.dir()),
        }
    };
    // Each loop's meeting with the line; `Err` is a loop that could not be read.
    let read_loop = |lp: &nacre_topo::Loop| -> Result<LoopOnLine, ()> {
        let pieces: Vec<Corner> = lp
            .half_edges
            .iter()
            .map(|he| corner_of(model, face, lp, he).map_err(|_| ()))
            .collect::<Result<_, ()>>()?;
        // The strip test's precondition: a piece off the class's plane is judged against a
        // distance that is not its own — the two-descriptions hazard `face_clears_footprint`
        // refuses on. Unread, not guessed.
        if !pieces.iter().all(|c| c.on_plane(coeffs)) {
            return Err(());
        }
        let sides: Vec<StripSide> = pieces
            .iter()
            .map(|c| c.strip_side(coeffs, foot, m, &line_r2))
            .collect();
        if sides.iter().all(|s| *s == StripSide::Plus)
            || sides.iter().all(|s| *s == StripSide::Minus)
        {
            return Ok(LoopOnLine::Clear);
        }
        if let [
            Corner::Round {
                centre,
                rho2,
                theta: Theta::Whole,
                ..
            },
        ] = pieces.as_slice()
        {
            if sides[0] != StripSide::Crosses {
                // On the line at one point, or nowhere: no positive length either way.
                return Ok(LoopOnLine::Clear);
            }
            // The chord: centre `c`, half-length `h` in the axis scale, `h² = (ρ² − d²)/(m·m)`
            // with `d² = u_c²/(w·w)`.
            let chord = (|| {
                let (u, v) = uv(centre, &[zero; 3], QuadVal::from_rat(zero))?;
                let (u, v) = (u.as_rat()?, v.as_rat()?);
                let ww = dot3_rat(&w, &w)?;
                let d2 = u
                    .checked_mul(u)?
                    .checked_mul(Rat::new(ww.denom(), ww.numer())?)?;
                let h2 = rho2.narrow()?.checked_sub(d2)?.checked_mul(inv_mm)?;
                Some(LoopOnLine::Chord(
                    QuadVal::new(v, Rat::from_int(-1), h2)?,
                    QuadVal::new(v, Rat::from_int(1), h2)?,
                ))
            })();
            return chord.ok_or(());
        }
        if pieces.iter().any(|c| matches!(c, Corner::Round { .. })) {
            return Err(()); // an arc in a loop the line reaches: no crossing value for it here
        }
        let k = pieces.len();
        let mut out: [Vec<QuadVal>; 2] = [Vec::new(), Vec::new()];
        for i in 0..k {
            let (a, b) = (&pieces[i], &pieces[(i + 1) % k]);
            let (sa, sb) = (sides[i], sides[(i + 1) % k]);
            let orient = |s: StripSide, flip: bool| match (s, flip) {
                (StripSide::Plus, false) | (StripSide::Minus, true) => Orient::Positive,
                (StripSide::Minus, false) | (StripSide::Plus, true) => Orient::Negative,
                _ => Orient::Zero,
            };
            let crosses: [bool; 2] = core::array::from_fn(|r| {
                nacre_geom::intersect::ray_straddle(orient(sa, r == 1), orient(sb, r == 1))
                    .is_some()
            });
            if !crosses[0] && !crosses[1] {
                continue;
            }
            // Where the step meets `u = 0`: from an end `(u, v)` along `d`, `v − u·(d·m)/((d·w)(m·m))`.
            let value = (|| {
                let d = step_dir(&lp.half_edges[i], a, b)?;
                let dw = dot3_rat(&d, &w)?;
                if dw == zero {
                    return None;
                }
                let kk = dot3_rat(&d, m)?
                    .checked_mul(Rat::new(dw.denom(), dw.numer())?)?
                    .checked_mul(inv_mm)?;
                let (u, v) = point_uv(a).or_else(|| point_uv(b))?;
                v.checked_sub(&u.checked_mul_rat(kk)?)
            })()
            .ok_or(())?;
            for r in 0..2 {
                if crosses[r] {
                    out[r].push(value);
                }
            }
        }
        Ok(LoopOnLine::Crossings(out))
    };
    let outer = read_loop(&face.outer).ok()?;
    if matches!(outer, LoopOnLine::Clear) {
        return Some(LineContact {
            touches: false,
            runs_through: false,
            contact_at: None,
            contact: Vec::new(),
            sides: Vec::new(),
        });
    }
    let holes: Vec<LoopOnLine> = face
        .inner
        .iter()
        .filter_map(|lp| read_loop(lp).ok())
        .collect();
    // The breakpoints along the line, in order: every crossing and chord end, and the span's ends.
    let mut points: Vec<QuadVal> = vec![QuadVal::from_rat(span[0]), QuadVal::from_rat(span[1])];
    for l in core::iter::once(&outer).chain(holes.iter()) {
        match l {
            LoopOnLine::Clear => {}
            LoopOnLine::Chord(lo, hi) => points.extend([*lo, *hi]),
            LoopOnLine::Crossings(c) => points.extend(c[0].iter().chain(c[1].iter()).copied()),
        }
    }
    let mut sorted: Vec<QuadVal> = Vec::with_capacity(points.len());
    for p in points {
        let mut at = sorted.len();
        for (i, q) in sorted.iter().enumerate() {
            match p.cmp_value(q)? {
                Ordering::Equal => {
                    at = usize::MAX;
                    break;
                }
                Ordering::Less => {
                    at = i;
                    break;
                }
                Ordering::Greater => {}
            }
        }
        if at != usize::MAX {
            sorted.insert(at, p);
        }
    }
    // Is the gap just above `lo` inside the loop's region, under reading `r`?
    let inside = |l: &LoopOnLine, lo: QuadVal, hi: QuadVal, r: usize| -> Option<bool> {
        Some(match l {
            LoopOnLine::Clear => false,
            LoopOnLine::Chord(a, b) => {
                a.cmp_value(&lo)? != Ordering::Greater && hi.cmp_value(b)? != Ordering::Greater
            }
            LoopOnLine::Crossings(c) => {
                let mut below = 0usize;
                for x in &c[r] {
                    if x.cmp_value(&lo)? != Ordering::Greater {
                        below += 1;
                    }
                }
                below % 2 == 1
            }
        })
    };
    let (s0, s1) = (QuadVal::from_rat(span[0]), QuadVal::from_rat(span[1]));
    let mut out = LineContact {
        touches: false,
        runs_through: false,
        contact_at: None,
        contact: Vec::new(),
        sides: Vec::new(),
    };
    for pair in sorted.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        if lo.cmp_value(&s0)? == Ordering::Less || hi.cmp_value(&s1)? == Ordering::Greater {
            continue;
        }
        let mut reading = [false; 2];
        for (r, slot) in reading.iter_mut().enumerate() {
            let mut inn = inside(&outer, lo, hi, r)?;
            for h in &holes {
                inn = inn && !inside(h, lo, hi, r)?;
            }
            *slot = inn;
        }
        let middle = || {
            lo.as_rat()
                .zip(hi.as_rat())
                .and_then(|(a, b)| a.checked_add(b)?.checked_mul(Rat::new(1, 2)?))
        };
        if reading[0] || reading[1] {
            out.contact.push([lo, hi]);
            out.sides.push(reading);
        }
        if (reading[0] || reading[1]) && !out.touches {
            out.touches = true;
            out.contact_at = middle();
        }
        if reading[0] && reading[1] && !out.runs_through {
            out.runs_through = true;
            out.contact_at = middle();
        }
    }
    Some(out)
}

/// The axis-parameter spans this cylinder's lateral faces occupy — one per face, in the scale
/// `axis_param_of_plane` produces.
///
/// **Empty means "unusable", never "nothing in the way".** A span this returns is a rectangle
/// the wall rule must miss, and "every one of no rectangles is missed" is a pass no caller should
/// hand out by default — so a face whose span could not be stated, and a cylinder with no lateral
/// face in the table at all, both come back empty and leave the axis unused. The second could be
/// argued safe (no lateral face, no band, nothing to protect), but that argument rests on the face
/// table being complete here, which is a separate premise from the one this function is about.
fn lateral_spans(faces: &[FaceRow], surf: Handle<Surface>) -> Vec<[nacre_exact::Rat; 2]> {
    lateral_footprints(faces, surf)
        .into_iter()
        .map(|f| f.span.expect("a listed footprint has a span"))
        .collect()
}

/// **The ranges of a cylinder class's lateral faces** ([`CylFaceInfo::stretch`]), one per face —
/// `None` if any face's range is not a pair of heights.
fn lateral_stretches(
    faces: &[FaceRow],
    surf: Handle<Surface>,
) -> Option<Vec<[nacre_exact::Rat; 2]>> {
    faces
        .iter()
        .filter_map(|row| match row {
            FaceRow::Cylinder(cf) if cf.surf == surf => Some(cf.stretch),
            _ => None,
        })
        .collect()
}

/// **The footprints of a cylinder class's lateral faces** — [`lateral_spans`]' rule with the
/// angular extent alongside: empty if any face's span could not be stated, else one
/// footprint per face.
fn lateral_footprints(faces: &[FaceRow], surf: Handle<Surface>) -> Vec<Footprint> {
    let mut out: Vec<Footprint> = Vec::new();
    for row in faces {
        let FaceRow::Cylinder(cf) = row else { continue };
        // Matched by **handle**, not by geometry: two operands may state the same cylinder
        // twice, and each statement gets its own rows, its own bands clipped to its own spans,
        // and its own turn through the gate's cylinder loop.
        if cf.surf != surf {
            continue;
        }
        let Some(span) = cf.footprint.span else {
            return Vec::new();
        };
        // The reader below asks "every vertex at or below `span[0]`, or every one at or above
        // `span[1]`", which is only the interval's outside while the pair is ordered — reversed,
        // that same phrasing reads "outside the *union* of two half-lines" and would clear a face
        // sitting squarely in the band. `lateral_axial_reach` orders it; this is where that is
        // relied on, so this is where it is said.
        debug_assert!(span[0] <= span[1], "a span is stated low end first");
        out.push(cf.footprint);
    }
    out
}

fn lateral_reach(
    def: &nacre_topo::CylinderDef,
    fp: &Footprint,
    d: &[nacre_exact::Rat; 3],
) -> Option<Reach> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let dm = nacre_exact::dot3_rat(d, &m)?;
    let base = nacre_exact::dot3_rat(d, &o)?;
    let (lo, hi) = if dm == zero {
        (base, base)
    } else {
        let [s0, s1] = fp.span?;
        let (a, b) = (s0.checked_mul(dm)?, s1.checked_mul(dm)?);
        (base.checked_add(a.min(b))?, base.checked_add(a.max(b))?)
    };
    let (lo_off, rho2_lo, hi_off, rho2_hi) = arc_extent(&fp.theta, r2, &m, d)?;
    Some(Reach {
        lo: lo.checked_add(lo_off)?,
        hi: hi.checked_add(hi_off)?,
        rho2_lo,
        rho2_hi,
    })
}

/// Per axis, the low and high ends of where a vertex's truth can be.
pub(crate) type CacheReach = ([f64; 3], [f64; 3]);

/// **Where a vertex's truth can be** — per axis, the interval its cache proves. The self-touch
/// sieve reads it for seam vertices, the gate for the corners of a face it asks to stand apart.
///
/// A `Realized` cache holds the truth within `coord ± bound`; the radius is widened by a few ulps
/// of both operands (and the smallest normal), since `coord − bound` computed in `f64` may round
/// inward by half an ulp of the result — which the widening outgrows on either side of a binade.
/// Widening only over-keeps. Any other variant proves nothing, and the answer is the whole line:
/// the sieve then keeps every candidate the vertex takes part in and leaves it to the exact test,
/// which is the only safe reading of "no bound". Measured: every seam point in the census and the
/// perf folds is `Realized`; the suite's 72 others (a history past the cost cap, a carrier with no
/// witness) are where this arm runs.
pub(crate) fn cache_reach(cache: &nacre_topo::PointCache) -> CacheReach {
    match *cache {
        nacre_topo::PointCache::Realized { coord, bound } => {
            let c = coord.as_array();
            let w: [f64; 3] = core::array::from_fn(|k| {
                bound[k].upper_f64() * (1.0 + 4.0 * f64::EPSILON) + round_slack(c[k])
            });
            (
                core::array::from_fn(|k| c[k] - w[k]),
                core::array::from_fn(|k| c[k] + w[k]),
            )
        }
        nacre_topo::PointCache::Ceiling { .. } | nacre_topo::PointCache::Unrealized { .. } => {
            ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3])
        }
    }
}

/// **A few ulps of `x`, and the smallest normal** — how far an `f64` computed near `x` may round
/// from the value it stands for. A bound read in `f64` is widened by it to hold the truth:
/// [`cache_reach`]'s radius, and the self-touch sieve's rational circle boxes.
pub(crate) fn round_slack(x: f64) -> f64 {
    x.abs() * (4.0 * f64::EPSILON) + f64::MIN_POSITIVE
}

/// **Does this cylinder leave only whole rims?** — every loop of every lateral face is one edge
/// that closes on itself (a whole circle, the spelling [`lateral_spans`] reads a rim by), and
/// every planar face whose outer loop is such a circle (a cap's disk, a blind bore's floor) holds
/// only such circles. Then the cylinder puts no ruling edge, no pierce corner and no arc on any
/// other face, and a class that cuts its lateral cuts only circles of it: on that class the
/// cylinder's solid is sectioned as if the bore were filled or the boss removed, and every
/// traced boundary still closes. A disk with a polygon in it would not — the disk's own trace
/// stops before its holes are read.
fn only_whole_rims(model: &Model, faces: &[FaceRow], surf: Handle<Surface>) -> bool {
    let whole = |lp: &nacre_topo::Loop| whole_rim(model, lp);
    let loops = |f: &nacre_topo::Face| {
        std::iter::once(&f.outer)
            .chain(f.inner.iter())
            .cloned()
            .collect::<Vec<_>>()
    };
    let mut rims = Vec::new();
    for row in faces {
        let FaceRow::Cylinder(cf) = row else { continue };
        if cf.surf != surf {
            continue;
        }
        let Some(fh) = cf.face else { return false };
        for lp in loops(model.face(fh)) {
            if !whole(&lp) {
                return false;
            }
            rims.push(lp.half_edges[0].edge);
        }
    }
    faces.iter().all(|row| {
        let FaceRow::Plane(fi) = row else { return true };
        let Some(fh) = fi.face else { return true };
        let f = model.face(fh);
        let disk = f.outer.half_edges.len() == 1 && rims.contains(&f.outer.half_edges[0].edge);
        !disk || loops(f).iter().all(whole)
    })
}

/// **Where a plane with no world statement stands against a cylinder** — what the gate's nameless
/// arm reads before anything else.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Stance {
    /// The whole cylinder lies strictly on this side of the plane, in
    /// [`crate::combinatorics::side_of`]'s convention.
    Beside(i8),
    /// The plane runs along the axis within the radius: it cuts the lateral in two rulings.
    Rulings,
    /// The plane touches the cylinder along one ruling, proved exactly before its turn.
    Grazes,
    /// Neither is proved.
    Unknown,
}

/// **Whether the whole cylinder stands on one side of a nameless class's plane, and if not, whether
/// the plane cuts it in rulings or touches it along one** ([`Stance`]). `Beside` is one question
/// with two witnesses, asked in this order, and both sides are minted where every side is
/// ([`WorkingPlane::misses`], the judge's filter, which needs only the class's witness triangle):
/// * every corner of the lateral faces' boxes ([`lateral_box`], stated exactly) strictly on one side
///   ([`WorkingCyl::clear_of`]) — every rim, ruling, chord, pierce point and arc of the cylinder is
///   then on that side;
/// * a plane along the axis — exact through the two motions (`Model::plane_cylinder_relation`; a
///   turn about the axis leaves the axis alone) — whose distance from the axis is beyond the radius
///   ([`Judge::plane_cuts`]): the whole infinite cylinder is on the side of any point of its axis.
///   This is what a box cannot say when the plane passes through the box beside the arc — a
///   bore's box at its corner, a slot end concave from the material.
///
/// Within the radius such a plane cuts the lateral in two rulings, the one shape the ruling road
/// draws. Asked before the class's faces are: a cylinder beside the class has nothing on it,
/// whatever else that cylinder records.
///
/// **A plane whose turn is about the cylinder's own axis is asked exactly, before the turn**
/// ([`Model::cylinder_before_chain`]): the turn leaves the cylinder in place, so the plane's name
/// stands to the cylinder carried back as the plane stands to it in the world, and one clearance
/// call answers three ways — beyond the radius (`Beside`, the second witness above with no interval
/// to straddle), within it (`Rulings`), and exactly on it (`Grazes` — the answer an interval never
/// gives, since it proves two quantities apart and never equal). Any other plane along the axis is
/// asked on the interval road.
fn nameless_stance(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    c: usize,
    cyl: &WorkingCyl,
) -> Stance {
    let wp = &jd.planes[c];
    if let Some(side) = cyl.clear_of(wp) {
        return Stance::Beside(side);
    }
    if model.plane_cylinder_relation(wp.surf, cyl.surf) != Some(nacre_exact::AxisRelation::Along) {
        return Stance::Unknown;
    }
    let beside_the_axis = || {
        let axis_point = nacre_judge::WitnessPoint::at(cyl.def.origin());
        wp.misses(std::slice::from_ref(&axis_point))
            .map_or(Stance::Unknown, Stance::Beside)
    };
    if let Some(before) = model.cylinder_before_chain(wp.surf, &cyl.def)
        && let Some(name) = model
            .surface_name(wp.surf)
            .and_then(|n| n.narrow().copied())
    {
        return match nacre_exact::point_plane_clearance_rat(&name, &before.origin(), before.r2()) {
            nacre_exact::Orient::Positive => beside_the_axis(),
            nacre_exact::Orient::Negative => Stance::Rulings,
            nacre_exact::Orient::Zero => Stance::Grazes,
        };
    }
    match jd.plane_cuts(c, &cyl.def) {
        Some(false) => beside_the_axis(),
        Some(true) => Stance::Rulings,
        None => Stance::Unknown,
    }
}

/// **Is the line `p ∩ x` off every lateral face of `def`?** — `p` a class with no world statement
/// along the axis that sections or grazes the cylinder, `x` another class recorded on the same
/// cylinder, `pieces` the cylinder's lateral faces ([`lateral_pieces`]). A named `x` not along the
/// axis meets `p` in a line that crosses the axis direction, which no ruling does. Otherwise both
/// run along the axis — a named `x` by its coefficients, a nameless one because it sections or
/// grazes the cylinder too — and the line is off when there is none, the two planes parallel
/// ([`Judge::planes_parallel`]: two classes of one table are proved not coplanar, so parallel ones
/// share no line), or when [`line_off_faces`] shows it off the faces from the two planes alone —
/// off the surface, or at a direction no face's arc holds. `false` where neither is shown, and for
/// a `Wide` name, which nothing here reads.
fn line_off_lateral(
    jd: &Judge<'_, WorkingPlane>,
    p: usize,
    x: usize,
    def: &nacre_topo::CylinderDef,
    pieces: &[Footprint],
) -> bool {
    if jd.planes[x].world.is_some() {
        let Some(w) = jd.world_rat(x) else {
            return false;
        };
        if nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], &def.dir()) != nacre_exact::Orient::Zero {
            return true;
        }
    }
    let meet = nacre_judge::cylinder::AxisLine::Meet(p, x);
    jd.planes_parallel(p, x) == Some(true)
        || line_off_faces(jd, def, None, Some(meet), pieces) == Some(true)
}

/// **Is the line of sectioned class `s` and grazing class `g` a ruling of cylinder `def`?** — asked
/// where both are rational: the two planes on one chain ([`Model::surface_motion`]), the cylinder
/// carried back before it ([`Model::cylinder_before_chain`]), their two names meeting on it
/// ([`nacre_exact::quad::CylinderMeet::OnRuling`]). The chain turns both planes and leaves the
/// cylinder in place, so the line it carries out is a ruling still. `false` where any step is not
/// shown.
fn graze_ruling(
    model: &Model,
    jd: &Judge<'_, WorkingPlane>,
    s: usize,
    g: usize,
    def: &nacre_topo::CylinderDef,
) -> bool {
    let (ss, sg) = (jd.planes[s].surf, jd.planes[g].surf);
    if model.surface_motion(ss) != model.surface_motion(sg) {
        return false;
    }
    let Some(before) = model.cylinder_before_chain(ss, def) else {
        return false;
    };
    let name = |p| model.surface_name(p).and_then(|n| n.narrow().copied());
    let (Some(ns), Some(ng)) = (name(ss), name(sg)) else {
        return false;
    };
    matches!(
        nacre_exact::quad::plane_plane_cylinder(
            &ns,
            &ng,
            &before.origin(),
            &before.dir(),
            before.r2()
        ),
        Some(nacre_exact::quad::CylinderMeet::OnRuling(_))
    )
}

/// **Do the faces on class `c` stand apart from every lateral box?** — every face's edges are
/// straight (a curved edge bulges past its corners) and the box of its corners' proven reaches
/// ([`cache_reach`]; a corner whose cache proves no bound makes the box the whole space) lies
/// strictly beyond each box on some axis. A planar face lies in the hull of its own corners, so
/// a face apart by its corners' box is apart. `false` when there are no lateral boxes.
fn faces_apart(
    model: &Model,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    c: usize,
    boxes: &[[[nacre_exact::Rat; 3]; 2]],
) -> bool {
    if boxes.is_empty() {
        return false;
    }
    faces.iter().zip(plane_ix).all(|(row, ix)| {
        let (FaceRow::Plane(fi), ClassIx::Plane(k)) = (row, ix) else {
            return true;
        };
        if *k != c {
            return true;
        }
        let Some(fh) = fi.face else { return false };
        let f = model.face(fh);
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
            for he in &lp.half_edges {
                let e = model.edge(he.edge);
                let straight = e
                    .surfaces
                    .iter()
                    .all(|&s| matches!(model.surface(s), Surface::Plane { .. }));
                if !straight {
                    return false;
                }
                for &v in &e.vertices {
                    let (rl, rh) = cache_reach(model.vertex_cache(v));
                    for k in 0..3 {
                        lo[k] = lo[k].min(rl[k]);
                        hi[k] = hi[k].max(rh[k]);
                    }
                }
            }
        }
        let beyond = |a: f64, b: nacre_exact::Rat, below: bool| {
            nacre_exact::Rat::try_from_f64(a).is_some_and(|a| if below { a < b } else { a > b })
        };
        boxes.iter().all(|[blo, bhi]| {
            (0..3).any(|k| beyond(hi[k], blo[k], true) || beyond(lo[k], bhi[k], false))
        })
    })
}

/// **A box the lateral face lies in** — `[low, high]` on each world axis, rational and outward.
///
/// The face's reach along each axis is [`lateral_reach`]'s, which knows the face's **arc**, not
/// only its circle: a fillet is a quarter of its cylinder, and the circle's box around it is twice
/// as wide on each axis — wide enough that a plane passing the quarter at a distance cuts the box
/// (a 37° rib beside the fillet plate's corner, worked by hand: the circle's box reached 1.4 past
/// the rib's wall plane where the quarter's box cleared it by 1.6). A reach end carries one radical,
/// `lo − √ρ²` or `hi + √ρ²`; it is rounded outward to a rational by [`sqrt_at_least`], so the box
/// can only be wider than the face, never narrower. `None` is an overflow or a face with no span.
pub(crate) fn lateral_box(
    def: &nacre_topo::CylinderDef,
    fp: &Footprint,
) -> Option<[[nacre_exact::Rat; 3]; 2]> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    let mut lo = [zero; 3];
    let mut hi = [zero; 3];
    for k in 0..3 {
        let mut d = [zero; 3];
        d[k] = Rat::from_int(1);
        let reach = lateral_reach(def, fp, &d)?;
        lo[k] = reach.lo.checked_sub(sqrt_at_least(reach.rho2_lo)?)?;
        hi[k] = reach.hi.checked_add(sqrt_at_least(reach.rho2_hi)?)?;
    }
    Some([lo, hi])
}

/// A rational `u ≥ √x`, near it — the `f64` root, then checked exactly (`u² ≥ x`) and widened
/// until the check holds. `None` for a negative `x` or an overflow.
pub(crate) fn sqrt_at_least(x: nacre_exact::Rat) -> Option<nacre_exact::Rat> {
    use nacre_exact::Rat;
    let zero = Rat::from_int(0);
    if x == zero {
        return Some(zero);
    }
    if x < zero {
        return None;
    }
    let mut u = Rat::try_from_f64(x.to_f64().sqrt())?;
    // One part in 2³⁰ above: the `f64` root is within an ulp, so one widening settles it, and the
    // loop is only the proof.
    let widen = Rat::new((1 << 30) + 1, 1 << 30)?;
    for _ in 0..8 {
        if u.checked_mul(u)? >= x {
            return Some(u);
        }
        u = u.checked_mul(widen)?;
    }
    None
}

impl WorkingCyl {
    /// **The side of `wp`'s plane every lateral face of this cylinder lies on**, strictly —
    /// [`WorkingPlane::misses`] over the corners of all the faces' boxes at once (one cylinder can
    /// hold two faces apart along its axis, and a plane between them has no one side). `None`
    /// when there are no boxes, a corner is not proved off the plane, or the corners straddle it.
    /// The gate's first witness for [`Silence::Beside`] ([`nameless_stance`]); a road past the gate
    /// reads the record ([`WorkingCyl::beside`]) instead.
    fn clear_of(&self, wp: &WorkingPlane) -> Option<i8> {
        if self.boxes.is_empty() {
            return None;
        }
        let corners: Vec<nacre_judge::WitnessPoint> = self
            .boxes
            .iter()
            .flat_map(|[lo, hi]| {
                (0..8).map(move |bits: usize| {
                    nacre_judge::WitnessPoint::at(std::array::from_fn(|k| {
                        if bits >> k & 1 == 0 { lo[k] } else { hi[k] }
                    }))
                })
            })
            .collect();
        wp.misses(&corners)
    }
}

/// **Do two reaches provably miss each other?** — one lies wholly beyond the other, at either end.
///
/// Each end is a rational base and a radical, so a gap has to clear **two** roots:
/// `b.lo − √b.rho2_lo > a.hi + √a.rho2_hi` is `g > √p + √q`, which
/// [`nacre_exact::exceeds_root_sum`] answers exactly. `None` is `Rat` overflow forming a gap.
///
/// **Closed against closed**: equality is one face's rim touching the other at a point, which is
/// not clear. The footprint's *span* is read **open** instead, so the disk arm in
/// [`Corner::reaches`] spells its own comparison rather than borrowing this one — the same algebra
/// with the boundary the other way, and folding them together would move one convention silently.
fn reaches_apart(a: &Reach, b: &Reach) -> Option<bool> {
    Some(
        nacre_exact::exceeds_root_sum(b.lo.checked_sub(a.hi)?, a.rho2_hi, b.rho2_lo)?
            || nacre_exact::exceeds_root_sum(a.lo.checked_sub(b.hi)?, b.rho2_hi, a.rho2_lo)?,
    )
}

/// Is the closed interval `[lo, hi]` (in the reach's own `d·p` units) disjoint from the reach? —
/// [`reaches_apart`] against the reach a rational interval **is**: no radical at either end. The
/// station a plane offers is such an interval of one point.
fn reach_clears(reach: &Reach, lo: nacre_exact::Rat, hi: nacre_exact::Rat) -> Option<bool> {
    let zero = nacre_exact::Rat::from_int(0);
    reaches_apart(
        &Reach {
            lo,
            hi,
            rho2_lo: zero,
            rho2_hi: zero,
        },
        reach,
    )
}

/// **Does every lateral face of this cylinder provably miss the plane `n·p + d = 0`?** — whatever
/// the plane's tilt: the oblique arm's whole face question, and the parallel arm's second. The
/// plane's station in `n·p` units is `−d`; a face clears when its reach along `n`
/// ([`lateral_reach`], which reads the face's arc) is disjoint from that one point. An empty span
/// list is "unusable" ([`lateral_spans`]'s doc), never clear; so is overflow.
fn laterals_clear_plane(
    def: &nacre_topo::CylinderDef,
    footprints: &[Footprint],
    coeffs: &[nacre_exact::Rat; 4],
) -> bool {
    if footprints.is_empty() {
        return false;
    }
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let Some(station) = nacre_exact::Rat::from_int(0).checked_sub(coeffs[3]) else {
        return false;
    };
    footprints.iter().all(|fp| {
        lateral_reach(def, fp, &n).and_then(|reach| reach_clears(&reach, station, station))
            == Some(true)
    })
}

/// **Does one class carry closed traces of cylinders no road nests?** — `row` is the class's
/// answer about each cylinder ([`class_meet`]) beside that cylinder's axis. Two closed traces on
/// one class are nested by [`crate::nesting::disk_in_disk`], which answers for parallel axes
/// alone. On a perpendicular class every circle is of a cylinder whose axis is the class's normal,
/// so they are all parallel; an oblique section beside a circle, or beside the section of another
/// axis, is a pair no road nests. A class across a cylinder's axis is counted by that relation
/// alone, whether or not the cylinder reaches the class — wider than the pairs that meet, never
/// narrower. Parallelism is transitive, so neighbours in the row suffice.
fn carries_unnested(row: &[([nacre_exact::Rat; 3], &ClassMeet)]) -> bool {
    let carried: Vec<&[nacre_exact::Rat; 3]> = row
        .iter()
        .filter(|(_, meet)| matches!(meet, ClassMeet::Sections | ClassMeet::Across))
        .map(|(dir, _)| dir)
        .collect();
    row.iter().any(|(_, meet)| **meet == ClassMeet::Sections)
        && carried
            .windows(2)
            .any(|w| !nacre_exact::parallel_rat(w[0], w[1]))
}

/// **How a plane's section stands to a lateral face's boundary** — [`boundary_stand`]'s answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BoundaryStand {
    /// Every boundary section lies apart from the plane's section, all on one side of it: the
    /// section misses the face's band.
    Clear,
    /// The section crosses a boundary section, is one, or lies between two: it reaches the face's
    /// band.
    Sections,
    /// The section touches a boundary section at one point — a point the two share with no
    /// crossing to name it by.
    Touch,
}

/// **Where the plane `coeffs` stands against a lateral face's boundary sections `boundary`**
/// ([`crate::planes::lateral_boundary_sections`]) — each pair asked the one question
/// ([`nacre_exact::quad::section_meet`]; for two apart, which is higher,
/// [`nacre_exact::quad::section_height_order`] — two circles compare as their heights). The face
/// lies between its boundary sections, so a section apart from all of them and on one side of all
/// misses it; one crossing, coinciding with or between them reaches it; a touch anywhere is
/// [`BoundaryStand::Touch`], which takes precedence. The gate's ([`oblique_section_meet`], and the
/// touch it refuses across the axis too) and the tracer's (`arrangement::trace_on_face`), so the
/// two cannot part. `None` is a plane along the axis among them.
pub(crate) fn boundary_stand(
    def: &nacre_topo::CylinderDef,
    boundary: &[[nacre_exact::Rat; 4]],
    coeffs: &[nacre_exact::Rat; 4],
) -> Option<BoundaryStand> {
    use nacre_exact::quad::{
        SectionMeet, SectionOrder, Toward, section_height_order, section_meet,
    };
    let mut side: Option<core::cmp::Ordering> = None;
    let mut reaches = false;
    for w in boundary {
        match section_meet(coeffs, w, def)? {
            SectionMeet::Touch => return Some(BoundaryStand::Touch),
            SectionMeet::Same | SectionMeet::Cross => reaches = true,
            SectionMeet::Apart => {
                // Apart, so ordered with no generator.
                let SectionOrder::Ordered(o) =
                    section_height_order(coeffs, w, def, None, Toward::At)?
                else {
                    return None;
                };
                reaches |= side.is_some_and(|s| s != o);
                side = Some(o);
            }
        }
    }
    Some(if reaches {
        BoundaryStand::Sections
    } else {
        BoundaryStand::Clear
    })
}

/// **Is this loop a whole circle** — one edge that closes on itself, the spelling a rim is read by
/// ([`lateral_spans`])?
fn whole_rim(model: &Model, lp: &nacre_topo::Loop) -> bool {
    lp.half_edges.len() == 1 && {
        let [a, b] = model.edge(lp.half_edges[0].edge).vertices;
        a == b
    }
}

/// **Does oblique class `c` section this cylinder on its own arrangement?** —
/// [`ClassMeet::Sections`]' proof: the plane touches no section a lateral face's boundary rides — a
/// rim's or an arc's, a circle's or an ellipse's ([`boundary_stand`] over
/// [`crate::planes::lateral_boundary_sections`] — a touch is a point the section and that boundary's
/// plane share with no crossing to name it by). Then the class draws the section — a whole ellipse,
/// or arcs of it cut by rims, rulings and the class's other segments — and the tracer reads which
/// from the faces' boundaries. Where the class's faces reach the lateral the section is a line the
/// lateral's chart cuts its cells at like any other; a result face bounded along it carries an
/// elliptical edge. A class whose plane is a boundary section of the face — an operand's own slanted
/// cap — reaches it too.
/// The touch test asks each boundary section whole, also where the face holds only an arc of it (a
/// fillet's quarter): wider than the faces, never narrower. `Some(Clear)` when every face's boundary
/// is missed, `Some(Sections)` when one is reached, `None` when the proof fails — a cylinder with no
/// lateral face row, a face with no readable boundary, a touch.
fn oblique_section_meet(
    model: &Model,
    faces: &[FaceRow],
    cyl: &WorkingCyl,
    fps: &[Footprint],
    coeffs: &[nacre_exact::Rat; 4],
) -> Option<ClassMeet> {
    // No lateral face to read is not a clear one (`lateral_spans`' «unusable»).
    if fps.is_empty() {
        return None;
    }
    let mut sectioned = false;
    for row in faces {
        let FaceRow::Cylinder(cf) = row else { continue };
        if cf.surf != cyl.surf {
            continue;
        }
        let boundary = lateral_boundary_sections(model, model.face(cf.face?), &cyl.def)?;
        match boundary_stand(&cyl.def, &boundary, coeffs)? {
            BoundaryStand::Touch => return None,
            BoundaryStand::Sections => sectioned = true,
            BoundaryStand::Clear => {}
        }
    }
    Some(if sectioned {
        ClassMeet::Sections
    } else {
        ClassMeet::Clear
    })
}

/// **Does a plane across the axis touch a boundary section of cylinder `cyl`'s lateral faces at a
/// point?** — the across arm's half of the touch [`oblique_section_meet`] refuses: a slanted rim's
/// highest or lowest point on the plane. Two planes across the axis never touch (their sections are
/// one or apart), so a face bounded by circles alone answers `false`; a face whose boundary cannot
/// be read is left to the tracer, which reads the same boundary and declines by its own name.
fn touches_a_boundary(
    model: &Model,
    faces: &[FaceRow],
    cyl: &WorkingCyl,
    coeffs: &[nacre_exact::Rat; 4],
) -> bool {
    faces.iter().any(|row| {
        let FaceRow::Cylinder(cf) = row else {
            return false;
        };
        cf.surf == cyl.surf
            && cf.face.is_some_and(|f| {
                lateral_boundary_sections(model, model.face(f), &cyl.def)
                    .and_then(|b| boundary_stand(&cyl.def, &b, coeffs))
                    == Some(BoundaryStand::Touch)
            })
    })
}

/// **Two statements of one cylinder surface** — parallel axes on one line, one radius — under two
/// handles: a `translate`d twin (its motion key differs) or a statement with another `ref_dir`.
/// Cylinders intern by their literal statement, so the class table holds both; the arrangement has
/// no name for two classes on one surface (their circles would coincide on every ⊥ class), so the
/// pair rule refuses them as the coincident pair they are. Overflow answers `true` — "not shown
/// distinct" refuses, it never lets a pair through.
fn same_surface(a: &nacre_topo::CylinderDef, b: &nacre_topo::CylinderDef) -> bool {
    if !nacre_exact::parallel_rat(&a.dir(), &b.dir()) || a.r2() != b.r2() {
        return false;
    }
    let (oa, ob) = (a.origin(), b.origin());
    let d: Option<[nacre_exact::Rat; 3]> = (|| {
        Some([
            ob[0].checked_sub(oa[0])?,
            ob[1].checked_sub(oa[1])?,
            ob[2].checked_sub(oa[2])?,
        ])
    })();
    let Some(d) = d else { return true };
    let zero = nacre_exact::Rat::from_int(0);
    match nacre_exact::cross3_rat(&d, &a.dir()) {
        Some(c) => c.iter().all(|x| *x == zero),
        None => true,
    }
}

/// **The face-level clearance for a cylinder pair**: every lateral face of `a` against every
/// lateral face of `b`, each pair asked by [`pieces_clear`]; the classes clear when **every** pair
/// does. Then no point of `b`'s faces lies on `a`'s, so the two classes share no face — the
/// proposition the arrangement needs, which the surface distance
/// ([`nacre_exact::cylinders_clear`]) is only one sufficient condition for.
fn lateral_faces_clear(faces: &[FaceRow], a: &WorkingCyl, b: &WorkingCyl) -> bool {
    let (fa, fb) = (lateral_pieces(faces, a.surf), lateral_pieces(faces, b.surf));
    fa.iter()
        .all(|x| fb.iter().all(|y| pieces_clear(&a.def, x, &b.def, y)))
}

/// **A cylinder class's lateral faces as the pieces a clearance asks about** — one footprint per
/// face, or, where [`lateral_footprints`] is empty («unusable», [`lateral_spans`]' doc), one
/// footprint of no span and the whole circle: the class still has faces, and that piece reads as
/// unbounded along anything not perpendicular to the axis. So «every piece is apart» is never an
/// empty truth. The door is for the readers that ask «apart» of the class's faces; the readers that
/// decline on an unusable list (the oblique arm, `faces_apart`, the cut record) read
/// [`lateral_footprints`] and keep that. Built where it is asked, not with the class — the gate's
/// own footprints are lazy for the reason [`wall_faces_clear`] gives.
pub(crate) fn lateral_pieces(faces: &[FaceRow], surf: Handle<Surface>) -> Vec<Footprint> {
    let v = lateral_footprints(faces, surf);
    if v.is_empty() {
        vec![Footprint {
            span: None,
            theta: Theta::Whole,
        }]
    } else {
        v
    }
}

/// **Does a line along cylinder `def`'s axis miss every lateral face of it?** — one predicate,
/// asked by the split of a segment's line that lies on the surface or that the intervals cannot
/// tell from one there, and by the gate of the line two recorded classes meet in
/// ([`line_off_lateral`]): such a line crosses no ruling (they are parallel) and no circle (its
/// class is ∥ the axis), so it touches the lateral only where a face holds its angle.
///
/// **The faces are asked, not the surface**, and only the line's angle: a face's rulings and arcs
/// stand at angles its arc holds, and two lines along the axis on one cylinder at different angles
/// never meet. The span is not read — the angle alone is what the population needs.
///
/// Two bodies. Holding the exact line, the half-space along its own direction ([`Theta::at_line`]
/// as a piece `span: None` × its angle, against each face — [`pieces_clear`]): `Some(true)` off
/// every face, `Some(false)` on one, `None` where the angle has no rational spelling. Holding only a
/// point of the line, or the two planes it is the meet of ([`AxisLine`]), its direction off the
/// axis against each face's arc ([`Judge::axis_line_off_arcs`] — a direction is the same along the
/// whole line, so any point of it will do; the meet of two planes is off as well when it is off the
/// surface): `Some(true)` where that proves the line off, `None` otherwise. A whole lateral
/// (`Theta::Whole`) is the whole turn, so a line of two planes is off it only off the surface, and a
/// face at one direction (`Theta::At`) is not shown clear. Where both bodies can answer, a
/// debug build checks the interval one never proves clear what the exact one does not.
///
/// [`AxisLine`]: nacre_judge::cylinder::AxisLine
pub(crate) fn line_off_faces(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    line: Option<&nacre_exact::quad::MeetLine>,
    interval_line: Option<nacre_judge::cylinder::AxisLine<'_>>,
    pieces: &[Footprint],
) -> Option<bool> {
    use nacre_exact::Rat;
    let interval = || -> bool {
        let Some(at) = interval_line else {
            return false;
        };
        let arcs: Option<Vec<[[Rat; 3]; 2]>> = pieces
            .iter()
            .map(|fp| match fp.theta {
                Theta::Arc(a) => Some([a.from, a.to]),
                // A whole lateral holds every direction: two ends pointing one way are the whole
                // turn whichever way that is, so `x` need not be radial, and its dot product with
                // itself does not overflow.
                Theta::Whole => {
                    let x = [Rat::from_int(1), Rat::from_int(0), Rat::from_int(0)];
                    Some([x, x])
                }
                Theta::At(_) => None,
            })
            .collect();
        arcs.is_some_and(|arcs| jd.axis_line_off_arcs(def, at, &arcs))
    };
    let Some(line) = line else {
        return interval().then_some(true);
    };
    let theta = Theta::at_line(&line.base(), &def.origin(), &def.dir())?;
    let seg = Footprint { span: None, theta };
    let clear = pieces.iter().all(|fp| pieces_clear(def, &seg, def, fp));
    debug_assert!(
        clear || !interval(),
        "the interval body proves a line off the faces the exact one puts it on"
    );
    Some(clear)
}

/// **Do piece `x` of cylinder `a`'s lateral and piece `y` of cylinder `b`'s share no point?** — the
/// one question asked of two pieces of cylinder laterals: a face, an edge (an arc is `[t, t]` × its
/// arc, a ruling its axis span × the whole circle), a line on a cylinder (`span: None` × its
/// angle). Who asks is design 「인구 게이트」's to list. A piece is a [`Footprint`]: its own extent, or —
/// where the extent cannot be stated — `span: None` (unbounded along the axis) or [`Theta::Whole`]
/// (the whole circle), which are the carrier's own extent on that axis and prove only what the
/// carrier would. No caller can hand it a carrier in place of a piece it has.
///
/// `true` is a proof, by a set `S` with `x ⊂ S` and `y ∩ S = ∅`; three kinds of `S` answer, and
/// any one suffices:
/// 1. **A half-space** — the two reaches part along one of [`separating_dirs`] ([`separated`]). The
///    two axes are one question, not two: [`separated`] at `d = m_a` and `d = m_b` asks one
///    piece's **span** against the other's **reach**, both ways round — along its own axis a
///    piece's reach *is* that span (`d ∥ m` leaves no radial term), as [`lateral_reach`] returns.
/// 2. **The other cylinder's open solid** — every corner of one piece's box ([`lateral_box`])
///    strictly inside the other's infinite cylinder. The box is convex and so is the open solid,
///    so the whole box — and the piece in it — is off the other's surface, and so off every piece
///    of it. This is what a piece lying across the inside of another needs: a pin nested crosswise
///    in a larger cylinder lies inside the larger lateral's convex hull, so no plane at all
///    separates the two. The face-level twin of [`nacre_exact::cylinders_nested`], the parallel
///    surfaces' containment. Strictly: a corner on the surface may be a piece point on it
///    — a rim touching the other's lateral from inside.
/// 3. **For parallel axes, the common cross-section** — the arcs miss ([`cross_sections_clear`]).
///
/// Sound for any piece inside the footprint: every reach and box is an enclosure of the footprint,
/// wider than the piece, never narrower. `None` along a direction, a box that cannot be stated, a
/// corner on the surface — each is "not shown", never clear.
pub(crate) fn pieces_clear(
    a: &nacre_topo::CylinderDef,
    x: &Footprint,
    b: &nacre_topo::CylinderDef,
    y: &Footprint,
) -> bool {
    separating_dirs(a, x, b, y)
        .iter()
        .any(|d| separated(a, x, b, y, d) == Some(true))
        || piece_inside(a, x, b)
        || piece_inside(b, y, a)
        || (nacre_exact::parallel_rat(&a.dir(), &b.dir())
            && cross_sections_clear(a, x, b, y) == Some(true))
}

/// **Is a line on cylinder `def` apart from a face of it?** — the half-space along the line's own
/// direction (its [`Theta::At`], [`separating_dirs`]' candidate), read three ways: `Some(true)` the
/// line is off the face's arc, `Some(false)` on it or on an end, `None` arithmetic. Against a face of
/// its own statement this witness is the whole of [`pieces_clear`]'s answer — the axis directions
/// read a span the line has not, and no face's box lies inside its own cylinder — so the readers
/// that must tell «on the face» from «could not say» (the [`SharedRuling`] record, the tangency
/// rows) ask it here, where [`pieces_clear`] keeps only its `Some(true)`. `None` too for a line
/// with no single direction.
fn line_apart(def: &nacre_topo::CylinderDef, line: &Footprint, face: &Footprint) -> Option<bool> {
    let Theta::At(v) = line.theta else {
        return None;
    };
    separated(def, line, def, face, &v)
}

/// **Is piece `fp` of cylinder `p`'s lateral strictly inside cylinder `q`'s open solid?** —
/// [`pieces_clear`]'s second witness: every corner of the piece's box ([`lateral_box`]) strictly
/// within `q`'s radius. `false` when the box cannot be stated.
fn piece_inside(p: &nacre_topo::CylinderDef, fp: &Footprint, q: &nacre_topo::CylinderDef) -> bool {
    lateral_box(p, fp).is_some_and(|[lo, hi]| {
        (0..8usize).all(|bits| {
            let c: [nacre_exact::Rat; 3] =
                std::array::from_fn(|k| if bits >> k & 1 == 0 { lo[k] } else { hi[k] });
            nacre_exact::cylinder_radial_side(&c, &q.origin(), &q.dir(), q.r2())
                == nacre_exact::Orient::Negative
        })
    })
}

/// **The separating directions a cylinder pair can state rationally** — and why, for two pieces of
/// faces, these are all of them.
///
/// What a footprint describes is a box in `(axis parameter × angle)`, and the boundary of such a
/// face is made of **cap planes** (normal = the axis), **rulings** (direction = the axis) and
/// **arcs**. A separating axis is a face normal or an edge×edge, so the rational candidates are
/// each axis and — for skew axes — the two rulings' cross product, the common perpendicular. What
/// is left out is the arcs' radial continuum, which no rational direction names; the tool for that
/// is [`cross_sections_clear`], and it needs the common cross-section chart that only **parallel**
/// axes have. A piece inside the other's convex hull has no separating plane at all, rational or
/// not; that is [`piece_inside`]'s, the other witness of [`pieces_clear`]. So a pair none of these
/// clears is refused as "not shown to clear", which is true, and completing the test further means
/// the same test on a richer candidate set.
///
/// **The third direction is the face-level twin of a rung above it.** With a whole circle and no
/// span, `d ⊥ m_a` and `d ⊥ m_b` make both reaches `d·o ± r|d|`, so "apart" reads
/// `|d·(o_b − o_a)| > (r_a + r_b)|d|` — [`nacre_exact::cylinders_clear`]'s skew branch, letter for
/// letter.
///
/// Parallel axes have no third direction (the cross product is zero, and the surface rung is the
/// distance rule with its own parallel branch); overflow forming it simply leaves the list short.
///
/// **A piece of one direction names one more: its own radial vector `v`** ([`Theta::At`], a line
/// along the axis on the cylinder). Along `d = v` — perpendicular to the axis, so no span is read —
/// the line reaches the single value `d·o + |v|²`, and a face of the same cylinder whose arc does
/// not hold `v` reaches only `d·o + max(v·from, v·to)` there, which is less: the arc's peak along
/// `v` is at an end, and both ends are other directions of the same radius. So the half-space
/// shows a line off a face's arc apart from the face, and a line on an end ruling (`v = from`)
/// touching it — each reach read about its own statement's axis, whichever way that runs.
fn separating_dirs(
    a: &nacre_topo::CylinderDef,
    x: &Footprint,
    b: &nacre_topo::CylinderDef,
    y: &Footprint,
) -> Vec<[nacre_exact::Rat; 3]> {
    let (ma, mb) = (a.dir(), b.dir());
    let mut out = vec![ma, mb];
    if !nacre_exact::parallel_rat(&ma, &mb) {
        if let Some(perp) = nacre_exact::cross3_rat(&ma, &mb) {
            out.push(perp);
        }
    }
    for fp in [x, y] {
        if let Theta::At(v) = fp.theta {
            out.push(v);
        }
    }
    out
}

/// **Do face `x` of `a` and face `y` of `b` provably miss each other along `d`?** — each face's
/// reach along the one direction ([`lateral_reach`]), asked whether they are disjoint. Sound for
/// any `d`: a reach is a bounding interval of the face's projection, and two sets whose
/// projections miss cannot share a point. `None` is "not proved" — an unstatable span, or
/// overflow.
fn separated(
    a: &nacre_topo::CylinderDef,
    x: &Footprint,
    b: &nacre_topo::CylinderDef,
    y: &Footprint,
    d: &[nacre_exact::Rat; 3],
) -> Option<bool> {
    reaches_apart(&lateral_reach(a, x, d)?, &lateral_reach(b, y, d)?)
}

/// **Two lateral faces on parallel axes: do their rims' arcs miss each other in the common
/// cross-section?**. The chart is `a`'s: `û₁` the reference direction's unit part ⊥
/// the axis, `û₂ = m̂ × û₁` — rational exactly when both norms are (`inv_sqrt_rat`; every prism
/// on a world frame), else `None`. Each face's arc is its angular extent, or the whole circle
/// when none is stated; the question is [`nacre_geom::mixed::arcs_share_a_point`]'s, and a
/// touch is not clear. Sound for any face of the class because the extent is the footprint's
/// bounding arc: wider than a notched face, never narrower.
fn cross_sections_clear(
    a: &nacre_topo::CylinderDef,
    x: &Footprint,
    b: &nacre_topo::CylinderDef,
    y: &Footprint,
) -> Option<bool> {
    use nacre_exact::Rat;
    use nacre_geom::mixed::ArcSpec;
    let zero = Rat::from_int(0);
    let (o, m, e) = (a.origin(), a.dir(), a.ref_dir());
    let (mm, em) = (
        nacre_exact::dot3_rat(&m, &m)?,
        nacre_exact::dot3_rat(&e, &m)?,
    );
    let mut e1 = [zero; 3];
    for k in 0..3 {
        e1[k] = mm.checked_mul(e[k])?.checked_sub(em.checked_mul(m[k])?)?;
    }
    let inv_e1 = nacre_exact::inv_sqrt_rat(nacre_exact::dot3_rat(&e1, &e1)?)?;
    let inv_m = nacre_exact::inv_sqrt_rat(mm)?;
    let scaled = |v: &[Rat; 3], k: Rat| -> Option<[Rat; 3]> {
        Some([
            v[0].checked_mul(k)?,
            v[1].checked_mul(k)?,
            v[2].checked_mul(k)?,
        ])
    };
    let u1 = scaled(&e1, inv_e1)?;
    let u2 = scaled(
        &nacre_exact::cross3_rat(&m, &e1)?,
        inv_m.checked_mul(inv_e1)?,
    )?;
    let chart = |p: &[Rat; 3]| -> Option<[Rat; 2]> {
        let mut v = *p;
        for k in 0..3 {
            v[k] = v[k].checked_sub(o[k])?;
        }
        Some([
            nacre_exact::dot3_rat(&v, &u1)?,
            nacre_exact::dot3_rat(&v, &u2)?,
        ])
    };
    let add2 = |p: [Rat; 2], q: [Rat; 2]| -> Option<[Rat; 2]> {
        Some([p[0].checked_add(q[0])?, p[1].checked_add(q[1])?])
    };
    let spec = |c: &nacre_topo::CylinderDef, fp: &Footprint| -> Option<ArcSpec> {
        let centre = chart(&c.origin())?;
        // The 2-D arc court works in `Rat`; a wide square is not judged here (`None`).
        let r2 = c.r2().narrow()?;
        let (start, end) = match &fp.theta {
            // One direction is no arc of this court ([`ArcSpec`]'s `start == end` is the whole
            // circle); it names its own separating direction ([`separating_dirs`]).
            Theta::At(_) => return None,
            Theta::Arc(arc) => (
                add2(
                    centre,
                    [
                        nacre_exact::dot3_rat(&arc.from, &u1)?,
                        nacre_exact::dot3_rat(&arc.from, &u2)?,
                    ],
                )?,
                add2(
                    centre,
                    [
                        nacre_exact::dot3_rat(&arc.to, &u1)?,
                        nacre_exact::dot3_rat(&arc.to, &u2)?,
                    ],
                )?,
            ),
            Theta::Whole => {
                // A whole circle's seam, `centre + (r, 0)`, is a rational point only for a
                // rational radius; without one this pair is not judged here (`None`, the
                // caller's honest road).
                let p = add2(centre, [c.radius_rat()?, zero])?;
                (p, p)
            }
        };
        Some(ArcSpec {
            centre,
            r2,
            start,
            end,
        })
    };
    nacre_geom::mixed::arcs_share_a_point(&spec(a, x)?, &spec(b, y)?).map(|meet| !meet)
}

#[cfg(test)]
#[path = "../tests/cyl_gate.rs"]
mod tests;
