use super::*;
/// **The population gate** — decides, exactly, whether this operand pair stays inside
/// the axis-perpendicular population the cylinder arrangement serves, and names the refusal
/// otherwise. All arithmetic is checked `Rat` on world-stated descriptions; anything the gate
/// cannot decide exactly is [`RejectReason::CylinderGateUndecided`] — a conservative honest
/// refusal, never a guess.
///
/// Per (plane class, cylinder) pair, with `n` the class's rational normal and `m`/`o`/`r` the
/// cylinder's raw axis/origin/radius:
/// - `n × m = 0` — a perpendicular cut. Passes, a cap seated flush on the other body included.
/// - `n · m = 0` — a wall parallel to the axis. It must provably miss the **rectangle** the
///   cylinder occupies in that plane. The infinite plane clearing the axis by more than `r`
///   (`(n·o + d)² > r²·|n|²`) settles it outright and decides most inputs; otherwise each face on
///   the class answers for itself, across the strip or along a lateral face's span
///   ([`face_clears_footprint`]). A face not shown to miss either rides the rulings road — the
///   wall's plane within the radius (`0 ≤ d < r`, through the axis or offset from it), the
///   pair **recorded and passed** — or, the plane exactly `r` from the axis, **passes and is
///   recorded as a [`Tangency`] instead**: a tangency divides nothing, so it earns no crossing.
///
///   ★★★★★ **One clearance call, three records, and that is the whole rule.** `Positive` clears
///   and writes nothing; `Negative` writes a crossing (two rulings); `Zero` writes a tangency (one
///   grazing line). The third arm does not refuse, and it does not `crossings.insert` either:
///   that would put the
///   pair on the ruling road, and three roads there spell "two distinct roots". Not recording it
///   keeps every one of those sentences true and the arrangement unchanged: ☑ of 21 measured
///   cells **15 assemble**, `validate` clean and volumes exact; the other 6 are the third-plane
///   population [`Tangency::line_in_another_plane`] names, which the arrangement refuses on its
///   own (`CoincidentNodes`).
/// - anything else — an oblique plane. It passes when every lateral face of the cylinder
///   provably misses the plane — the face's reach along `n` against the plane's station
///   ([`lateral_reach`], [`oblique_plane_clears`]) — and is otherwise recorded and refused as
///   [`RejectReason::ObliqueCylinderCut`] (an ellipse). ★ Every one of the
///   gate's four sites speaks about faces; a gusset whose slanted plane runs past a plate's holes
///   is not refused here for the plane alone.
///
/// Per cylinder pair **of different owners** (two classes of one valid solid keep their faces
/// apart by construction and are not asked): axes clear of each other (`dist > r₁+r₂`, whatever
/// their orientation) pass outright; one surface stated under two handles ([`same_surface`]) is
/// refused as the coincident pair it is; otherwise the pair passes when the **faces** of one class
/// provably miss the other's — every lateral face's axis span against the other faces' reach
/// along that axis ([`lateral_faces_clear`], either direction; for parallel axes one cylinder
/// strictly inside the other clears outright, [`nacre_exact::cylinders_nested`], and otherwise
/// the spans alone decide) — and a pair whose faces cannot be shown to miss is recorded and refused as
/// [`RejectReason::CylinderPairContact`] (M6b, where the quartic intersection curve lives).
#[allow(clippy::type_complexity)]
pub(crate) fn cylinder_gate(
    model: &Model,
    cyl_surfs: &[Handle<Surface>],
    geom: &[WorkingPlane],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    n_a: usize,
) -> Result<
    (
        Vec<WorkingCyl>,
        std::collections::HashSet<(usize, usize)>,
        Vec<Tangency>,
    ),
    BoolError,
> {
    // ★ Every question below is a **sign**, and the scalar layer answers signs totally: the
    // local checked-`Rat` closures this used to carry declined on overflow, which put a width
    // limit inside `CylinderGateUndecided` and made that name say less than it claimed.
    use nacre_exact::Orient;
    let undecided = || reject(RejectReason::CylinderGateUndecided);

    let mut crossings = std::collections::HashSet::new();
    let mut tangencies: Vec<Tangency> = Vec::new();
    // ★ **Two more records, read once each**: the (plane, cylinder) pairs the oblique
    // arm could not show apart, and the cylinder pairs the pair rule could not. Today their one
    // reader is the refusal below each loop; the day the ellipse road and the
    // cylinder–cylinder road (M6b) arrive, that reader becomes a hand-over — the graduation
    // `crossings` and `tangencies` already made from a refusal to a record.
    // They stay local: nothing downstream reads them yet, and a field no one reads is not built.
    let mut oblique: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let mut cyls = Vec::with_capacity(cyl_surfs.len());
    for &surf in cyl_surfs {
        // ★ **The world statement or nothing.** A moved cylinder's def is written before its
        // motion, and every test below compares it against world planes. A chain that is a pure
        // rational translation carries it out exactly ([`world_cylinder_def`] — the same door the
        // face rows take); anything else (a rotation, a frame node, overflow) has no world
        // description here and is refused rather than measured across two frames.
        let Some(def) = world_cylinder_def(model, surf) else {
            return Err(undecided());
        };
        let nacre_geom::Surface::Cylinder(cache) = model.surface_cache(surf) else {
            unreachable!("push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache")
        };
        // ★ **One surface, two solids.** Cylinders intern by their exact statement, so two operands
        // whose laterals coincide arrive as one class carrying rows of both — the coaxial pair of
        // equal radius, by another spelling. The chart reads a class as one solid's surface
        // (`cyl_chart::chart_of`), so this is refused here, by the name the pair rule gives that
        // pair, instead of asserting inside the chart: two identical
        // circle prisms fused would reach that assertion. Otherwise the side
        // that owns the class is written on it — the pair loop reads it.
        let owned = |rows: &[FaceRow]| {
            rows.iter()
                .any(|r| matches!(r, FaceRow::Cylinder(cf) if cf.surf == surf))
        };
        let owner = match (owned(&faces[..n_a]), owned(&faces[n_a..])) {
            (true, true) => return Err(reject(RejectReason::CylinderPairContact)),
            (true, false) => SolidSide::A,
            (false, true) => SolidSide::B,
            // A class is named by some face row; a table with none describes nothing.
            (false, false) => return Err(undecided()),
        };
        cyls.push(WorkingCyl {
            surf,
            def,
            realized: *cache,
            owner,
        });
    }

    // ★ **One row per `cyl_surfs` entry, in order** — the loop above either pushes or returns, so
    // this table is indexed by the very `ClassIx::Cyl` number that named the surface. Consumers
    // index it directly (`merge_circles`), which is only sound while that holds.
    debug_assert_eq!(
        cyls.len(),
        cyl_surfs.len(),
        "the cylinder class table is index-aligned with the class numbering"
    );

    for (ci, cyl) in cyls.iter().enumerate() {
        let (o, m, r2) = (cyl.def.origin(), cyl.def.dir(), cyl.def.r2());
        // The footprint's second axis, gathered **lazily and at most once** per cylinder: it does
        // not depend on the plane class the loop below walks, but almost no boolean ever asks for
        // it — the plane-level test decides first. ★ Measured before this was made lazy: one cut
        // over a plate with 16 bores built the table 16 times and read it 0, which is exactly the
        // shape `wall_faces_clear` warns about two doc comments below.
        let mut footprints: Option<Vec<Footprint>> = None;
        for (c, wp) in geom.iter().enumerate() {
            // ★ **A world description or nothing** — the question this loop asks is geometric
            // (does this class's plane clear that cylinder), and both sides have to speak about
            // the world. `rotated` is not that question: a plane whose truth carries a
            // translation is *judged* through its chain and still has exact world coefficients,
            // and that is the population the rulings road serves.
            let Some(coeffs) = wp.world_rat else {
                return Err(undecided());
            };
            let n = [coeffs[0], coeffs[1], coeffs[2]];
            // A perpendicular cut is the circle population, and it passes — **including a cap
            // seated flush on the other body's face**. That seating used to be refused, and the
            // refusal was wider than anything it could name: what makes a seated circle hard is
            // its boundary meeting the counterpart's boundary, and a boundary is either an edge
            // on a plane (that plane is parallel to the axis → the wall rule below, or oblique →
            // the oblique rule; both are already conservative because the wall rule judges the
            // *infinite* plane, not the face) or another cylinder's rim (two circles can only
            // overlap when the axes stand closer than r₁+r₂ → the pair rule below). So what the
            // seated rule turned away was exactly the population whose circle lies wholly inside
            // the counterpart's face — measured across the family (through hole, blind hole, boss
            // fuse, common, drilling a pocket floor): every one exact and validating clean.
            //
            // ★ **Nothing downstream has to catch a degenerate seating, because this gate still
            // does** — by the two rules below rather than by a rule about seating. A seated
            // circle can only reach the counterpart's boundary through a plane parallel to the
            // axis (which must clear the radius, cross it, or touch it) or an oblique one (refused
            // outright) or another cylinder's rim (the pair rule), so a tangency or a crossing is
            // **named** before the arrangement ever sees it. ★ Named, not refused: since the
            // tangent arm opened, a touch passes with a [`Tangency`] row and the *verdict* is
            // `assembly::tangency_reject`'s. What the sentence guarantees is unchanged — no
            // degenerate seating reaches the arrangement unnamed.
            if !nacre_exact::parallel_rat(&n, &m) {
                if nacre_exact::dot_sign_rat(&n, &m) != Orient::Zero {
                    // ★ **The oblique arm asks the faces** — the fourth of the gate's
                    // four sites to speak about faces rather than surfaces. The plane's station
                    // `n·p = −d` against every lateral face's reach along `n` ([`lateral_reach`],
                    // the question its doc was written for): if every face provably misses the
                    // infinite plane, no face of that plane's class can meet the lateral, and the
                    // caps are the plane–plane arrangement's business. A pair not shown to miss is
                    // recorded; the refusal reads the record once, after this loop.
                    let fps = footprints.get_or_insert_with(|| lateral_footprints(faces, cyl.surf));
                    if !oblique_plane_clears(&cyl.def, fps, &coeffs) {
                        oblique.insert((c, ci));
                    }
                    continue;
                }
                // A parallel wall must provably miss the lateral surface. ★ **The question is
                // about the wall's *faces*, not its plane** — the uniform-slab theorem this feeds
                // says so in its own words ("the other operand's **boundary** does not meet the
                // open cylinder slab"). Judging the infinite plane is a cheaper *sufficient*
                // condition, so it is asked first and still decides most inputs; when it fails,
                // the faces on this class get to answer for themselves. Refusing on the plane
                // alone turned away a whole family the engine serves — a boss standing far away
                // whose wall plane, extended, happens to pass through a hole.
                //
                // ★★ What a face is asked is whether it misses the **rectangle** this cylinder
                // occupies in that plane: the strip across, the lateral face's span along. See
                // [`face_clears_footprint`].
                if nacre_exact::point_plane_clearance_rat(&coeffs, &o, r2) != Orient::Positive {
                    // ★ The face-level test reads each face's own vertices, which are realized
                    // world coordinates — so it needs no frame guard of its own; `coeffs` above is
                    // already the world description (a class without one never reaches here). The
                    // guard this replaces refused every moved class outright, which is what kept a
                    // translated body out of the cylinder roads.
                    let spans: Vec<[nacre_exact::Rat; 2]> = footprints
                        .get_or_insert_with(|| lateral_footprints(faces, cyl.surf))
                        .iter()
                        .map(|f| f.span.expect("a listed footprint has a span"))
                        .collect();
                    if !wall_faces_clear(model, faces, plane_ix, c, &coeffs, &o, &m, r2, &spans)? {
                        // ★ **The record-and-pass arm**: a wall whose plane runs **within**
                        // the radius — any
                        // `0 ≤ d < r`, the through-axis wall included — and whose faces did not
                        // clear the strip is recorded and passed; the tracer's ruling/chord arms
                        // fire only on recorded pairs and read the *faces* (a ruling piece comes
                        // from a face's own cycles, a cap's chord is the true meet of that class's
                        // plane with the disk), so a recorded face that misses the lateral
                        // contributes nothing. The record is «may meet; the tracer decides».
                        //
                        // ★ The offset (`0 < d < r`) keeps no refusal: with the region
                        // emitter, bosses and bores, rational and irrational
                        // rulings, axis inside or outside the plate, a split bore — every one
                        // assembles, validates clean and answers the exact volume oracle. The
                        // arithmetic (`plane_plane_cylinder`'s roots, `ruling_side`, the chart's
                        // circular order) never assumed the diameter; only the vocabulary did.
                        //
                        // ★★ **Nothing keeps a refusal here.** Lifting it assembles
                        // a volume-correct solid and `validate` cannot see the contact
                        // (`Cut` returns `Ok`, `validate` clean, no net sees
                        // the line), which is why the answer is a judge and not a fence. A
                        // tangency is a double
                        // root and three roads spell two distinct roots, so a
                        // `crossings.insert` here would put the pair on the ruling
                        // road and assemble nothing; not recording it never reaches those roads.
                        //
                        // ★★★★★ **A tangency is a graze, not a crossing** — so it passes, and it
                        // is **not** recorded. `crossings`' proposition is "the plane runs *within*
                        // the radius", which a tangent plane does not; the three roads gated on
                        // that record assert it in a `debug_assert` (`arrangement`' circular-
                        // hole arm, `chord_on_class`, `rulings_on_class`), and all three stay true
                        // because this pair never reaches them. The arrangement then sees nothing here, which is right: the line
                        // divides no cell of this plane and stations no sector of the chart.
                        //
                        // What the pair can still do is pinch the *result*, and that is a question
                        // about the operation — so the geometry is stated in a row and
                        // `assembly::tangency_reject` asks `keep` — beside `self_touch_reject`,
                        // where the grouping can also say whether the pieces share a solid.
                        if nacre_exact::point_plane_clearance_rat(&coeffs, &o, r2) == Orient::Zero {
                            tangencies.extend(tangency_rows(
                                model, faces, plane_ix, geom, n_a, c, ci, &coeffs, &cyl.def,
                                cyl.surf, cyl.owner,
                            ));
                        } else {
                            crossings.insert((c, ci));
                        }
                    } else if nacre_exact::point_plane_clearance_rat(&coeffs, &o, r2)
                        == Orient::Negative
                        && {
                            // ★ Only a footprint that can be **stated** proves a cut; an
                            // unstatable one (a rotated class, no span) proves nothing either way
                            // and keeps the old silence — recording on it drew rulings no face
                            // has (measured, the rigid-motion oracle).
                            let fps = footprints
                                .get_or_insert_with(|| lateral_footprints(faces, cyl.surf));
                            !fps.is_empty() && !oblique_plane_clears(&cyl.def, fps, &coeffs)
                        }
                    {
                        // ★ **The class cuts this lateral face though no face of the other solid touches
                        // it**: the refusal question is face against face and it is clear, but
                        // the arrangement on this class still needs the ruling — every face of *this* solid
                        // that crosses the class ends on it (a cap's section ends where its rim arc meets
                        // the plane), and without the ruling that end dangles and the walk doubles back
                        // (measured: a gusset beside a filleted plate, its side plane within the fillet's
                        // radius). So the record is «the plane cuts the face», not «the faces touch»; the
                        // tangent plane (clearance `Zero`) is the tangency's own record.
                        crossings.insert((c, ci));
                    }
                }
            }
        }
    }

    // The oblique record's one reader today (see the arm above).
    if !oblique.is_empty() {
        return Err(reject(RejectReason::ObliqueCylinderCut));
    }

    // ★★★★★ **The pair rule asks classes of different owners whether their faces share a point,
    // and writes what it cannot prove.** The proposition is "the two classes share no face". Three
    // things decide it, in order of cost. Two classes of one solid share none by construction — a
    // valid solid's faces meet only along their edges — so those pairs are not asked (the
    // fillets of one plate, the two half cylinders of one slot). The distance between the axes
    // exceeding the radius **sum** proves it for the two infinite surfaces, whatever their
    // orientation, and decides most inputs. ★ It is not the whole rule — that would make a fact
    // about surfaces read as a fact about faces: a stud fused through a cube and then a second stud
    // across it would be refused because their *axes* cross, though the first stud's remaining
    // faces
    // sit past `|z| = 0.5` and the second's whole surface within `|z| = 0.2`. So a pair the
    // distance cannot clear asks the faces themselves — the same question the plane–cylinder arm
    // asks per face (`face_clears_footprint`), spelled for a lateral face's reach against the
    // other's along every direction the pair can state ([`lateral_faces_clear`] over
    // [`separating_dirs`] — each axis, and for skew axes the common perpendicular, which is where
    // a fillet beside a crosswise drill is seen apart). Parallel axes take the same door, after
    // one more surface fact — one infinite cylinder strictly **inside** the other never meets it
    // ([`nacre_exact::cylinders_nested`]: a pin in a bore, a smaller pin stacked on a boss) —
    // and the reach along a parallel axis has no radial term, so what is left is the axis spans.
    // ★ Except **one surface under
    // two handles** — parallel axes on one line, one radius ([`same_surface`]): a `translate`d twin
    // or a restatement with another `ref_dir`. The arrangement has no name for two classes on one
    // surface (their circles coincide on every ⊥ class), so that pair is refused as the coincident
    // pair it is; interning cylinders by geometry is a later cell's.
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
    // The pair record's one reader today — the place that hands the record to the
    // cylinder–cylinder road when it exists.
    //
    // ★★★★★ **An obligation the cell that opens this refusal inherits**. The arrangement
    // leans on this line for a proposition of its own: a class may carry a circle (⊥ one
    // cylinder) and rulings (∥ another) at once, and that is safe **only while the two never
    // meet** — neither split cuts the other's kind, so a crossing they both walked past would be
    // a node no road mints, and the point is `plane ∩ cylinder ∩ cylinder`, a nested radical with
    // no name. Two such edges meeting means two lateral faces share a point, which is exactly what
    // this refusal denies. `ClassEdges::of` carries the matching `debug_assert` and the derivation;
    // when this population is admitted, that net has to become a **shipped** check there.
    if !cyl_pairs.is_empty() {
        return Err(reject(RejectReason::CylinderPairContact));
    }
    // The rulings road's record rides out beside the table (`PlaneSetup::crossings`): the
    // pairs the record-and-pass arm above admitted without a clearance proof. Empty for every
    // population outside the rulings road.
    Ok((cyls, crossings, tangencies))
}

/// **Does every face on plane class `c` provably miss this cylinder?** — the boundary question
/// [`cylinder_gate`]'s wall rule asks once the cheaper plane-level one has failed.
///
/// `Err` is the honest "could not decide exactly"; `Ok(false)` means some face was not shown to
/// clear, which is the wall refusal's whole content.
///
/// ★ **No table is built for this.** The scan runs only on the class that failed the plane test —
/// rare — so walking the face rows there costs nothing on the common path and allocates nothing.
/// A per-class face table computed for every boolean would be built
/// for everyone and read by almost no one. ★★ The `spans` this takes is held to the same rule: the
/// caller builds it on first use, not per cylinder — measured, an ordinary cut over a 16-bore
/// plate wants it **zero** times.
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
        let Some(fh) = fi.face else {
            return Err(reject(RejectReason::CylinderGateUndecided));
        };
        seen += 1;
        if !face_clears_footprint(model, model.face(fh), coeffs, o, m, r2, spans)? {
            return Ok(false);
        }
    }
    // A class exists because faces made it, so finding none is a wiring failure rather than an
    // input — and "every one of no faces clears" is a pass this must not hand out by default.
    if seen == 0 {
        return Err(reject(RejectReason::CylinderGateUndecided));
    }
    Ok(true)
}

/// **A wall plane exactly `r` from a cylinder's axis, and everything a verdict needs about it.**
///
/// ★★★★★ **A tangency is a graze, not a crossing.** The plane meets the lateral in one line and
/// *divides nothing*, so it earns no [`PlaneSetup::crossings`] record — that set's proposition is
/// "the plane runs **within** the radius", and the three roads gated on it assert exactly that
/// (☑ `crossings.contains` has three readers, each with the matching `debug_assert`). The
/// arrangement therefore never sees this pair, which is why opening the gate needs no arrangement
/// change at all (measured: of 21 cells, the 15 outside the third-plane population assemble with
/// `validate` clean and exact volumes).
///
/// What such a plane can still do is **pinch the result**: near the tangent line the material
/// splits into three regions — the lens inside the cylinder, the **two** wedges between the
/// parabola and the plane, and the far half-space — and whether the kept ones hang together is a
/// question only [`crate::draft::keep`] can answer. So this is a *record*, not a verdict:
/// the gate states the geometry, the operation decides.
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
    /// **The wall face has vertices on both sides of the tangent line** — so the contact is a
    /// *segment*, not a corner grazing it at a point. A point tangency is a
    /// valid solid, so without this the verdict would accuse a shape it cannot convict. ★ Needed
    /// because [`nacre_exact::cylinder_strip_side`] answers `StripSide::Inside` for `U = 0`,
    /// which its own doc calls *"merely conservative"* on an empty strip — and a tangent plane's
    /// strip is exactly that, so "not clear of the strip" says nothing on its own here.
    ///
    /// ★★ **The axis span is deliberately *not* part of this** — see [`face_straddles_line`]: it is
    /// already guaranteed, and asking again with the open-interval reading throws away every corner
    /// of a face flush with the cap planes (☑ measured on the frozen `bore-slab` shape).
    pub(crate) straddles: bool,
    /// **Another plane class holds this tangent line.** Such a plane is a *secant*: substituting
    /// `u = c·v` into `(u+r)² + v² = r²` gives `v·(v(1+c²) + 2cr) = 0`, so it meets the cylinder in
    /// this ruling *and* one more, runs within the radius, and is recorded as a crossing. The local
    /// picture is then **six** regions, not three, the two wedges can take different `keep`s, and
    /// the record below cannot speak.
    ///
    /// ★ A whole **disk** wall face clearing the tangent line is read by the footprint reader
    /// and writes no row at all, so what reaches here is the abstention's true subject.
    /// Abstain rather than answer. (☑ Both constructed members of
    /// this population reach `CoincidentNodes` in the arrangement anyway — two samples are not a
    /// population claim, so the abstention stands.)
    pub(crate) line_in_another_plane: bool,
    /// Nothing above could be stated exactly (a face with no rational world description, an
    /// overflow). Same answer as `line_in_another_plane`, different cause.
    pub(crate) undecided: bool,
    /// The tangency point, taken at the **middle of the lateral face's own span** rather than at
    /// the axis origin — `RejectWhere`'s doc asks for a witness that is actually there.
    pub(crate) witness: Point3,
}

/// **How a class's rational coefficients relate to a face's stored plane** — `+1` when the two
/// describe the same direction, `-1` when they oppose. The same reading as
/// `arrangement::world_rat_sense`, asked of a face rather than a class root: both must be nonzero
/// in the component compared, because `raw` is `f64` and a component it rounds to zero would hand
/// back a sign with nothing behind it.
fn rel_to_stored(coeffs: &[nacre_exact::Rat; 4], plane: &Plane) -> Option<i8> {
    let raw = plane.coefficients();
    let zero = nacre_exact::Rat::from_int(0);
    let i = (0..4).find(|&i| coeffs[i] != zero && raw[i] != 0.0)?;
    Some(if (coeffs[i] > zero) == (raw[i] > 0.0) {
        1
    } else {
        -1
    })
}

/// Sign of `n·p + d` for a rational point — the side of the plane `p` is on. `None` on overflow.
fn plane_side_of_rat(coeffs: &[nacre_exact::Rat; 4], p: &[nacre_exact::Rat; 3]) -> Option<i8> {
    let mut acc = coeffs[3];
    for k in 0..3 {
        acc = acc.checked_add(coeffs[k].checked_mul(p[k])?)?;
    }
    Some(match acc.cmp(&nacre_exact::Rat::from_int(0)) {
        core::cmp::Ordering::Greater => 1,
        core::cmp::Ordering::Less => -1,
        core::cmp::Ordering::Equal => 0,
    })
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
    let mut den = nacre_exact::Rat::from_int(0);
    for k in 0..3 {
        num = num.checked_add(n[k].checked_mul(o[k])?)?;
        den = den.checked_add(n[k].checked_mul(n[k])?)?;
    }
    // `Rat` has no division: the reciprocal is the exact inverse and `Rat::new` refuses `0`,
    // which is precisely the degenerate normal this must not divide by.
    let q = num.checked_mul(nacre_exact::Rat::new(den.denom(), den.numer())?)?;
    let mut out = [nacre_exact::Rat::from_int(0); 3];
    for k in 0..3 {
        out[k] = o[k].checked_sub(q.checked_mul(n[k])?)?;
    }
    Some(out)
}

/// **Does any *other* plane class hold this tangent line?** — the exact condition under which the
/// three-region model above stops being complete. Two rational questions per class: the line's
/// direction lies in the plane (`n·m = 0`), and its base point is on it.
fn line_lies_in_another_class(
    geom: &[WorkingPlane],
    wall: usize,
    base: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
) -> bool {
    for (k, wp) in geom.iter().enumerate() {
        if k == wall {
            continue;
        }
        let Some(w) = wp.world_rat else { continue };
        if nacre_exact::dot_sign_rat(&[w[0], w[1], w[2]], m) != nacre_exact::Orient::Zero {
            continue;
        }
        if plane_side_of_rat(&w, base) == Some(0) {
            return true;
        }
    }
    false
}

/// **The rows a tangent `(class, cylinder)` pair writes** — one per (wall face, lateral face) pair
/// that did not clear the footprint. A face that *did* clear cannot reach the tangent line, so it
/// contributes nothing; that is why collecting here neither widens nor narrows the rule.
///
/// ★★★★★ **This walk cannot fail.** [`wall_faces_clear`] short-circuits on its first non-clearing
/// face and raises on shapes it cannot read; walking further and *propagating* those raises would
/// change which reason a tangency wears. Every failure here becomes `undecided` on the row
/// instead — the verdict then abstains, which is the same answer with an honest name.
#[allow(clippy::too_many_arguments)]
pub(super) fn tangency_rows(
    model: &Model,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    geom: &[WorkingPlane],
    n_a: usize,
    c: usize,
    ci: usize,
    coeffs: &[nacre_exact::Rat; 4],
    def: &nacre_topo::CylinderDef,
    surf: Handle<Surface>,
    owner: SolidSide,
) -> Vec<Tangency> {
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let side_of = |i: usize| {
        if i < n_a { SolidSide::A } else { SolidSide::B }
    };
    // The line's base and whether a third plane holds it — one answer for the whole pair.
    let base = tangency_foot(coeffs, &o);
    let line_in_another_plane = base
        .as_ref()
        .is_some_and(|b| line_lies_in_another_class(geom, c, b, &m));
    // Which side of the wall plane the cylinder is on — the lens' side. ★ A tangency puts the
    // whole cylinder on one side, and `0` would mean the axis lies *in* the plane, which needs
    // `r = 0`; `CylinderDef` refuses that at construction (`NonPositiveRadius`, "positive by
    // construction"). Filtered anyway rather than assumed: a `Some(0)` here would compare unequal
    // to every material side and hand the verdict a silent `false`, and this file's own rule is
    // that an answer it cannot state becomes `undecided`, never a default.
    let lens_side = base
        .as_ref()
        .and_then(|_| plane_side_of_rat(coeffs, &o))
        .filter(|&s| s != 0);
    // ★ Built **once**, not once per face — the same warning `wall_faces_clear`'s caller carries
    // (a cut over a plate with 16 bores once built this table 16 times and read it 0).
    let spans = lateral_spans(faces, surf);
    // The lateral faces of this cylinder, with their own orientation and span.
    let laterals: Vec<(usize, &CylFaceInfo)> = faces
        .iter()
        .enumerate()
        .filter_map(|(i, row)| match row {
            FaceRow::Cylinder(cf) if cf.surf == surf => Some((i, cf)),
            _ => None,
        })
        .collect();
    let mut out = Vec::new();
    for (fi_ix, row) in faces.iter().enumerate() {
        let (ClassIx::Plane(k), FaceRow::Plane(fi)) = (plane_ix[fi_ix], row) else {
            continue;
        };
        if k != c {
            continue;
        }
        // ★ **Same-solid pairs write no row**. A wall tangent to its own solid's
        // cylinder is that solid's smooth edge — a fillet — not a contact between operands. The
        // judge would skip such a row anyway (`straddles` is false for a face that ends on the
        // line), but a row it could not decide would refuse the boolean by a name that is about
        // nothing. The same sentence the pair loop stands on: a valid solid's faces do not meet.
        if matches!(
            (side_of(fi_ix), owner),
            (SolidSide::A, SolidSide::A) | (SolidSide::B, SolidSide::B)
        ) {
            continue;
        }
        let (cleared, unread) = fi
            .face
            .map(|fh| {
                let got = face_clears_footprint(model, model.face(fh), coeffs, &o, &m, r2, &spans);
                #[cfg(feature = "tangency-trace")]
                if got.is_err() {
                    let f = model.face(fh);
                    let arcs = f
                        .outer
                        .half_edges
                        .iter()
                        .filter(|he| {
                            !matches!(model.edge_curve(he.edge), nacre_geom::Curve::Line(_))
                        })
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
                // ★★★★★ **A face this could not read is «unknown», not «innocent».** The doc
                // above has always said every failure here becomes `undecided` on the row — and
                // it did not: the failure was folded to «did not clear», a row was written, and
                // the row then carried a `straddles` that had **never been computed** into the
                // verdict. `!straddles` acquits, so an unreadable wall face was quietly cleared
                // of pinching the result. Measured: every abstaining row in the corpus is such a
                // face, and only `line_in_another_plane` kept the acquittal from being reached.
                (got.as_ref().is_err(), got.unwrap_or(false))
            })
            // No face row at all is the same kind of silence.
            .map(|(unread, cleared)| (cleared, unread))
            .unwrap_or((false, true));
        if cleared {
            continue; // a face that clears the footprint cannot reach the tangent line
        }
        // The wall face's material side, w.r.t. the class's coefficients. ★ The outward direction
        // is `n_out = orient_sign · plane.normal()` (`FaceInfo::n_out`, "the single source of
        // outward"), material lies opposite it, so the answer is `−orient_sign · rel` where `rel`
        // relates the *stored* plane to the class's rational one. ★★ `world_rat` is **not** that
        // relation — measured: the seed plane `x = 0` carries `world_rat = (1,0,0,0)` while its
        // stored normal is `−x`. The sign that does relate them is spelled once already, in
        // `arrangement::world_rat_sense`, and this is that spelling asked of a face.
        let mat_side = rel_to_stored(coeffs, &fi.plane).map(|rel| -fi.orient_sign * rel);
        // A property of the wall face and the line, not of which lateral is paired with it.
        let straddles = fi
            .face
            .map(|fh| face_straddles_line(model, model.face(fh), coeffs, &o, &m, r2))
            .unwrap_or(false);
        for (cy_ix, cf) in &laterals {
            let witness = base.as_ref().zip(cf.footprint.span).and_then(|(b, span)| {
                let mid = span[0]
                    .checked_add(span[1])?
                    .checked_mul(nacre_exact::Rat::new(1, 2)?)?;
                let mut p = [0.0f64; 3];
                for j in 0..3 {
                    p[j] = b[j].checked_add(mid.checked_mul(m[j])?)?.to_f64();
                }
                Some(Point3::from_array(p))
            });
            let undecided = unread
                || base.is_none()
                || lens_side.is_none()
                || mat_side.is_none()
                || witness.is_none()
                || cf.def.is_none();
            #[cfg(feature = "tangency-trace")]
            #[allow(clippy::print_stderr)]
            {
                eprintln!(
                    "TROW liap={line_in_another_plane} straddles={straddles} undecided={undecided}"
                );
            }
            out.push(Tangency {
                wall: c,
                cyl: ci,
                wall_solid: side_of(fi_ix),
                cyl_solid: side_of(*cy_ix),
                lens_in_wall_solid: lens_side.is_some() && lens_side == mat_side,
                cyl_orient: cf.orient_sign,
                straddles,
                line_in_another_plane,
                undecided,
                witness: witness.unwrap_or(Point3::from_array([f64::NAN; 3])),
            });
        }
    }
    // ★ A tangency of a solid with **itself** is not this operation's business: for one solid to
    // own both the wall face and a lateral tangent to it along that line, its own material would
    // have to be the two wedges alone or the lens and the far side alone — and both of those *are*
    // the pinch. Such an operand is invalid before any boolean runs, so dropping the row discards
    // no information about a valid input.
    out.retain(|t| t.wall_solid != t.cyl_solid);
    out
}

/// **Does the face reach across the tangent line?** — the proof that the contact is a *segment*
/// rather than a corner grazing the line at one point (a valid tangency).
/// `StripSide::Plus`/`Minus` are the two sides of the zero-width strip a tangent plane cuts, and
/// a face with a corner on each has the line running through its hull.
///
/// ★ **The axis overlap needs no test here** — the caller only asks this of a face that already
/// failed [`face_clears_footprint`], and that function returns `false` exactly when the face
/// clears *neither* axis: not across the strip **and**, for every lateral span, not wholly at or
/// below its start nor wholly at or above its end. The second half is the overlap, already proved.
/// Re-testing it with [`Corner::reaches`] would be worse than redundant: the span is read **open**
/// (`planes`'s own note), so a face whose corners sit exactly on the cap planes — a slab cut
/// flush with a bore's own height, the frozen `bore-slab` shape — would have every corner dropped
/// and the straddle read as absent. ☑ Measured: that is exactly what happened.
fn face_straddles_line(
    model: &Model,
    face: &Face,
    coeffs: &[nacre_exact::Rat; 4],
    o: &[nacre_exact::Rat; 3],
    m: &[nacre_exact::Rat; 3],
    r2: &nacre_exact::BigRat,
) -> bool {
    use nacre_exact::StripSide;
    let (mut plus, mut minus) = (false, false);
    for he in &face.outer.half_edges {
        // ★ **The same reader the clearance road uses.** This used to take rational vertices and
        // skip everything else, so a face ringed by tangent corners straddled nothing as far as
        // this could tell — and «does not straddle» acquits. A piece it still cannot read is not
        // a silent acquittal either: such a face fails the clearance call as well, and its row is
        // marked undecided there.
        let Ok(corner) = corner_of(model, face, he) else {
            continue;
        };
        if !corner.on_plane(coeffs) {
            continue;
        }
        match corner.strip_side(coeffs, o, m, r2) {
            StripSide::Plus => plus = true,
            StripSide::Minus => minus = true,
            // A corner sitting *on* the line spans nothing, and neither does a piece that reaches
            // the strip without crossing out the far side.
            StripSide::Inside => {}
            // ★ **A piece with width can be the whole straddle by itself.** Two corners on
            // opposite sides is one way for the line to run through the face; a single piece whose
            // interior lies on both sides is the other, and it is the same fact.
            StripSide::Crosses => {
                plus = true;
                minus = true;
            }
        }
    }
    plus && minus
}

/// The axis-parameter spans this cylinder's lateral faces occupy — one per face, in the scale
/// `axis_param_of_plane` produces.
///
/// ★★ **Empty means "unusable", never "nothing in the way".** A span this returns is a rectangle
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

/// **The footprints of a cylinder class's lateral faces** — [`lateral_spans`]' rule with the
/// angular extent alongside: empty if any face's span could not be stated, else one
/// footprint per face.
fn lateral_footprints(faces: &[FaceRow], surf: Handle<Surface>) -> Vec<Footprint> {
    let mut out: Vec<Footprint> = Vec::new();
    for row in faces {
        let FaceRow::Cylinder(cf) = row else { continue };
        // ★ Matched by **handle**, not by geometry: two operands may state the same cylinder
        // twice, and each statement gets its own rows, its own bands clipped to its own spans,
        // and its own turn through the gate's cylinder loop.
        if cf.surf != surf {
            continue;
        }
        let Some(span) = cf.footprint.span else {
            return Vec::new();
        };
        // ★ The reader below asks "every vertex at or below `span[0]`, or every one at or above
        // `span[1]`", which is only the interval's outside while the pair is ordered — reversed,
        // that same phrasing reads "outside the *union* of two half-lines" and would clear a face
        // sitting squarely in the band. `lateral_t_range` orders it; this is where that is
        // relied on, so this is where it is said.
        debug_assert!(span[0] <= span[1], "a span is stated low end first");
        out.push(cf.footprint);
    }
    out
}

pub(super) fn lateral_reach(
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
    let (lo_off, rho2_lo, hi_off, rho2_hi) = arc_extent(fp.theta.as_ref(), r2, &m, d)?;
    Some(Reach {
        lo: lo.checked_add(lo_off)?,
        hi: hi.checked_add(hi_off)?,
        rho2_lo,
        rho2_hi,
    })
}

/// **Do two reaches provably miss each other?** — one lies wholly beyond the other, at either end.
///
/// Each end is a rational base and a radical, so a gap has to clear **two** roots:
/// `b.lo − √b.rho2_lo > a.hi + √a.rho2_hi` is `g > √p + √q`, which
/// [`nacre_exact::exceeds_root_sum`] answers exactly. `None` is `Rat` overflow forming a gap.
///
/// ★ **Closed against closed**: equality is one face's rim touching the other at a point, which is
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
pub(super) fn reach_clears(
    reach: &Reach,
    lo: nacre_exact::Rat,
    hi: nacre_exact::Rat,
) -> Option<bool> {
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

/// **Does every lateral face of this cylinder provably miss the plane `n·p + d = 0`?** — the oblique
/// arm's face question. The plane's station in `n·p` units is `−d`; a face clears when
/// its reach along `n` ([`lateral_reach`]) is disjoint from that one point. An empty span list is
/// "unusable" ([`lateral_spans`]'s doc), never clear; so is overflow.
fn oblique_plane_clears(
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
/// lateral face of `b`, each pair asked along [`separating_dirs`]. A pair clears when **some**
/// direction separates them; the classes clear when **every** pair does. Then no point of `b`'s
/// faces lies on `a`'s, so the two classes share no face — the proposition the arrangement needs,
/// which the surface distance ([`nacre_exact::cylinders_clear`]) is only one sufficient
/// condition for.
///
/// ★ The two axes used to be spelled here as two questions of different shapes — one class's
/// **span** against the other's **reach**, asked both ways round. They are the one question
/// [`separated`] asks at `d = m_a` and `d = m_b`: along its own axis a face's reach *is* that
/// span (`d ∥ m` leaves no radial term), so the interval the old spelling built by hand is what
/// [`lateral_reach`] returns. One rule, a list of directions, and the asymmetry is gone.
///
/// Empty [`lateral_spans`] is "unusable" (its doc), never clear: it becomes one footprint of
/// unknown span, which the reach reads as unbounded along anything not perpendicular to the axis.
/// `None` from any direction is not clear.
fn lateral_faces_clear(faces: &[FaceRow], a: &WorkingCyl, b: &WorkingCyl) -> bool {
    // A class whose spans cannot be stated still has faces to ask about, with an unbounded
    // reach along anything not perpendicular to its axis: one footprint that says so.
    let listed = |surf: Handle<Surface>| -> Vec<Footprint> {
        let v = lateral_footprints(faces, surf);
        if v.is_empty() {
            vec![Footprint {
                span: None,
                theta: None,
            }]
        } else {
            v
        }
    };
    let (fa, fb) = (listed(a.surf), listed(b.surf));
    let parallel = nacre_exact::parallel_rat(&a.def.dir(), &b.def.dir());
    let dirs = separating_dirs(&a.def, &b.def);
    fa.iter().all(|x| {
        fb.iter().all(|y| {
            let apart = dirs
                .iter()
                .any(|d| separated(&a.def, x, &b.def, y, d) == Some(true));
            apart || (parallel && cross_sections_clear(a, x, b, y) == Some(true))
        })
    })
}

/// **The separating directions a cylinder pair can state rationally** — and why these are all of
/// them.
///
/// What a footprint describes is a box in `(axis parameter × angle)`, and the boundary of such a
/// face is made of **cap planes** (normal = the axis), **rulings** (direction = the axis) and
/// **arcs**. A separating axis is a face normal or an edge×edge, so the rational candidates are
/// each axis and — for skew axes — the two rulings' cross product, the common perpendicular. What
/// is left out is the arcs' radial continuum, which no rational direction names; the tool for that
/// is [`cross_sections_clear`], and it needs the common cross-section chart that only **parallel**
/// axes have. So a pair the three cannot separate is refused as "not shown to clear", which is
/// true, and completing the test further means the same test on a richer candidate set.
///
/// ★★ **The third direction is the face-level twin of a rung above it.** With a whole circle and
/// no span, `d ⊥ m_a` and `d ⊥ m_b` make both reaches `d·o ± r|d|`, so "apart" reads
/// `|d·(o_b − o_a)| > (r_a + r_b)|d|` — [`nacre_exact::cylinders_clear`]'s skew branch, letter
/// for letter. The faces answer the same question the infinite surfaces do, with their own extent.
///
/// Parallel axes have no third direction (the cross product is zero, and the surface rung is the
/// distance rule with its own parallel branch); overflow forming it simply leaves the list short.
pub(super) fn separating_dirs(
    a: &nacre_topo::CylinderDef,
    b: &nacre_topo::CylinderDef,
) -> Vec<[nacre_exact::Rat; 3]> {
    let (ma, mb) = (a.dir(), b.dir());
    let mut out = vec![ma, mb];
    if !nacre_exact::parallel_rat(&ma, &mb) {
        if let Some(perp) = nacre_exact::cross3_rat(&ma, &mb) {
            out.push(perp);
        }
    }
    out
}

/// **Do face `x` of `a` and face `y` of `b` provably miss each other along `d`?** — each face's
/// reach along the one direction ([`lateral_reach`]), asked whether they are disjoint. Sound for
/// any `d`: a reach is a bounding interval of the face's projection, and two sets whose
/// projections miss cannot share a point. `None` is "not proved" — an unstatable span, or
/// overflow.
pub(super) fn separated(
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
/// the axis, `û₂ = m̂ × û₁` — rational exactly when both norms are (`inv_sqrt_exact`; every prism
/// on a world frame), else `None`. Each face's arc is its angular extent, or the whole circle
/// when none is stated; the question is [`nacre_geom::mixed::arcs_share_a_point`]'s, and a
/// touch is not clear. Sound for any face of the class because the extent is the footprint's
/// bounding arc: wider than a notched face, never narrower.
fn cross_sections_clear(
    a: &WorkingCyl,
    x: &Footprint,
    b: &WorkingCyl,
    y: &Footprint,
) -> Option<bool> {
    use nacre_exact::Rat;
    use nacre_geom::mixed::ArcSpec;
    let zero = Rat::from_int(0);
    let (o, m, e) = (a.def.origin(), a.def.dir(), a.def.ref_dir());
    let (mm, em) = (
        nacre_exact::dot3_rat(&m, &m)?,
        nacre_exact::dot3_rat(&e, &m)?,
    );
    let mut e1 = [zero; 3];
    for k in 0..3 {
        e1[k] = mm.checked_mul(e[k])?.checked_sub(em.checked_mul(m[k])?)?;
    }
    let inv_e1 = nacre_exact::inv_sqrt_exact(nacre_exact::dot3_rat(&e1, &e1)?)?;
    let inv_m = nacre_exact::inv_sqrt_exact(mm)?;
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
    let spec = |c: &WorkingCyl, fp: &Footprint| -> Option<ArcSpec> {
        let centre = chart(&c.def.origin())?;
        // The 2-D arc court works in `Rat`; a wide square is not judged here (`None`).
        let r2 = c.def.r2().narrow()?;
        let (start, end) = match &fp.theta {
            Some(arc) => (
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
            None => {
                // A whole circle's seam, `centre + (r, 0)`, is a rational point only for a
                // rational radius; without one this pair is not judged here (`None`, the
                // caller's honest road).
                let p = add2(centre, [c.def.radius_exact()?, zero])?;
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
