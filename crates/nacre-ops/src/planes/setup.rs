use super::*;
/// The minimal per-op plane table two solids share: the
/// concatenated plane list (`a`'s then `b`'s), the face→index map, and each solid's
/// [`combinatorics::EdgeFaces`]. Built once and shared: indices into the returned `planes`/`surf_ix`
/// are common to both solids, so a vertex of `a` and a face of `b` compose in one index space.
/// Destructure it with `..` (`let PlaneSetup { planes: faces_tab, geom: planes, plane_ix, .. } = …`):
/// the tables here grow as the arrangement learns to say "plane" and "face" in different index
/// spaces, and a positional tuple made every one of those steps touch all ~25 call sites.
///
/// The plane classes (`canon`) are computed here to build `geom`/`plane_ix` and then dropped — the
/// dense `plane_ix` is the only face→plane map anything downstream needs, so the sparse union-find
/// output does not escape.
pub(crate) struct PlaneSetup {
    pub(crate) planes: Vec<FaceRow>,
    pub(crate) surf_ix: HashMap<Handle<Face>, usize>,
    pub(crate) inc_a: combinatorics::EdgeFaces,
    pub(crate) inc_b: combinatorics::EdgeFaces,
    /// Where `a`'s faces end and `b`'s begin in `planes`. The concatenation always created this
    /// boundary; it was just never written down, so every later "whose face is this?" had to
    /// rebuild it.
    pub(crate) n_a: usize,
    /// The arrangement's planes, densely indexed — see [`dense_planes`].
    pub(crate) geom: Vec<WorkingPlane>,
    /// `plane_ix[face]` is that face's class — a plane index into `geom`, or a cylinder class
    /// ([`ClassIx`]).
    pub(crate) plane_ix: Vec<ClassIx>,
    /// Whose faces each plane class carries — see [`class_owners`].
    pub(crate) class_owner: Vec<Option<SolidSide>>,
    /// The cylinder classes, in [`ClassIx::Cyl`] numbering order — empty for an all-planar
    /// boolean. Filled by the population gate, which is also what refuses the interactions this
    /// milestone does not build.
    pub(crate) cyls: Vec<WorkingCyl>,
    /// The gate's carried answer for the rulings road: `(plane class, cylinder class)`
    /// pairs allowed through **without** a clearance proof — see
    /// [`combinatorics::TraceInput::crossings`]. ★ It used to say "always empty while the wall rule
    /// refuses that population" — the gate-opening cell arrived, and the `crossings.insert` below
    /// fills it for a wall whose plane holds the axis exactly.
    pub(crate) crossings: std::collections::HashSet<(usize, usize)>,
    /// **The tangent `(wall face, lateral face)` pairs** — the graze twin of `crossings`, and the
    /// reason they are two sets rather than one: a recorded *crossing* says "this plane runs within
    /// the radius, so it cuts the lateral in two rulings", and a tangency says the opposite ("it
    /// touches along one line and cuts nothing"). Merging them would make the record's own
    /// proposition false and break the four `debug_assert`s that lean on it. See [`Tangency`].
    pub(crate) tangencies: Vec<Tangency>,
    /// How this operation judges, and where its evidence goes — the two facts that belong to the
    /// operation rather than to any one plane. The caller pairs them with a table to make a
    /// [`Judge`].
    pub(crate) standard: Standard,
    pub(crate) notes: Notes,
}

/// Which sub-phase of [`plane_index_setup`] a [`Watch`] charges — spike instrumentation, and only
/// in a test build (see `arrangement::phase`).
pub(crate) enum Sub {
    TriPt3,
    Std,
    Collect,
    Edges,
    Classes,
    Dense,
}

/// Times the scope it is charged from, or does nothing at all in a release build.
pub(crate) struct Watch(#[cfg(test)] std::time::Instant);

impl Watch {
    pub(crate) fn new() -> Self {
        Watch(
            #[cfg(test)]
            std::time::Instant::now(),
        )
    }
    #[allow(unused_variables)]
    pub(crate) fn charge(self, which: Sub) {
        #[cfg(test)]
        {
            use crate::arrangement::phase;
            let c = match which {
                Sub::TriPt3 => &phase::S_TRIPT3,
                Sub::Std => &phase::S_STD,
                Sub::Collect => &phase::S_COLLECT,
                Sub::Edges => &phase::S_EDGES,
                Sub::Classes => &phase::S_CLASSES,
                Sub::Dense => &phase::S_DENSE,
            };
            c.fetch_add(
                self.0.elapsed().as_nanos() as u64,
                std::sync::atomic::Ordering::Relaxed,
            );
        }
    }
}

pub(crate) fn plane_index_setup(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<PlaneSetup, BoolError> {
    let (mut setup, cyl_surfs) = plane_index_setup_inner(model, a, b)?;
    // ★ The cylinder door: the population gate decides **by name** what stands in the way (an
    // oblique cut, a wall touching the lateral surface, an undecidable pair), and what passes now
    // goes on to be arranged.
    if !cyl_surfs.is_empty() {
        let (cyls, crossings, tangencies) = cylinder_gate(
            model,
            &cyl_surfs,
            &setup.geom,
            &setup.planes,
            &setup.plane_ix,
            setup.n_a,
        )?;
        setup.cyls = cyls;
        setup.crossings = crossings;
        setup.tangencies = tangencies;
    }
    Ok(setup)
}

pub(crate) fn plane_index_setup_inner(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
) -> Result<(PlaneSetup, Vec<Handle<Surface>>), BoolError> {
    let t = Watch::new();
    let mut planes = collect_planes(model, a)?;
    let n_a = planes.len();
    planes.extend(collect_planes(model, b)?);
    t.charge(Sub::Collect);
    let t = Watch::new();
    let standard = standard_for(&planes);
    t.charge(Sub::Std);
    let notes = Notes::new();
    if standard.prec > JUDGE_PREC_CAP {
        return Err(reject(RejectReason::PrecisionBudget {
            needed: standard.prec,
            cap: JUDGE_PREC_CAP,
        }));
    }
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in planes.iter().enumerate() {
        // Synthetic faces are appended later, after this table is built; every entry here is real.
        surf_ix.insert(pi.face().expect("collect_planes yields real faces"), i);
    }
    let t = Watch::new();
    let inc_a = combinatorics::edge_faces(model, a, &surf_ix)?;
    let inc_b = combinatorics::edge_faces(model, b, &surf_ix)?;
    t.charge(Sub::Edges);
    // One judging context for the whole operation: the witnesses, the standard they are held to,
    // and where the evidence goes. The face table judges first (it is what *defines* the plane
    // classes), then the dense plane table inherits the same three.
    let t = Watch::new();
    let canon = plane_classes(&Judge::new(&planes, standard, &notes));
    t.charge(Sub::Classes);
    let t = Watch::new();
    let (geom, plane_ix, cyl_surfs) = dense_planes(&planes, &canon);
    let class_owner = class_owners(&plane_ix, n_a, geom.len());
    t.charge(Sub::Dense);
    Ok((
        PlaneSetup {
            planes,
            surf_ix,
            inc_a,
            inc_b,
            n_a,
            geom,
            plane_ix,
            class_owner,
            cyls: Vec::new(),
            crossings: std::collections::HashSet::new(),
            tangencies: Vec::new(),
            standard,
            notes,
        },
        cyl_surfs,
    ))
}

/// **How precisely this operation's rotated definitions must be realized.**
///
/// The judges' error radius is `C · 2⁻ᵖʳᵉᶜ`, and `C` belongs to the model — it grows about one
/// bit per turn of rotation history and with the coordinate magnitudes. A fixed precision
/// therefore decides, silently, how long a model's history may be: at 256 bits a solid turned 245
/// times stops building, with a reject that names a symptom rather than the cause. So the
/// precision is read off the model instead.
///
/// The target is the **coincidence precision**: two things closer than this are treated as
/// coincident, and the kernel will only say so once it has *proved* the separation is below it.
/// Its default is derived rather than chosen —
///
/// - `output_precision = scale · 2⁻⁵²`, the finest distinction the `f64` coordinates this kernel
///   emits can carry. Below it nothing survives export, so distinguishing is meaningless.
/// - `coincidence_precision = output_precision · 2⁻¹²⁸`, two whole words further down. Erring low
///   only costs bits, while erring high merges features that were genuinely apart, so the
///   asymmetry says push it down; and a word is the natural unit because astro-float allocates
///   whole words anyway.
///
/// `scale` is the largest coordinate magnitude in either operand, taken over the whole table so
/// the result does not depend on traversal order (replay must reproduce it exactly).
///
/// The precision that reaches the target is then [`nacre_judge::judge_precision`]'s to compute;
/// [`JUDGE_PREC_CAP`] is where the kernel stops and says so instead, and [`CLIMB_HEADROOM`] is
/// what a single hard judgement may spend on top of it.
///
/// **Only the coincidence limit is a candidate for a setting.** Everything else here — the
/// precision, the cap, the headroom — is derived from it and from the model, because a bit count
/// means a different physical thing in every model ("256 bits" is `1e-76` for a solid turned once
/// and `1e+15` for one turned three hundred times).
/// The plane data of a row, `None` for a cylinder — the kind filter the plane-only sweeps
/// share.
#[inline]
pub(crate) fn plane_of(r: &FaceRow) -> Option<&FaceInfo> {
    match r {
        FaceRow::Plane(p) => Some(p),
        FaceRow::Cylinder(_) => None,
    }
}

// A cylinder row contributes no witness points to the standard — structurally right, not an
// omission: an axis-aligned-grade rational cylinder is the `rotated == false` case (its exact
// def needs no high-precision realization); the rotated-cylinder story is CIP's.
fn standard_for(rows: &[FaceRow]) -> Standard {
    // ★ **A face that was never moved contributes exactly nothing, so it is not asked.**
    //
    // Its `tri_pt3` are `WitnessPoint::at_nearest` of the plane's own rational points and the chain
    // is empty: an f64-representable point (base = `mantissa · 2^exp`, a power-of-two denominator
    // and a numerator within `Rat`'s 127 bits) realizes as an *exact* interval, and a decimal one
    // carries the ½-ulp bound it was stated with — nothing a replay could add. The realization has
    // no rotation error to report (`an_exact_point_demands_no_precision` in `nacre-judge`, and the
    // const assert at `TRIAL_PREC` that keeps it true).
    //
    // So the loop below used to spend a full high-precision replay per point to compute a zero —
    // measured, an axis-aligned 60-fin fold did that 24,120 times for 15.7ms and a `worst` of
    // exactly `Mag::ZERO`. The same shape was removed one level down when a stated zero replaced
    // `WitnessPoint::at`'s measurement for these points ("nine BigFloat operations to compute a
    // zero").
    //
    // `max` over the empty set is `Mag::ZERO`, which is the right answer for a model with no
    // rotation history — `precision_for` reads that as "nothing to size" and returns `TRIAL_PREC`.
    let worst = worst_trial(
        rows.iter()
            .filter_map(plane_of)
            .filter(|p| p.rotated)
            .flat_map(|p| p.tri_pt3.iter()),
    );
    // `scale`, by contrast, is every point's business: it is the model's size, and an unmoved face
    // is as far from the origin as any other.
    standard_from(
        rows.iter()
            .filter_map(plane_of)
            .flat_map(|p| p.tri_pt3.iter()),
        worst,
    )
}

/// **How deep a model may be before the operation is rejected instead.**
///
/// Not a resolution limit — the arithmetic is correct at any depth — but a **cost** limit, so it
/// is set from measured cost. A judgement's realization is quadratic-ish in the precision, and the
/// cap is placed where a single boolean's judging stays in the seconds rather than the minutes:
/// 4096 bits covers a rotation history of roughly four thousand turns (measured: `C` grows one bit
/// per turn), which is far past any real model, and a model that does exceed it is told *why*
/// rather than handed a wrong answer or an unbounded wait.
pub(crate) const JUDGE_PREC_CAP: usize = 4096;

/// **How thin a witness the kernel will still judge**, expressed as the bits a single judgement
/// may ask for *beyond* what the model itself needed.
///
/// This is a **separate budget from [`JUDGE_PREC_CAP`], and it has to be.** Sharing one absolute
/// ceiling would mean a deeply-turned model — already near the cap — leaves a hard judgement no
/// room at all, so the same sliver would be judged in a fresh model and abandoned in a turned one.
/// The model's depth and a judgement's difficulty are different quantities; only the second
/// belongs here.
///
/// It has a physical reading. A judgement's uncertainty is `(C / |cofactor|) · 2⁻ᵖʳᵉᶜ`, and the
/// model already chose `prec` so that `C · 2⁻ᵖʳᵉᶜ` clears the coincidence limit; what is left is
/// `log₂(1 / cofactor)` — the **thinness of the witness**, a needle triangle or three planes that
/// almost share a line. Two words says: a witness up to `2¹²⁸` (≈ 3·10³⁸) times more degenerate
/// than the model's own size is still judged to the end.
///
/// Two words, and not a measured number, for the same reason the coincidence limit is two words
/// below the output resolution: the error is asymmetric. Too small abandons a judgement that had
/// an answer; too large only spends bits. And measurement says there is nothing to tune — across
/// the rotation corpus and models turned 100 and 800 times, **no judgement asked for even one bit
/// beyond the model's own precision** (measured with the headroom forced to zero).
pub(crate) const CLIMB_HEADROOM: usize = 128;

/// A judging context over a hand-built table, for fixtures.
///
/// The standard is the derived default for a unit-scale model, and the collector is leaked so a
/// fixture is a one-liner — a handful of `Vec`s per test run, and nothing reads them. A fixture
/// that *does* want the evidence builds its own [`Notes`] and calls [`Judge::new`].
#[cfg(test)]
pub(crate) fn test_judge<W>(planes: &[W]) -> Judge<'_, W> {
    let notes: &'static Notes = Box::leak(Box::new(Notes::new()));
    Judge::new(
        planes,
        Standard {
            prec: 256,
            coincidence: Mag::pow2(-180),
            scale: Mag::of(1.0),
            cap: 256 + CLIMB_HEADROOM,
        },
        notes,
    )
}

/// [`standard_for`] over a bare set of definitions, with no plane table — **a test helper.**
///
/// It once served witness selection in `rotated_vertex`; that consumer is gone, and splitting
/// [`worst_trial`] out of [`standard_from`] is what surfaced it. Kept because a fixture that asks
/// "what precision does *this* point demand" wants exactly the two halves in order, and spelling
/// them out at every call site says less than the name does.
#[cfg(test)]
pub(crate) fn standard_for_points<'a>(
    pts: impl IntoIterator<Item = &'a WitnessPoint> + Clone,
) -> Standard {
    standard_from(pts.clone(), worst_trial(pts))
}

/// **The realization depth this set of definitions demands** — `max` over their trial bounds.
///
/// Split from [`standard_from`] because the two halves have nothing in common but the answer: this
/// one is **all of the cost** (a full high-precision replay per point), and the other is f64
/// arithmetic on already-known numbers. Keeping them apart is what lets a caller that already knows
/// this maximum skip straight to the second half.
///
/// **This is where a boolean spends most of what is left after the arrangement went parallel**
/// (measured: 76% of setup, and setup is 43% of the largest booleans once the trace is off the
/// critical path). Each point's trial realization is independent and they combine by **maximum**,
/// which is associative and exact — so evaluating them across cores cannot move the answer the way
/// a reassociated sum would.
fn worst_trial<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>) -> Mag {
    let pts: Vec<&WitnessPoint> = pts.into_iter().collect();
    let bounds = crate::par::map_range(pts.len(), |i| nacre_judge::trial_bound(pts[i]));
    bounds
        .into_iter()
        .fold(Mag::ZERO, |w, b| if w.lt(b) { b } else { w })
}

/// The standard for points whose worst trial bound is already known: `scale` off the f64
/// coordinates, the coincidence limit derived from it, and the precision that reaches it.
fn standard_from<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>, worst: Mag) -> Standard {
    let mut scale = 1.0f64;
    for p in pts {
        for c in p.coord() {
            scale = scale.max(c.abs());
        }
    }
    let scale = Mag::of(scale);
    let output_precision = scale.times(Mag::pow2(-52));
    let coincidence = output_precision.times(Mag::pow2(-128));
    let prec = nacre_judge::precision_for(worst, coincidence);
    Standard {
        prec,
        coincidence,
        scale,
        // Relative to this model's own depth — see [`CLIMB_HEADROOM`]. The absolute ceiling that
        // leaves is `JUDGE_PREC_CAP + CLIMB_HEADROOM`, since a model deeper than the first is
        // rejected before any judging starts.
        cap: prec + CLIMB_HEADROOM,
    }
}
