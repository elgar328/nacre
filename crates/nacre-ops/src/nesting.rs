//! **Is this cell of a plane class inside that one?** — the one engine, and the one place a
//! witness for it is named.
//!
//! **Why a module of its own.** The question is answered from two roads — the arrangement's nesting
//! (`arrangement::nest_cells`, `innermost_host`) and the coplanar merge's ownership
//! (`assembly::unify_coplanar_faces`) — and when each road had its own copy, the copies drifted:
//! five spellings of «what point may I use as a witness», four of «what does that point say about
//! the target», three of «try them until one decides, and name the failure honestly», and two of
//! the four-arm dispatch itself; the shortest — the arm that asks a **ring** against a **disk** —
//! was one witness kind wide. A plate with every corner filleted has a cap ring of eight tangencies
//! and no three-plane corner at all, so such an arm would find nothing and refuse, and the plate
//! could not enter any boolean.
//!
//! Nothing here is new geometry. There is **one** of each predicate, behind a module wall, so
//! another copy cannot be written outside it.

use super::*;
use crate::combinatorics::{self, NodeId};
use crate::planes::*;
use crate::tolerant::Judge;

/// **The circle a ring *is***, when every one of its edges is an arc of one cylinder and they
/// chain the whole way round — `None` otherwise.
///
/// **Its clauses beyond "the first edge is an arc" are guards, and none of them fired** over the
/// suite (88 acceptances, 0 rejections past that first test).
///
/// **Those 88 are not 88 questions.** Read against another measurement — *zero whole-circle
/// rings as an engine question's source over 27,000* — the two look like a contradiction. They
/// are **different populations**: this
/// counts every ring the supply asks about, and 85% of the probe's rows are not engine questions
/// at all. Both are true, and together they say the thing that matters: a ring that
/// *is* a whole circle is **never** the source of an engine question in this corpus. ⇒ giving
/// [`Cell::Ring`]'s whole-circle arm a rim
/// would be building for an empty population. **It is not built, and this is
/// why.** They are kept because each states a
/// proposition proved somewhere else — two circles on one class cannot meet (the gate), a ring is
/// a chain (measured over 153,798 rings), a mixed sense would retrace one arc — and a guard that
/// stops holding is how a producer change is meant to surface here rather than two layers down.
///
/// **A cut circle's cell is still a disk, and must be asked a disk's question.** The
/// dispatch calls a cell a `Disk` only where its half-edge is literally a `Circle`; a circle a wall
/// has **split into arcs** arrives as a ring instead, whose corners are pierce points — so the
/// three-plane names come out
/// **empty** and the road refuses with a name about rays it never cast. The shape is the same
/// circle either way, and this is what says so.
///
/// **Why the chain is the whole circle.** The gate proves every pair of cylinders clear — their
/// surfaces by more than the radius sum, parallel surfaces one inside the other, or their faces
/// ([`crate::planes::pieces_clear`]; the gate is [`crate::planes::cylinder_gate`], else
/// `CylinderPairContact`) — so two
/// circles on one class **cannot meet** — arcs that chain head-to-tail therefore all ride the same
/// circle, and running one way round (`ccw` all equal — a mixed pair would retrace one arc) closes
/// it exactly once. Measured: **every** ring with an empty probe list is either such a chain — and
/// each of those the circle's own centre answers — or a wall panel with no circle at all, which
/// refuses by its own name ([`crate::RejectReason::RingHasNoWitness`]). Against refusing every
/// empty probe list as `NoClearRay`, in cells: of the crossing census's 44, **8** build and **20**
/// are renamed, the rest being the other road's (`boolean`).
/// Counted as *cells*, not as raises — a traced boolean runs twice in `debug`, and
/// `reject_census`'s own note forbids reading raise counts as populations.
fn ring_own_circle<'a>(ring: &'a [combinatorics::RingEdge]) -> Option<&'a nacre_topo::CylinderDef> {
    let n = ring.len();
    if n < 2 {
        return None;
    }
    let arc = |e: &'a combinatorics::RingEdge| match &e.carrier {
        combinatorics::Carrier::Arc(ac) => Some(&**ac),
        _ => None,
    };
    let first = arc(&ring[0])?;
    for (i, e) in ring.iter().enumerate() {
        let a = arc(e)?;
        if a.cyl != first.cyl || a.ccw != first.ccw || e.to != ring[(i + 1) % n].node {
            return None;
        }
    }
    Some(&first.def)
}

/// **Is a rational point inside a ring?** — `Ok(None)` when *this point* cannot answer, `Err` when
/// a **value** could not be formed exactly.
///
/// **The two are different facts and the split is the point of this function.** With one
/// witness (a disk against a ring) a point *on* the ring and a class with no rational description
/// can share a refusal; a caller with **several** witnesses must not read them the same way: a
/// point on the ring has the **next witness as its remedy**, a value that cannot be formed does
/// not — the lesson `point_in_component`'s doc states for the road one dimension up.
///
/// **And it dispatches to the strongest road it can.** A ring the rational chart can name takes
/// [`nacre_geom::intersect::point_in_ring_2d_rat`], which is **half-open in y** and so decides
/// even where a ring corner sits on the ray; only a ring with pierce corners or arc steps falls
/// to `point_in_mixed_ring`, which abstains there. Anything that hands this a witness gets the
/// better answer for free.
fn rational_point_in_ring(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    p: &[nacre_exact::Rat; 3],
    ring: &[combinatorics::RingEdge],
) -> Result<Option<bool>, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let coeffs = jd.world_rat(wc).ok_or_else(undecided)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    let center = *p;
    // The class's rational chart — the one copy of that rule ([`combinatorics::Chart2dRat`]);
    // parity is affine-invariant, so the basis need not be orthonormal.
    // **A ring the chart road cannot name takes the mixed road**: pierce corners have no
    // rational coordinates and arc steps no straight chart
    // image, so the parity walks the ring step by step in ℚ(√c) instead. Rings the chart road
    // can name still take it — the mixed arm runs on exactly the population the chart road
    // cannot name.
    // A three-plane corner with no rational meet (a class with no world statement among its
    // three) is a value this road cannot form, not this witness's circumstance — the mixed walk
    // below would abstain on it, and an abstention sends the caller to its next witness, which
    // meets the same corner.
    if ring.iter().any(|e| {
        combinatorics::three_plane_name(e.node).is_some()
            && combinatorics::node_coords_rat(jd, e.node).is_none()
    }) {
        return Err(undecided());
    }
    if combinatorics::ring_is_mixed(ring) {
        return Ok(combinatorics::point_in_mixed_ring(
            jd, cyls, &coeffs, &center, ring,
        ));
    }
    let chart = combinatorics::Chart2dRat::of_normal(&n).ok_or_else(undecided)?;
    let p2 = chart.project(&center).ok_or_else(undecided)?;
    let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
    let ring2 = chart.ring(jd, &nodes).ok_or_else(undecided)?;
    Ok(
        match nacre_geom::intersect::point_in_ring_2d_rat(p2, &ring2) {
            nacre_geom::intersect::RingSide::Inside => Some(true),
            nacre_geom::intersect::RingSide::Outside => Some(false),
            // On the boundary this witness says nothing — the caller's next one may.
            nacre_geom::intersect::RingSide::OnBoundary => None,
        },
    )
}

/// **Whether the disk cylinder `a` cuts on class `wc` lies inside the disk `b` cuts there** — one
/// rational inequality on the plane: `r_a < r_b` and `(r_b − r_a)² > |c_b − c_a|²`. A disk lies
/// inside another iff its rim does, and the rims of two classes never meet (the gate proves the
/// faces clear, or the solid is valid), so the two centres and radii decide — concentric or not.
///
/// **Written once because it is asked from two directions**: the arrangement's nesting
/// ([`cell_in_cell`], a pin's trace under a boss's cap) and the coplanar merge (a region whose
/// outer bound is a circle asking which circle holes it owns).
pub(crate) fn disk_in_disk(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    a: &nacre_topo::CylinderDef,
    b: &nacre_topo::CylinderDef,
) -> Result<bool, BoolError> {
    let undecided = || reject(RejectReason::WitnessNotRational);
    let pa = combinatorics::circle_centre_rat(jd, wc, a).ok_or_else(undecided)?;
    let pb = combinatorics::circle_centre_rat(jd, wc, b).ok_or_else(undecided)?;
    // **Radii decide only where both sections are circles**. On an oblique class each
    // is an ellipse with semi-minor `r` and semi-major `r/|cos θ|`, and comparing radii is then
    // wrong in both directions: it says "not inside" for a pair separated **along the major axis**
    // that really nests (`ra 1`, `rb 2`, offset `4/3` on the `(3,4,0)` axis: `3·(4/3)/5 + 1 ≤ 2`),
    // and it happens to be right when the offset is concentric or along the minor axis. Refusing
    // the whole oblique population buys the first and sells the second — and the sold half is the
    // trade this kernel's DNA asks for (`honest-reject > silent-wrong`).
    // Asked **before** the `ra >= rb` return below: that return is sound for ellipses too (the
    // semi-minor is exactly `r`, so `ra < rb` is necessary for nesting), so putting the guard
    // after it would leave half the population unguarded. And asked **after** the two centres, so
    // the class is known to have coefficients — no unreachable "what if it has none" arm.
    let coeffs = jd.world_rat(wc).ok_or_else(undecided)?;
    let n = [coeffs[0], coeffs[1], coeffs[2]];
    if !combinatorics::class_carries_circle(&n, &a.dir())
        || !combinatorics::class_carries_circle(&n, &b.dir())
    {
        return Err(reject(RejectReason::ObliqueCircleClass));
    }
    let (ra2, rb2) = (a.r2(), b.r2());
    if ra2 >= rb2 {
        return Ok(false);
    }
    // `dist < r_b − r_a` with the radii as squares: the difference is never formed — the scalar
    // door reads `√rb² > √dist² + √ra²` through the root-sum identity, in integers. Both centres
    // lie on the class plane the axis is normal to, so the distance from `pb` to `a`'s axis is
    // the distance between the centres, which is what this always compared.
    Ok(
        nacre_exact::cylinders_nested(&pa, &a.dir(), ra2, &pb, rb2)
            == nacre_exact::Orient::Negative,
    )
}

/// **Which cell a question is about** — a loop of the class's arrangement, or a whole disk.
///
/// The two roads that ask this question hold their cells differently (the arrangement by index
/// into its walk's rings and circles, the coplanar merge by `Bound`), so the engine takes neither:
/// it takes the shape itself.
#[derive(Clone, Copy)]
pub(crate) enum Cell<'a> {
    Ring(&'a [combinatorics::RingEdge]),
    Disk(&'a nacre_topo::CylinderDef),
}

/// **Where a witness is, and in what vocabulary** — the 2D sibling of [`combinatorics::Probe`],
/// which does the same job one dimension up (a component's witness carries a ray direction too,
/// which is why the two are not one type).
///
/// A `Named` witness needs no coordinates, and that is not a convenience: a **rotated** class has
/// no rational description, so its points have no rational coordinates at all — the ray road is
/// the only road there, and it runs on names.
#[derive(Clone, Copy, Debug)]
enum Where {
    Named([usize; 3]),
    Coord([nacre_exact::Rat; 3]),
}

/// **A point that can decide a nesting question — and how much it decides.**
///
/// **The distinction is the whole reason this is a type.** The two loops do not
/// cross, so a point of `a`'s **boundary** settles the question outright: `∂a` is connected and
/// misses `∂b`, so all of `∂a` is on the side that point is on. A point of `a`'s **interior**
/// settles something weaker — «this point is in `b`» — which is *not* the same claim, because
/// `b` may sit inside `a` and then `a`'s interior points are in `b` while `a` is not. Carried as
/// a hand-written second clause on one arm («and the ring is not inside the circle»), that is one
/// copy of the arm away from being left out; as a type, every witness says which kind it is.
enum Witness {
    On(Where),
    In(Where),
    /// A witness this cell *has* but whose value could not be formed exactly — checked-`Rat`
    /// overflow, a class with no rational description. Not the same as having none.
    Unformed,
}

/// What one witness said about the target.
enum Said {
    In,
    Out,
    /// This witness cannot rank itself against the target — it sits on the target's boundary, or
    /// every ray from it grazes. **The next witness is the remedy**, which is why this is not an
    /// error.
    Abstain,
    /// The value needed to ask could not be formed. The next witness may still answer, but if
    /// none does, the honest name is the arithmetic one and not "no clear ray".
    Unformed,
    /// The question reads a plane class the gate proved silent about the cylinders
    /// ([`crate::planes::WorkingCyl::is_silent`]) — a class with no world statement — and nothing
    /// here could answer it without one. Not [`Said::Unformed`]: the same arithmetic also fails on
    /// a `Wide` name or an overflow, and those are width, where this is the gate's own
    /// undecidedness carried to the place that met it.
    Silent,
    /// **The road cannot read the target** — the ray road against a ring with a pierce corner,
    /// cast from a plane with a world statement (no cylinder table there), or a two-edge lens. Not
    /// [`Said::Abstain`]: that is one witness's circumstance and the next witness of the same road
    /// is its remedy; here every witness whose name holds no class the gate passed silent would say
    /// the same, so the rest of those are skipped and the other road (coordinates) answers. A
    /// witness whose name holds one never says this: [`ask`] answers it by the ring's box. For the
    /// refusal's name it counts as an abstention.
    Declines,
}

impl Said {
    /// **What the ray road's refusal means for this witness** — the one place its errors are read.
    ///
    /// From a ray plane with a world statement the road asks a pierce corner's side with no
    /// cylinder table (`every_ray`'s `silent_table`), so from such a plane a ring **with a pierce
    /// corner** is unreadable: `PierceVertexUnnamed` there is the ring's, and the road declines.
    /// Without one, the same name comes from an arc whose ends do not settle it against the ray's
    /// plane — both on it, both on one side (none or two crossings), or one on it (none or one):
    /// on that plane the road has no cylinder to side or count the arc with (a plane the gate
    /// passed silent answers from the gate's record — by boxes where it sections the cylinder).
    /// That is a plane the witness chose, so the next witness
    /// may answer, and it is an abstention. A ring of fewer than three edges is a lens
    /// when it is mixed (a chord and its arc — an ordinary shape the ray road cannot walk) and a
    /// producer defect when it is all straight, which is raised under its own name. Anything else
    /// (`RingNaming`: a ring whose names the road cannot read) is not a witness's circumstance and
    /// is raised rather than retried on the next witness.
    fn of_ray_error(e: BoolError, has_pierce: bool, mixed: bool) -> Result<Said, BoolError> {
        let BoolError::Rejected { reason, .. } = &e else {
            return Err(e);
        };
        match reason {
            RejectReason::NoClearRay | RejectReason::PointOnRing => Ok(Said::Abstain),
            RejectReason::PierceVertexUnnamed if has_pierce => Ok(Said::Declines),
            RejectReason::PierceVertexUnnamed => Ok(Said::Abstain),
            RejectReason::DegenerateRing if mixed => Ok(Said::Declines),
            _ => Err(e),
        }
    }
}

/// **A circle's own witnesses: four points on its rim, then its centre.**
///
/// The rim points are `centre ± r·û₁` and `centre ± r·û₂` over the cylinder's rational unit
/// cross-section frame ([`nacre_exact::cyl_unit_frame`]). They are **boundary** witnesses, and
/// that is the whole point: a boundary witness that lands strictly inside settles the question at
/// once, where the centre — an interior point — has to ask the converse, and a disk **cannot
/// answer** a converse (`inside`'s `may_ask_interior` skips its one witness, so the reverse
/// descent has always ended in [`crate::RejectReason::RingHasNoWitness`]).
///
/// **This is not new geometry — it is the point the kernel already mints.** For a whole circle
/// `crate::construct` states `ref_dir` as `(vertex − centre)/radius`, so `centre + radius·ref_dir`
/// **is** that ring's one vertex, exactly, by construction — the point a seam vertex
/// (`Vertex::OnSeam`) names. What the arrangement drops is the *node*: an
/// uncut circle carries no `NodeId`, which is why the corner supply above finds nothing and why
/// this rebuilds the point from the statement instead of reading it.
///
/// **Four, and no pair is droppable.** The parity ray runs along the class chart's `e₁`, and a
/// rim point shares the centre's ray line exactly when `û ∥ e₁`. `û₁ ⊥ û₂` span the plane ⊥ the
/// axis, so **at most one** of the two pairs can lie on that line — which also means at least one
/// pair always leaves it. Dropping either pair would leave a class with no fresh ray. (The next
/// rung, if four ever graze, is a Pythagorean direction — `((a²−b²)û₁ + 2ab·û₂)/(a²+b²)` is exact
/// and unit with no new square root — and it buys a **new ray line**, which is the thing that is
/// scarce, not a new radical.)
///
/// **Soundness rests on `û·û = 1` exactly**, not on the frame being cheap:
/// `cylinder_radial_side` reduces to `‖m‖²·r²·(û·û − 1)` for any `‖m‖`, so an inexact unit vector
/// would make a rim point a `Witness::On` that is not on the boundary — and `inside` returns from
/// a boundary `In` **without** the converse. Hence the exact `inv_sqrt_rat` road and the
/// post-conditions below, asserted where the point is made rather than where it is used
/// ([`combinatorics::conjugate_midpoint`]'s rule).
///
/// **The premise that makes these boundary points is guarded, not assumed.** A rim point
/// lies on the *cell* — not merely on the cylinder — because the class plane is ⊥ the axis:
/// `circle_centre_rat` answers for a tilted axis too (it declines only on `n·m = 0`), and
/// `û ⊥ m` does not give `û ⊥ n` unless `m ∥ n`. An `On` that is not on the boundary returns
/// from `inside` **without** the converse, which would be a silent wrong answer rather than a
/// refusal, so the rim is offered only to a class that carries the circle
/// ([`combinatorics::class_carries_circle`]). On an oblique class
/// `ask`'s radial arm is **right** (it reads the solid cylinder, and an
/// elliptical section is exactly the plane's points within `r` of the axis), while `disk_in_disk`
/// is wrong there — not for the centre distance (both
/// centres lie *on* the plane) but for the **radii**, which are not an ellipse's width.
///
/// A centre that cannot be formed keeps its own name ([`Witness::Unformed`] →
/// [`crate::RejectReason::WitnessNotRational`]); a frame that cannot be formed suppresses the
/// **rim only**, because "this class has no rational description" and "this cylinder's frame has
/// an irrational norm" are different sentences and only the first is that reason's.
fn rim_and_centre<'a>(
    jd: &Judge<'a, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> Vec<Witness> {
    let Some(centre) = combinatorics::circle_centre_rat(jd, wc, def) else {
        return vec![Witness::Unformed];
    };
    let mut out = Vec::with_capacity(5);
    // **The rim is a boundary point only where the class is ⊥ to the axis**. Off that
    // the section is an ellipse: `û ⊥ axis` does not give `û ⊥ n`, so `centre ± r·û` leaves the
    // class plane and a `Witness::On` there is a lie `inside` answers **without the converse** —
    // a silent wrong answer, not a refusal. The centre stays: it is the ellipse's centre either
    // way, and `ask`'s radial arm reads the solid cylinder, which is right for any section.
    // The guard is total ([`combinatorics::class_carries_circle`]) and fires nowhere today; the
    // gate is what empties its population, and this is what says so where the point is made.
    // `is_some_and` folds a *second* cause in — a class with no rational description also loses
    // its rim. That population is empty too: such a class passes the gate only as silent about
    // every cylinder (`WorkingCyl::is_silent`), and a silent class carries no circle. If it ever stops
    // being empty the cell falls to `RingHasNoWitness` rather than to a wrong answer, which is the
    // direction this kernel takes.
    let carries = jd
        .world_rat(wc)
        .is_some_and(|c| combinatorics::class_carries_circle(&[c[0], c[1], c[2]], &def.dir()));
    // Rim witnesses at `centre + r·u` need the radius itself; a cell whose squared radius has no
    // rational root is simply not offered them — it keeps its centre, exactly as when a witness's
    // arithmetic leaves `i128` below.
    if carries
        && let Some((u1, u2)) = nacre_exact::cyl_unit_frame(&def.dir(), &def.ref_dir())
        && let Some(r) = def.radius_rat()
    {
        let neg = |v: &[nacre_exact::Rat; 3]| -> Option<[nacre_exact::Rat; 3]> {
            let z = nacre_exact::Rat::from_int(0);
            Some([
                z.checked_sub(v[0])?,
                z.checked_sub(v[1])?,
                z.checked_sub(v[2])?,
            ])
        };
        let at = |u: &[nacre_exact::Rat; 3]| -> Option<[nacre_exact::Rat; 3]> {
            let mut p = centre;
            for k in 0..3 {
                p[k] = p[k].checked_add(r.checked_mul(u[k])?)?;
            }
            Some(p)
        };
        for d in [Some(u1), neg(&u1), Some(u2), neg(&u2)] {
            // A point whose arithmetic left `i128` is simply not offered: the cell still has its
            // centre, and calling it `Unformed` would rename the refusal.
            let Some(p) = d.and_then(|d| at(&d)) else {
                continue;
            };
            debug_assert_eq!(
                nacre_exact::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.r2()),
                nacre_exact::Orient::Zero,
                "a rim witness is on its own rim"
            );
            // The guard above makes this unreachable **for its own cause** — a class that
            // is not ⊥ never gets here. It stays because it also catches a broken frame or a
            // centre that is not the axis' meet, which no guard above asks about.
            debug_assert!(
                jd.world_rat(wc).is_none_or(|c| {
                    let n = [c[0], c[1], c[2]];
                    nacre_exact::dot3_rat(&n, &p)
                        .and_then(|d| d.checked_add(c[3]))
                        .is_none_or(|v| v == nacre_exact::Rat::from_int(0))
                }),
                "a rim witness is on the class it is a witness of"
            );
            out.push(Witness::On(Where::Coord(p)));
        }
    }
    out.push(Witness::In(Where::Coord(centre)));
    out
}

/// **How many boundary witnesses a disk offers** — the count alone, so a test can read it without
/// [`Witness`] and [`Where`] leaving this module. `circle_centre_rat` lives in `combinatorics`
/// rather than here for exactly that reason (*"opening that module's witness atoms … the shape
/// those atoms were made private to prevent"*), and a test's lock is not worth undoing it.
#[cfg(test)]
pub(crate) fn rim_witness_count(
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    def: &nacre_topo::CylinderDef,
) -> usize {
    rim_and_centre(jd, wc, def)
        .iter()
        .filter(|w| matches!(w, Witness::On(_)))
        .count()
}

/// **Every witness `cell` can offer, in one order, once.**
///
/// Boundary witnesses first and in this order: the corners' three-plane **names** (the ray road,
/// and the only road a rotated class has), then the corners' rational coordinates, then the
/// rational **pierce** corners (a fillet tangency), then a whole chord's **midpoint**,
/// then a rational point **inside** an edge whose ends are two solves' roots. An interior witness
/// last, and only where there is one: a disk's centre, and a ring that *is* a circle.
///
/// **A circle has no corner — and the conclusion drawn from that is false.**
/// Every supply above reads the *arrangement*: names, nodes, chords, edge interiors. An uncut
/// circle carries none of those — it arrives here as a `CylinderDef` and nothing else, and the
/// contributions merged onto it upstream state no `NodeId` — which is
/// true — and it is tempting to take the next step and conclude that such a cell has no
/// boundary point that can be named at all. That step is wrong. A circle's boundary points are
/// **geometric**, and they are exactly rational: [`rim_and_centre`].
///
/// **Lazy on purpose.** The corpus asks this question 33,781 times in one census pass, and every
/// kind after the first is a `pierce_meet` solve. The name that decides is almost always the first
/// one offered (measured: 23,581 of 24,449 ring questions), so nothing after it should be built.
fn witnesses<'a>(
    jd: &'a Judge<'a, WorkingPlane>,
    cyls: &'a [crate::planes::WorkingCyl],
    wc: usize,
    cell: Cell<'a>,
) -> Box<dyn Iterator<Item = Witness> + 'a> {
    // A disk's interior witness, wherever the disk came from — the cell itself, or a ring that
    // turned out to be a whole circle. Written once: this module is where each supply has one
    // spelling.
    let centre_of =
        |def: &nacre_topo::CylinderDef| match combinatorics::circle_centre_rat(jd, wc, def) {
            Some(c) => Witness::In(Where::Coord(c)),
            None => Witness::Unformed,
        };
    match cell {
        // A disk's boundary is a circle, and a circle's rim points are exact — the
        // supply is not the centre alone.
        // The wrapper does **not** buy laziness here and is not there for it: `inside` is the
        // only caller and its `for` always polls once, so the closure always runs. It is the
        // shape the ring arm needs when the two supplies fold into one, kept the same on both
        // sides so that fold is a move and not a rewrite. What this arm costs is five points
        // where the centre alone is one — a cost that has not been measured.
        Cell::Disk(def) => {
            Box::new(std::iter::once_with(move || rim_and_centre(jd, wc, def)).flatten())
        }
        Cell::Ring(r) => {
            let named = r.iter().filter_map(|e| {
                combinatorics::three_plane_name(e.node).map(|t| Witness::On(Where::Named(t)))
            });
            let coords = r.iter().flat_map(move |e| {
                combinatorics::edge_witness_points(jd, cyls, e)
                    .map(|p| Witness::On(Where::Coord(p)))
            });
            // A ring that is a whole circle does have boundary points anyone can name
            // exactly. It has four — the same rim
            // points [`rim_and_centre`] gives a `Cell::Disk`. This arm still offers only the
            // centre, although the two arms of one supply reading differently is exactly what
            // this module exists to prevent.
            //
            // **Deliberately: the answer to «give it the rim» is «do not».** The two numbers
            // that look like a contradiction — `ring_own_circle`'s 88 acceptances against *zero*
            // whole-circle rings as an engine question's source over 27,000 — are **different
            // populations** (85% of the probe's rows are not engine questions at all).
            // Both hold, and together they say a ring that *is* a whole circle never sources an
            // engine question here: the rim would be built for nobody. Its centre answers
            // meanwhile, as an interior witness, and that is enough for every row measured.
            // And it is eager: `Option::map` runs `centre_of` before this iterator is polled,
            // so the module's own "lazy on purpose" rule is already broken here. Folding both
            // arms onto `rim_and_centre` fixes that in the same move.
            let centre = ring_own_circle(r).map(centre_of).into_iter();
            Box::new(named.chain(coords).chain(centre))
        }
    }
}

/// **What this witness says about that cell** — the one place a point is ranked against a target.
fn ask(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    w: &Where,
    target: Cell<'_>,
) -> Result<Said, BoolError> {
    let disk = |side: nacre_geom::intersect::RingSide| match side {
        nacre_geom::intersect::RingSide::Inside => Said::In,
        nacre_geom::intersect::RingSide::Outside => Said::Out,
        // On the rim: this witness says nothing, the next may — one reading for the ring and
        // the disk alike.
        nacre_geom::intersect::RingSide::OnBoundary => Said::Abstain,
    };
    Ok(match (w, target) {
        (Where::Named(t), Cell::Ring(rb)) => {
            #[expect(
                clippy::disallowed_methods,
                reason = "one witness's answer; the retry over witnesses is `cell_inside`'s loop"
            )]
            let inside = combinatorics::point_in_ring(jd, cyls, wc, *t, rb);
            match inside {
                Ok(true) => Said::In,
                Ok(false) => Said::Out,
                // What the road's refusal means is read in one place ([`Said::of_ray_error`]); only
                // this witness's circumstance is an abstention, and anything else is raised.
                Err(e) => match Said::of_ray_error(
                    e,
                    rb.iter()
                        .any(|x| combinatorics::pierce_name(x.node).is_some()),
                    combinatorics::ring_is_mixed(rb),
                )? {
                    // A corner on a class with no world statement has no rational coordinates,
                    // so no coordinate witness stands behind it when the ray road cannot place it.
                    // The one question its interval answers is whether it stands outside the
                    // ring's box — the coordinate arm's own question below, asked of the same box.
                    Said::Declines | Said::Abstain if crate::planes::silent_among(cyls, t) => {
                        match jd.three_hp(*t, jd.standard.prec) {
                            Some(x) if outside_ring_box(jd, cyls, wc, &x, rb) => Said::Out,
                            _ => Said::Silent,
                        }
                    }
                    said => said,
                },
            }
        }
        (Where::Coord(p), Cell::Ring(rb)) => match rational_point_in_ring(jd, cyls, wc, p, rb) {
            Ok(Some(true)) => Said::In,
            Ok(Some(false)) => Said::Out,
            Ok(None) => Said::Abstain,
            // A **value** that could not be formed — the class has no rational description. Not
            // an abstention: no coordinate witness will do better on this class. A ring whose
            // corners a silent class names has no rational corners, and is still asked whether
            // the point stands outside the ring's box — the one question intervals answer.
            Err(_) if combinatorics::names_a_silent_class(cyls, rb) => {
                let prec = jd.standard.prec;
                if outside_ring_box(
                    jd,
                    cyls,
                    wc,
                    &p.map(|x| nacre_exact::HpBounded::of_rat(x, prec)),
                    rb,
                ) {
                    Said::Out
                } else {
                    Said::Silent
                }
            }
            Err(_) => Said::Unformed,
        },
        (Where::Named(t), Cell::Disk(def)) => {
            // **A corner on a class the gate proved apart from this disk's cylinder is outside
            // the disk, exactly** — here, and not in the shared reading below: the corner lies on
            // a face of that class (a class's segments are its faces met by the cell's plane),
            // those faces stand apart from the lateral boxes, and the disk lies in the box its rim
            // does. A `Sectioned` class never gets here: its cylinder has no whole rim, so no disk.
            // The argument reads «a corner named by the class lies on one of its faces», which a
            // cell's corner keeps and a point off the class's faces does not; an alias that put the
            // class into the name of a point off its faces would break it unseen. A class proved
            // `Beside` the cylinder needs no such premise and is read by the shared rule.
            let apart = cyls
                .iter()
                .find(|cy| cy.def == *def)
                .is_some_and(|cy| t.iter().any(|&k| cy.is_silent(k) && cy.beside(k).is_none()));
            if apart {
                return Ok(Said::Out);
            }
            match combinatorics::named_point_in_disk(jd, cyls, *t, def) {
                Some(side) => disk(side),
                None => Said::Unformed,
            }
        }
        (Where::Coord(p), Cell::Disk(def)) => disk(combinatorics::point_in_disk(p, def)),
    })
}

/// **Is `p` outside the box that holds the ring?** — then it is outside the ring. The box is
/// [`combinatorics::ring_box_points_hp`]'s — every corner whatever names it, and every arc's own
/// box — in the judge's interval arithmetic, and `p` is an interval too (a coordinate witness exactly,
/// a corner on a class with no world statement as its witness triangles' meet). A comparison whose
/// interval straddles zero proves nothing, and the answer is then `false`.
fn outside_ring_box(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    p: &[nacre_exact::HpBounded; 3],
    rb: &[combinatorics::RingEdge],
) -> bool {
    let Some(points) = combinatorics::ring_box_points_hp(jd, cyls, wc, rb) else {
        return false;
    };
    let prec = jd.standard.prec;
    (0..3).any(|k| {
        let signs: Vec<Option<bool>> = points
            .iter()
            .map(|x| p[k].sub(&x[k], prec).sign())
            .collect();
        signs.iter().all(|s| *s == Some(true)) || signs.iter().all(|s| *s == Some(false))
    })
}

pub(crate) fn cell_inside(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    a: Cell<'_>,
    b: Cell<'_>,
) -> Result<bool, BoolError> {
    inside(jd, cyls, wc, a, b, true)
}

fn inside(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    a: Cell<'_>,
    b: Cell<'_>,
    // **The converse is asked one level deep and no further.** An interior witness of `a` that
    // lands in `b` leaves the question open, and `b`'s **boundary** witnesses settle it; if `b`
    // has none, the honest answer is the refusal below — not a second descent, which two cells
    // that each offer only an interior point (two rings that are both whole circles) would turn
    // into an infinite one.
    may_ask_interior: bool,
) -> Result<bool, BoolError> {
    // **The list is one list, and it is walked to the end.** Two road-choices could cut it
    // short, and both are statements about a *road* rather than about the question:
    //   · skipping the ray road for a **mixed** target in advance — the ray road says so itself
    //     when asked ([`Said::Declines`]), and from then on its remaining named witnesses are
    //     skipped where they would put the same question to the same road: every one whose name
    //     holds no class the gate passed silent (measured over the suite, the ignored sweep, the
    //     census and the kit: of 39,040 skips, the 38,992 without a silent class would each have
    //     declined again, and the 48 with one, which are asked, each answered);
    //   · answering a plain ring target on names **alone** once it names any probe, which refuses
    //     a ring whose every ray grazed with coordinates still in hand.
    // Neither is a fact about *whether `a` is inside `b`*, and the premise that lets them go is
    // measured, not assumed: where both roads can answer they agree
    // (`the_two_roads_never_disagree`).
    let (mut asked, mut unformed, mut declined, mut silent) = (false, false, false, false);
    #[cfg(test)]
    let mut tally = nesting_probe::DeclineTally::default();
    for w in witnesses(jd, cyls, wc, a) {
        let (p, interior) = match w {
            Witness::Unformed => {
                unformed = true;
                continue;
            }
            Witness::On(p) => (p, false),
            Witness::In(p) => {
                if !may_ask_interior {
                    continue;
                }
                (p, true)
            }
        };
        // `ask` answers a name holding a class the gate passed silent by the ring's box, so only
        // the others would repeat the road's refusal.
        let repeats =
            declined && matches!(&p, Where::Named(t) if !crate::planes::silent_among(cyls, t));
        if repeats {
            continue;
        }
        #[cfg(test)]
        if repeats {
            tally.named_after += 1;
        }
        #[cfg(test)]
        if declined && matches!(p, Where::Named(_)) {
            tally.silent_after += 1;
        }
        match ask(jd, cyls, wc, &p, b)? {
            Said::Unformed => unformed = true,
            Said::Silent => silent = true,
            Said::Abstain => asked = true,
            Said::Declines => {
                asked = true;
                declined = true;
                #[cfg(test)]
                {
                    tally.declines += 1;
                }
            }
            // A point of `a` outside `b`: whether it is a boundary point or an interior one, `a`
            // is not inside `b`.
            Said::Out => return Ok(false),
            Said::In if !interior => return Ok(true),
            Said::In => return Ok(!inside(jd, cyls, wc, b, a, false)?),
        }
    }
    Err(reject(if silent {
        RejectReason::CylinderGateUndecided
    } else if unformed {
        RejectReason::WitnessNotRational
    } else if asked {
        RejectReason::NoClearRay
    } else {
        RejectReason::RingHasNoWitness
    }))
}

/// **The arrangement's adapter onto [`cell_inside`]** — its cells arrive as indices into the
/// walk's rings and circle defs, and two questions are settled before a witness is asked.
///
/// **The indices are the caller's; the shapes are not.** [`Cell`] takes the shape itself
/// rather than either road's container, and this adapter follows the same rule one level up: it
/// is handed the exact statement of each circle, never the vessel the arrangement merged it in.
///
/// `Ok(None)` is *not comparable*: two polygons that **share a node** touch rather than nest (a
/// shared node is a split point, and that also excludes a contour's own `+1` partner, which
/// carries the same ring). It is a **ring-ring** rule and stays one: the disk arms have never
/// asked it, and the coplanar merge — the engine's other road — works with rings that do share
/// nodes. Whether the asymmetry is right is an open question of its own.
///
/// **The engine's precondition holds here structurally**: these cells are faces of one class's
/// DCEL, and the split gave every crossing a vertex — so two of their loops cannot cross, they can
/// only nest or share nodes.
///
/// **Asked from two directions here and a third elsewhere.**
/// `crate::arrangement::cells::nest_cells` asks it of (contour, `+1` cell) to find hosts and
/// `innermost_host` of (host, host) to order them; the coplanar merge asks the same question of its
/// own bounds, through this engine rather than a copy of this dispatch.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cell_in_cell(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    wc: usize,
    rings: &[Vec<combinatorics::RingEdge>],
    circle_defs: &[&nacre_topo::CylinderDef],
    circle_ix: &[Option<usize>],
    a: usize,
    b: usize,
) -> Result<Option<bool>, BoolError> {
    // **Adjacency stays a ring–ring rule.** A shared node is a split
    // point, so the two loops touch rather than one wrapping the other (that also excludes a
    // contour's own `+1` partner, which carries the same ring). It is *not* asked of the disk arms
    // today, and hoisting it into the engine would answer «not comparable» where they answer — and
    // would change the merge road, whose rings do share nodes. The asymmetry is left where it is
    // and written down; whether it is right is an open question.
    // **The route is one named decision.** Spelled as three fall-throughs with the instrument
    // below sitting above all of them, **85% of its rows describe questions no witness is asked
    // for** (measured: 189,662 rows against 28,000 engine questions over the lib suite). A row
    // that says "the source had N witnesses" about a
    // question answered by a shared node, or by two radii, is a lie the audit tests then read.
    let route = if circle_ix[a].is_none() && circle_ix[b].is_none() {
        let (ra, rb) = (&rings[a], &rings[b]);
        if ra.iter().any(|e| rb.iter().any(|f| f.node == e.node)) {
            Route::SharedNode
        } else {
            Route::Engine
        }
    } else if let (Some(ca), Some(cb)) = (circle_ix[a], circle_ix[b]) {
        Route::DiskPair { ca, cb }
    } else {
        Route::Engine
    };
    #[cfg(test)]
    if nesting_probe::on() {
        let inv = |i: usize| -> (usize, usize, usize, usize, usize, bool) {
            if circle_ix[i].is_some() {
                return (0, 0, 0, 0, 0, false);
            }
            let r = &rings[i];
            (
                r.iter()
                    .filter(|e| combinatorics::three_plane_name(e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::node_coords_rat(jd, e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::pierce_coords_rat(jd, cyls, e.node).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::conjugate_midpoint(jd, cyls, e).is_some())
                    .count(),
                r.iter()
                    .filter(|e| combinatorics::pierce_ends_between(jd, cyls, e).is_some())
                    .count(),
                ring_own_circle(r).is_some(),
            )
        };
        let (named, coords, pierce, chord, edge, circle) = inv(a);
        // `inv` speaks a *ring's* vocabulary, so a disk reads all zeros there and
        // **no row would describe the arm whose whole supply is its own**. Its rim is its
        // supply, so the row says how many it offered.
        let rim = match circle_ix[a] {
            Some(c) => rim_and_centre(jd, wc, circle_defs[c])
                .iter()
                .filter(|w| matches!(w, Witness::On(_)))
                .count(),
            None => 0,
        };
        // Where both roads can answer, ask both — the premise under test.
        let roads = (|| {
            if circle_ix[a].is_some() || circle_ix[b].is_some() {
                return None;
            }
            let (ra, rb) = (&rings[a], &rings[b]);
            // Through the engine's own per-witness door, so the comparison is between the two
            // **roads** and not between two hand-written loops: the first name that decides
            // against the first coordinate that decides.
            let decided = |w: Where| match ask(jd, cyls, wc, &w, Cell::Ring(rb)) {
                Ok(Said::In) => Some(true),
                Ok(Said::Out) => Some(false),
                Ok(Said::Abstain | Said::Unformed | Said::Silent | Said::Declines) | Err(_) => None,
            };
            let by_name = ra
                .iter()
                .filter_map(|e| combinatorics::three_plane_name(e.node))
                .find_map(|t| decided(Where::Named(t)))?;
            let by_coord = ra
                .iter()
                .flat_map(|e| combinatorics::edge_witness_points(jd, cyls, e))
                .find_map(|p| decided(Where::Coord(p)))?;
            Some((by_name, by_coord))
        })();
        nesting_probe::push(nesting_probe::Row {
            a_disk: circle_ix[a].is_some(),
            b_disk: circle_ix[b].is_some(),
            b_mixed: circle_ix[b].is_none() && combinatorics::ring_is_mixed(&rings[b]),
            named,
            coords,
            pierce,
            chord,
            edge,
            circle,
            rim,
            route,
            roads,
        });
    }
    // **Two disks are not a witness question at all**: a disk lies inside another iff
    // its rim does, and the rims of two classes never meet, so the radii and the centre distance
    // decide it exactly. Answering `None` here leaves two disk cells unnested and *silently*
    // keeps both operands' caps: a pin stacked on a boss fuses into **two** untouched bodies.
    match route {
        Route::SharedNode => Ok(None),
        Route::DiskPair { ca, cb } => {
            disk_in_disk(jd, wc, circle_defs[ca], circle_defs[cb]).map(Some)
        }
        Route::Engine => {
            let cell = |i: usize| match circle_ix[i] {
                Some(c) => Cell::Disk(circle_defs[c]),
                None => Cell::Ring(&rings[i]),
            };
            cell_inside(jd, cyls, wc, cell(a), cell(b)).map(Some)
        }
    }
}

/// **Which road answers a nesting question** — named because the instrument has to say so, and
/// because a rule spelled at the pierce and again at the probe is two rules
/// (`nesting_probe::Row::route`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Route {
    /// Two loops that share a node touch rather than nest: not comparable, no witness asked.
    SharedNode,
    /// Two disks: radii and centre distance decide, exactly ([`disk_in_disk`]). It carries the
    /// two circle indices rather than letting the arm re-derive them — the decision already knows
    /// them, and a second `let … else { unreachable!() }` would be the rule spelled twice.
    DiskPair { ca: usize, cb: usize },
    /// The witness engine ([`cell_inside`]).
    Engine,
}

/// **What every nesting question is asked with** (test-only).
///
/// **A row is not a population, and for two reasons — both measured.**
/// 1. **Most rows are about questions no witness was asked for.** This pushes from
///    [`cell_in_cell`], which answers three ways ([`Route`]): a shared node, two radii, or the
///    engine. Sampled over the lib suite at the moment the 28,000th engine question ran, the
///    adapter had been entered **189,662** times — so roughly **six of every seven rows** describe
///    an offer nobody read. Filter on [`Row::route`] before counting anything,
///    and never read a raw row count as a population (the rule `reject_census` states for raises,
///    one layer in).
/// 2. **The merge road is invisible.** `boolean` calls [`cell_inside`] directly, so those
///    questions never pass here at all — sampled at **~1% of engine questions in the lib suite and
///    ~3% in census**, which is why this is a footnote rather than the headline.
///
/// **The defect this instrument watches for is invisible to a reject count**: a corpus that never
/// puts a **ring with no three-plane corner** against a **disk** never fires an arm whose witness
/// supply is truncated to that one kind. So the instrument measures the **population**: for each
/// question, what the source cell had to offer and which road could answer.
///
/// Recorded only for work a test asked to own (`ledger::owned`) — the road-agreement column runs
/// *both* roads where both are available, which is work no ordinary test should pay for.
#[cfg(test)]
pub(crate) mod nesting_probe {
    use crate::ledger::Ledger;

    /// One question, as it arrived.
    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Row {
        pub a_disk: bool,
        pub b_disk: bool,
        pub b_mixed: bool,
        /// `a`'s witnesses by kind. **The rule is one place**: `chord` and `edge` are
        /// two arms of [`combinatorics::edge_interior_points`] rather than two producers chained by
        /// four different spellings, and the breakdown is kept because a *count per kind* is not a
        /// spelling of the rule.
        pub named: usize,
        pub coords: usize,
        pub pierce: usize,
        pub chord: usize,
        pub edge: usize,
        pub circle: bool,
        /// A **disk's** supply, which the five ring counts above cannot describe: how many rim
        /// witnesses it offered. `0` for a ring, including one that is a whole circle.
        pub rim: usize,
        /// **Which road answered.** Only [`super::Route::Engine`] rows have a witness asked of
        /// them; the counts above describe an offer nobody read on the other two. 85% of the rows
        /// this probe pushes are of that kind (189,662 against 28,000 engine questions, lib suite),
        /// which is why a population read off them unfiltered is wrong.
        pub route: super::Route,
        /// Where **both** roads could answer, what each said — the direct measurement of the
        /// engine's premise, that any witness gives the same answer.
        pub roads: Option<(bool, bool)>,
    }

    /// One row per question, owned by the test whose work asked it.
    pub(crate) static ROWS: Ledger<Row> = Ledger::new();

    /// Whether to record here at all. A question is cheap; describing its **supply** is not, and
    /// the suite asks tens of thousands of them, so only a test that came for these rows pays.
    pub(crate) fn on() -> bool {
        crate::ledger::owned_thread()
    }

    pub(crate) fn push(r: Row) {
        ROWS.push(r);
    }

    /// Per question that the ray road declined: how many witnesses it declined for, how many
    /// named witnesses the skip in [`super::inside`] should have passed over were still asked
    /// (`0` is the skip), and how many named witnesses holding a class the gate passed silent were
    /// asked after the refusal. A row is pushed as the question leaves, by whichever exit.
    pub(crate) static DECLINES: Ledger<(usize, usize, usize)> = Ledger::new();

    #[derive(Default)]
    pub(crate) struct DeclineTally {
        pub(crate) declines: usize,
        pub(crate) named_after: usize,
        pub(crate) silent_after: usize,
    }

    impl Drop for DeclineTally {
        fn drop(&mut self) {
            if self.declines > 0 && on() {
                DECLINES.push((self.declines, self.named_after, self.silent_after));
            }
        }
    }
}

#[cfg(test)]
#[path = "tests/nesting.rs"]
mod tests;
