use super::*;
/// Why a circle's chord on a class could not be stated — [`chord_nodes`]' two refusals, which its
/// two callers read differently: a **circular outer** falls silent on missing coefficients (the
/// disk's old silence for a class it cannot read) and declines the rest, a **circular hole**
/// declines both.
pub(super) enum ChordFail {
    /// The class `wc` has no rational world coefficients (a rotated class).
    NoCoefficients,
    /// No cylinder statement, no coefficients for the face's own class, or the meet is not a pair
    /// of roots — a piece this road cannot state.
    Unstatable,
}

/// **The two ends of the chord a recorded wall class `wc` cuts on a circle of a planar face** (the
/// face lies in class `fc`, the circle rides cylinder `cyl`), as flip nodes of
/// [`trace_transversal_face`]'s parity sweep — a diameter when the wall runs through the axis,
/// any chord within the radius otherwise. `Ok(None)` for a pair the gate did not record:
/// the wall clears the circle (the gate's plane test), or cuts an **irrational** chord this road
/// cannot state yet — today's silence, and the same trigger the rulings road reads.
///
/// ★ **Written once because a disk's outer and a bored face's hole are the same circle asked the
/// same question.**
///
/// The cylinder's statement comes from any face row on its class — the same table `merge_circles`
/// indexes — so a caller that traces with no cylinder table (the test shims) is served too.
pub(super) fn chord_nodes(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    wc: usize,
    fc: usize,
    cyl: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
) -> Result<Option<[Node; 2]>, ChordFail> {
    use nacre_exact::quad::CylinderMeet;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(ChordFail::NoCoefficients)?;
    if !crossings.contains(&(wc, cyl)) {
        return Ok(None);
    }
    let def = faces
        .iter()
        .zip(plane_ix)
        .find_map(|(row, ix)| match (row, ix) {
            (FaceRow::Cylinder(cf), ClassIx::Cyl(k)) if *k == cyl => cf.def.as_ref(),
            _ => None,
        })
        .ok_or(ChordFail::Unstatable)?;
    // ★ The record is the crossing statement: a listed pair's plane runs within the
    // radius, so `L` enters the circle at one root and leaves at the other. The premise is the
    // gate's, asserted rather than re-derived. (A *tangent* plane has one root and is never
    // listed; see `planes::Tangency`.)
    debug_assert_eq!(
        nacre_exact::point_plane_clearance_rat(&w, &def.origin(), def.r2()),
        nacre_exact::Orient::Negative,
        "a recorded pair's plane runs within the radius"
    );
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(ChordFail::Unstatable)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let Some(CylinderMeet::Pair { .. }) =
        nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2)
    else {
        return Err(ChordFail::Unstatable);
    };
    let node = |root| Node {
        id: NodeId::pierce(wc, fc, cyl, root),
        pin: combinatorics::EndPin::Cylinder,
        flip: true,
        run: None,
        flanks_differ: false,
        single_touch: false,
        graze_above: None,
    };
    Ok(Some([
        node(nacre_topo::QuadRoot::Lo),
        node(nacre_topo::QuadRoot::Hi),
    ]))
}

/// **The planar scan's crossing on a ruling** — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a cylinder's ruling, named as the pierce node
/// `wc ∩ fc ∩ cyl` at the root that *is* this ruling.
///
/// A ruling edge lies in the face's own plane `fc` (a plane holding a ruling runs through the
/// axis — or is tangent, which the gate passes but does **not** record, so no ruling of it ever
/// reaches here), so the pair `{wc, fc}` cuts the cylinder in two
/// points, one on each of `fc`'s two rulings, and `(cyl, side)` — the identity
/// [`crate::combinatorics::Wall::Ruling`] carries, measured by
/// [`crate::combinatorics::ruling_side`] against `fc` when the
/// ring was named — says which. The same predicate asked of each root picks it. `wc` must be ⊥
/// to the axis for the class to cross a ruling in a point at all (∥ contains it; a tilt is
/// refused at the gate).
///
/// ★ Solved in the caller's order `(wc, fc)`, as [`rulings_on_class`] and [`chord_nodes`]
/// spell it — [`NodeId::pierce`] canonicalizes the pair and the root together, so this is the one
/// name every road gives the point. The lateral face names the same point when its hole ring
/// crosses the class ([`cycle_on_class`]'s crossing arm restates the hole corner's own root to
/// `wc` by the axis senses of the two ⊥ classes); the two spellings agree because the roots of
/// `{⊥, wall}` are ordered along `ε·k·(m̂ × n_wall)`, so «which root» and «which side of `wall`»
/// are the same question — a derivation the exact-volume rows of the crossing census check.
///
/// `Err(CurvedRingWall)` is every shape this does not state: no exact description, a wall that
/// is not parallel to the axis, no pair of roots, or both roots on one side (the two rulings of
/// a wall within the radius — through the axis or offset from it — are symmetric about
/// the plane through the axis with normal `m × n̂`, so `ruling_side` tells them apart).
pub(crate) fn crossing_on_ruling(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    side: i8,
) -> Result<NodeId, DeclineKind> {
    use nacre_exact::quad::CylinderMeet;
    let no = DeclineKind::CurvedRingWall;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(no)?;
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(no)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    if !nacre_exact::parallel_rat(&[w[0], w[1], w[2]], &m) {
        return Err(no);
    }
    // ★ A **tangent** wall (`side == 0`) has one ruling and one root — `Double`.
    let meet = nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2);
    if side == 0 {
        return match meet {
            Some(CylinderMeet::Tangent { .. }) => {
                Ok(NodeId::pierce(wc, fc, cyl, nacre_topo::QuadRoot::Double))
            }
            _ => Err(no),
        };
    }
    let Some(CylinderMeet::Pair { line, s }) = meet else {
        return Err(no);
    };
    let mut found = None;
    for (root, sv) in [
        (nacre_topo::QuadRoot::Lo, &s[0]),
        (nacre_topo::QuadRoot::Hi, &s[1]),
    ] {
        if combinatorics::ruling_side(&v, def, (&line, sv)) == Some(side) {
            if found.is_some() {
                return Err(no); // both roots on one side: not a pair of rulings
            }
            found = Some(root);
        }
    }
    Ok(NodeId::pierce(wc, fc, cyl, found.ok_or(no)?))
}

/// **The planar scan's crossing on an arc** — the point where the class line `L = wc ∩ fc`
/// leaves the face across an edge riding a circle of `cyl`, named as the pierce node
/// `wc ∩ fc ∩ cyl` at the root that lies **inside** the arc.
///
/// The arc lies in the face's own plane `fc` (a cap's ⊥ plane), so the pair `{wc, fc}` cuts the
/// cylinder in two points on the arc's circle, and the crossing is the one the travelled arc
/// `a → b` strictly contains ([`theta_between`], the ruling sweep's containment — the arc's CCW
/// pair is `ccw ? [a, b] : [b, a]`). Exactly one must: none is a crossing the walk mis-read, and
/// both is an arc meeting the line twice — a shape the gate keeps out (a > π arc against its own
/// diameter's class), refused rather than guessed.
///
/// ★ The same canonical name the lateral's ruling sweep gives this point as a station
/// (`crossing_on_ruling(fc, wc, …)` — [`NodeId::pierce`] folds the pair and root together), so
/// [`merge_coincident`] reads one point, not two.
#[allow(clippy::too_many_arguments)]
pub(crate) fn crossing_on_arc(
    jd: &Judge<'_, WorkingPlane>,
    def: &nacre_topo::CylinderDef,
    wc: usize,
    fc: usize,
    cyl: usize,
    ccw: bool,
    a: NodeId,
    b: NodeId,
) -> Result<NodeId, DeclineKind> {
    use nacre_exact::quad::CylinderMeet;
    let no = DeclineKind::CurvedRingWall;
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(no)?;
    let v = combinatorics::class_coeffs_rat(jd, fc).ok_or(no)?;
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    let Some(CylinderMeet::Pair { s, .. }) =
        nacre_exact::quad::plane_plane_cylinder(&w, &v, &o, &m, r2)
    else {
        return Err(no);
    };
    let _ = s;
    let (lo, hi) = if ccw { (a, b) } else { (b, a) };
    let mut found = None;
    for root in [nacre_topo::QuadRoot::Lo, nacre_topo::QuadRoot::Hi] {
        let node = NodeId::pierce(wc, fc, cyl, root);
        if theta_between(jd, cyl, def, lo, hi, node).map_err(|_| no)? {
            if found.is_some() {
                return Err(no); // both roots inside: the arc meets the line twice
            }
            found = Some(node);
        }
    }
    found.ok_or(no)
}

/// The scan's crossings on rulings, realized — the second road to the sign
/// [`crossing_on_ruling`] chooses by. One entry per crossing named in this binary: the point,
/// how far it sits from the class plane, the face plane and the cylinder (all should be 0), and
/// the ruling side read off the realization beside the side the ring carried.
/// **The disk-side rule's premise, watched where the rule is applied** (`emit_faces`' arc labels).
///
/// A cut circle's per-arc label must be the cell on the arc's **disk** side, and that side is read
/// off the cells themselves: a cell cannot straddle the circle — the circle is an arrangement edge
/// — so any rational corner of a cell with a definite radial side names the side the whole cell is
/// on. This watches the premise rather than the conclusion: **no cell has corners on both sides**,
/// and **an arc's two cells are never on the same side**. Both would make the answer a coin toss,
/// and neither is checkable from the label afterwards (a holed lateral's two sectors can carry a
/// literally identical label — see [`ArcLabel`]).
///
/// ★ It replaced a rule derived from the class's **stored** frame (`axis_up`), which two classes
/// with identical stored *and* canonical normals were measured to disagree about — the same plane,
/// the same circle, the same two arcs, opposite half-edges. That rule had been set by measuring
/// two fixtures and generalising; this one asks the geometry every time.
/// **Which side of a circle a cell lies on, from its own corners** — the geometry that watches
/// the disk-side rule. A cell cannot straddle the circle (the circle is an arrangement edge), so
/// any rational corner with a definite radial side names the side the whole cell is on; a pierce
/// corner sits *on* the circle and says nothing, and so does a three-plane corner that happens to
/// land there.
#[cfg(test)]
fn corner_sides<'a>(
    jd: &'a Judge<'_, WorkingPlane>,
    edges: &'a ClassEdges<'_>,
    cell: &'a Cell,
    def: &'a nacre_topo::CylinderDef,
) -> impl Iterator<Item = bool> + 'a {
    let (o, m, r2) = (def.origin(), def.dir(), def.r2());
    cell.half_edges.iter().filter_map(move |&h| {
        let p = combinatorics::node_coords_rat(jd, edges.origin(h))?;
        match nacre_exact::cylinder_radial_side(&p, &o, &m, r2) {
            nacre_exact::Orient::Negative => Some(true),
            nacre_exact::Orient::Positive => Some(false),
            nacre_exact::Orient::Zero => None,
        }
    })
}

/// **Which side of a circle a cell lies on, by its rational corners — when they agree.** Not by
/// the *first* corner: a cell can have corners on both sides of a circle it borders. ★ A fillet
/// does — the cap face's cell is bounded by a **quarter** of the circle and reaches far beyond
/// it, so its corners lie outside while the arc bounds it from the disk side. Such a cell has no
/// single side and the witness abstains; a cell all of whose corners agree answers.
#[cfg(test)]
pub(super) fn cell_side(
    jd: &Judge<'_, WorkingPlane>,
    edges: &ClassEdges<'_>,
    cell: &Cell,
    def: &nacre_topo::CylinderDef,
) -> Option<bool> {
    let mut sides = corner_sides(jd, edges, cell, def);
    let first = sides.next()?;
    sides.all(|s| s == first).then_some(first)
}
