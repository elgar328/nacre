use super::*;
/// A recorded chord [`chord_stations`] cannot state: no cylinder statement, or a meet that is not
/// a pair of roots (a class with no narrow world name and no interval pair).
pub(super) struct ChordUnstatable;

/// **The two ends of the chord a wall class `wc` cuts on a circle of a planar face** (the face lies
/// in class `fc`, the circle rides cylinder `cyl`), as the two stations where the circle crosses
/// the class line ([`sweep_line`]) — a diameter when the wall runs through the axis, any chord
/// within the radius otherwise. The wall is one whose arrangement draws the cylinder's section: a
/// cut, a class with no world statement that sections the cylinder along its axis on its own
/// arrangement, or a named class oblique to the axis that draws its ellipse
/// ([`crate::planes::CylRecord::Sections`]) — the chord is part of that section, as the rulings and
/// the ellipse's arcs are, and ends where they do. `Ok(None)` for any other pair: the wall clears
/// the circle (the gate's plane test; for an oblique class, a line missing that rim circle), or
/// cuts an **irrational** chord this road cannot state yet — today's silence, and the same trigger
/// the rulings road reads.
///
/// **The record is asked before the wall's coefficients**, as [`rulings_on_class`] asks it: an
/// unrecorded pair says nothing whatever the wall's description.
///
/// **Written once because a disk's outer and a bored face's hole are the same circle asked the
/// same question.**
///
/// The record is the cylinder table's
/// ([`WorkingCyl::draws_rulings`](crate::planes::WorkingCyl::draws_rulings)), so a caller that
/// traces with no cylinder table draws no chord; the cylinder's statement comes from any face row
/// on its class — the same table `merge_circles` indexes.
pub(super) fn chord_stations(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    wc: usize,
    fc: usize,
    cyl: usize,
    cyls: &[crate::planes::WorkingCyl],
) -> Result<Option<[ClassLineStation; 2]>, ChordUnstatable> {
    let def_of = || {
        faces
            .iter()
            .zip(plane_ix)
            .find_map(|(row, ix)| match (row, ix) {
                (FaceRow::Cylinder(cf), ClassIx::Cyl(k)) if *k == cyl => cf.def.as_ref(),
                _ => None,
            })
    };
    let node = |id: NodeId| (LineEvent::Toggle, (id, combinatorics::EndPin::Cylinder));
    // **An oblique section's chord is the rim's own two crossings.** A circle a face carries is a
    // rim of a lateral face, and the class line meets it where the section crosses that rim —
    // the same two pierce names the tracer ends the ellipse's arcs at. A line missing the circle
    // leaves no chord; a touch the gate refused (`planes::rim_against_plane`), and the radius
    // premise below — a distance from the axis — means nothing for a plane that crosses the axis.
    if cyls.get(cyl).is_some_and(|cy| cy.sections_oblique(wc)) {
        let def = def_of().ok_or(ChordUnstatable)?;
        return match combinatorics::rim_crossings(jd, def, wc, fc, cyl) {
            Some(Some([lo, hi])) => Ok(Some([node(lo), node(hi)])),
            Some(None) => Ok(None),
            None => {
                debug_assert!(
                    false,
                    "an oblique section's class line crosses or misses a rim circle"
                );
                Err(ChordUnstatable)
            }
        };
    }
    if !cyls.get(cyl).is_some_and(|cy| cy.draws_rulings(wc)) {
        return Ok(None);
    }
    let def = def_of().ok_or(ChordUnstatable)?;
    // The record is the crossing statement: a recorded pair's plane runs within the radius, so `L`
    // enters the circle at one root and leaves at the other. The premise is the gate's, asserted
    // rather than re-derived — exactly for a cut, whose class has world coefficients; for a
    // sectioned class only as «the intervals do not say beyond», since the gate may have proved it
    // exactly before the class's turn, closer to the radius than the intervals part. (A *tangent*
    // plane has one root and is neither; see `planes::Tangency`.)
    #[cfg(debug_assertions)]
    match jd.world_rat(wc) {
        Some(w) => debug_assert_eq!(
            nacre_exact::point_plane_clearance_rat(&w, &def.origin(), def.r2()),
            nacre_exact::Orient::Negative,
            "a recorded pair's plane runs within the radius"
        ),
        None => debug_assert_ne!(
            jd.plane_cuts(wc, def),
            Some(false),
            "a sectioned class's plane runs within the radius"
        ),
    }
    // One meet for both: a cut pair's roots exactly, a sectioned one's as interval boxes — the same
    // two `NodeId::pierce` roots either way.
    let nacre_judge::cylinder::LineMeet::Pair(_) = jd.line_meet(wc, fc, def, false) else {
        return Err(ChordUnstatable);
    };
    Ok(Some([
        node(NodeId::pierce(
            wc,
            fc,
            cyl,
            crate::combinatorics::QuadRoot::Lo,
        )),
        node(NodeId::pierce(
            wc,
            fc,
            cyl,
            crate::combinatorics::QuadRoot::Hi,
        )),
    ]))
}

/// **The radial side of each rational corner of a cell** — what [`cell_side`] reads to watch the
/// disk-side rule. A pierce corner sits *on* the circle and says nothing, and so does a three-plane
/// corner that happens to land there.
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
/// the *first* corner: a cell can have corners on both sides of a circle it borders. A fillet
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
