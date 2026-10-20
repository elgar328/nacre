use super::*;
/// A recorded chord [`chord_stations`] cannot state: no cylinder statement, no coefficients for the
/// face's own class, or a meet that is not a pair of roots.
pub(super) struct ChordUnstatable;

/// **The two ends of the chord a recorded wall class `wc` cuts on a circle of a planar face** (the
/// face lies in class `fc`, the circle rides cylinder `cyl`), as the two stations where the circle
/// crosses the class line ([`sweep_line`]) — a diameter when the wall runs through the axis,
/// any chord within the radius otherwise. `Ok(None)` for a pair the gate did not record:
/// the wall clears the circle (the gate's plane test), or cuts an **irrational** chord this road
/// cannot state yet — today's silence, and the same trigger the rulings road reads.
///
/// **The record is asked before the wall's coefficients**, as [`rulings_on_class`] asks it: an
/// unrecorded pair says nothing whatever the wall's description, and a recorded one always has a
/// world name — the gate records only classes it could state.
///
/// **Written once because a disk's outer and a bored face's hole are the same circle asked the
/// same question.**
///
/// The cylinder's statement comes from any face row on its class — the same table `merge_circles`
/// indexes — so a caller that traces with no cylinder table (the test shims) is served too.
pub(super) fn chord_stations(
    jd: &Judge<'_, WorkingPlane>,
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    wc: usize,
    fc: usize,
    cyl: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
) -> Result<Option<[ClassLineStation; 2]>, ChordUnstatable> {
    use nacre_exact::quad::CylinderMeet;
    if !crossings.contains(&(wc, cyl)) {
        return Ok(None);
    }
    let w = combinatorics::class_coeffs_rat(jd, wc).ok_or(ChordUnstatable)?;
    let def = faces
        .iter()
        .zip(plane_ix)
        .find_map(|(row, ix)| match (row, ix) {
            (FaceRow::Cylinder(cf), ClassIx::Cyl(k)) if *k == cyl => cf.def.as_ref(),
            _ => None,
        })
        .ok_or(ChordUnstatable)?;
    // The record is the crossing statement: a listed pair's plane runs within the
    // radius, so `L` enters the circle at one root and leaves at the other. The premise is the
    // gate's, asserted rather than re-derived. (A *tangent* plane has one root and is never
    // listed; see `planes::Tangency`.)
    debug_assert_eq!(
        nacre_exact::point_plane_clearance_rat(&w, &def.origin(), def.r2()),
        nacre_exact::Orient::Negative,
        "a recorded pair's plane runs within the radius"
    );
    let Some(CylinderMeet::Pair { .. }) = combinatorics::class_line_meets_cylinder(jd, wc, fc, def)
    else {
        return Err(ChordUnstatable);
    };
    let node = |root| {
        (
            LineEvent::Toggle,
            (
                NodeId::pierce(wc, fc, cyl, root),
                combinatorics::EndPin::Cylinder,
            ),
        )
    };
    Ok(Some([
        node(nacre_topo::QuadRoot::Lo),
        node(nacre_topo::QuadRoot::Hi),
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
