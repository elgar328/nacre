use super::*;
/// **The lateral faces of every cylinder class, emitted from the chart** — the regions road
/// ([`regions::walk`]). Per class: the chart, its cells, each cell read off the
/// neighbouring classes ([`Chart::read_cell`]), then the walk. What is refused here by name:
/// a class with no row or rows of both solids, a θ order that cannot be formed, and a
/// **present cell whose chamber could not be read** (`CylinderGateUndecided` — two speaking
/// sides disagreed, or no side spoke); the walk names its own.
///
/// No sign is derived here: chambers come from [`Chart::read_cell`], the keep rule is
/// `read_cell::keep_for`, the ruling identity is [`Chart::ruling_name`], the rim's nodes are the
/// ones the cleaned plane faces hold (`rims` — the labels are still read against the split,
/// `Curved::split_rims`), and the winding is `seam_step`'s (`classify_cycles`).
pub(crate) fn emit_lateral(
    kind: crate::BoolKind,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
    rims: &crate::draft::HeldRims,
) -> Result<Vec<LocalFace>, BoolError> {
    let mut out: Vec<LocalFace> = Vec::new();
    for k in 0..cyls.len() {
        let def = &cyls[k].def;
        // One class is one solid's surface (`cyl_rows` refuses a class with no row; two solids
        // share no handles) — stated as the refusal `cyl_rows` names for the wiring failure.
        let mut sides = rows.iter().filter(|r| r.class == k).map(|r| r.side);
        let Some(side) = sides.next() else {
            return Err(reject(RejectReason::CylinderGateUndecided));
        };
        if sides.any(|s| s != side) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
        let chart = chart_of(jd, cyls, k, plane_faces, curved)?;
        let lines = Lines::of(jd, k, def, curved)?;
        // The θ order the split already formed cannot fail to form again; if it does, the point's
        // own `(line, s)` could not be stated — that name's sentence.
        let cells = chart
            .cells(jd, k, def)
            .ok_or_else(|| reject(RejectReason::WitnessNotRational))?;
        let reads: Vec<CellRead<'_>> = cells
            .iter()
            .map(|c| chart.read_cell(jd, k, def, kind, side, c, curved, &lines, rows))
            .collect::<Result<_, _>>()?;
        // A present cell whose chamber could not be read: two of its speaking sides disagreed,
        // or none spoke — `chamber`'s own refusal.
        if reads.iter().any(|r| r.present && r.emit.is_none()) {
            return Err(reject(RejectReason::CylinderGateUndecided));
        }
        let walked = regions::walk(
            jd, k, def, kind, side, &chart, &lines, &cells, &reads, curved, rims,
        )?;
        out.extend(walked.faces);
    }
    Ok(out)
}
