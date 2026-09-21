use super::*;
/// Both operands' traces on plane class `wc`, merged into one `Trace` (segments keep their
/// `solid` tag).
pub(super) fn trace_on_class(
    input: &combinatorics::TraceInput,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    plane_ix: &[ClassIx],
    aliases: &Aliases,
) -> Trace {
    let mut out = Trace::default();
    for (side, which) in [SolidSide::A, SolidSide::B].into_iter().enumerate() {
        trace_one(
            &input.faces[side],
            which,
            wc,
            jd,
            cyls,
            faces,
            plane_ix,
            &input.crossings,
            aliases,
            &mut out,
        );
    }
    out
}

/// Test shims: trace straight from the two solids, deriving [`combinatorics::TraceInput`] on the
/// spot. Production derives it once per boolean (`trace_result_faces`) because it is the same for
/// every class; a test that traces a single class should not have to say so.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_on_class_of(
    model: &Model,
    a: Handle<Solid>,
    b: Handle<Solid>,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc_a: &combinatorics::EdgeFaces,
    inc_b: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
    crossings: std::collections::HashSet<(usize, usize)>,
) -> Trace {
    let input = combinatorics::trace_input(
        model,
        [(a, inc_a), (b, inc_b)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
        &[],
        crossings,
    );
    trace_on_class(&input, wc, jd, cyls, faces, plane_ix, &Aliases::default())
}

/// One solid only — the second operand slot is filled with the same solid, whose loops are
/// identical, and only `faces[0]` is read.
///
/// ★ It carries the cylinder table: a disk cap's chord is two pierce-pinned nodes of
/// the parity sweep, and ordering them along the line asks the cylinder's statement.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_one_of(
    model: &Model,
    solid: Handle<Solid>,
    which: SolidSide,
    wc: usize,
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[crate::planes::WorkingCyl],
    faces: &[FaceRow],
    surf_ix: &HashMap<Handle<Face>, usize>,
    inc: &combinatorics::EdgeFaces,
    plane_ix: &[ClassIx],
    crossings: std::collections::HashSet<(usize, usize)>,
    out: &mut Trace,
) {
    let input = combinatorics::trace_input(
        model,
        [(solid, inc), (solid, inc)],
        surf_ix,
        faces.len(),
        jd,
        plane_ix,
        cyls,
        crossings,
    );
    trace_one(
        &input.faces[0],
        which,
        wc,
        jd,
        cyls,
        faces,
        plane_ix,
        &input.crossings,
        &Aliases::default(),
        out,
    );
}

/// A coincident-merged arrangement edge on a plane class. When the cut plane is a solid's **cap
/// plane**, every side wall traces the same edge twice — once as the cap's boundary (`Seated`) and
/// once as the wall's top edge (`Transversal`) — because the cap's rim *is* the wall's top edge.
/// These are one geometric edge, and the arrangement (and its DCEL) needs them as one: two
/// coincident `(fp, s)` edges at a vertex would collapse in `angular_order`'s zero bucket.
///
/// The merge keeps the geometry (`wall`, `end`) as one and **preserves every contribution** rather
/// than deciding a single `kind`: which label rule a coincident edge follows (a cap-rim edge flips
/// only the below-bit, seated-style) is verified in the label brick, not guessed here.
#[derive(Clone, Debug)]
pub(crate) struct MergedSeg {
    pub wall: usize,
    /// Read by the next brick (crossings + split); kept here so the merged edge carries its
    /// geometry, not just its contributions.
    pub end: [NodeId; 2],
    /// What pins each endpoint on this edge's line — see [`Seg::end_h`] for the plane case and
    /// [`combinatorics::EndPin`] for why the arc split needed a second arm.
    pub end_h: [combinatorics::EndPin; 2],
    /// Every `(solid, kind)` that produced this one geometric edge. Length 1 when nothing was
    /// coincident.
    pub merged: Vec<(SolidSide, SegKind)>,
    /// The travel sense from `end[0]` to `end[1]`, stated by whoever cut the piece — see
    /// [`combinatorics::Carrier::Plane`]. `None` everywhere except a sub-segment the arc or ruling
    /// split cut, whose ends it already knows the order of.
    pub sense: Option<i8>,
}

/// Merge segments that are the **same geometric edge** — same `wall` and same endpoint-triple set
/// (direction-independent) — into one `MergedSeg`, collecting their contributions. Partial overlap
/// (same `wall`, *different* extent) is left for the per-wall interval overlay in
/// `split_at_crossings` to resolve into non-overlapping sub-segments with unioned contributions.
pub(super) fn merge_coincident(
    jd: &Judge<'_, WorkingPlane>,
    segs: &[Seg],
    wc: usize,
    aliases: &Aliases,
) -> Vec<MergedSeg> {
    // Key an edge by (wall, sorted endpoint pair) — both folded onto their canonical names first,
    // because two producers can describe one edge with a different wall *and* different endpoint
    // names when planes are concurrent, and it takes both folds for the two keys to coincide.
    let key = |s: &Seg| -> (usize, [NodeId; 2]) {
        let (mut a, mut b) = (aliases.canon_point(s.end[0]), aliases.canon_point(s.end[1]));
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        (aliases.canon_wall(wc, s.wall), [a, b])
    };
    let mut order: Vec<(usize, [NodeId; 2])> = Vec::new();
    let mut groups: HashMap<(usize, [NodeId; 2]), MergedSeg> = HashMap::new();
    for s in segs {
        let k = key(s);
        groups
            .entry(k)
            .or_insert_with(|| {
                order.push(k);
                {
                    // The canonical name of the line, so every producer on it agrees. A handle
                    // stays as recorded while the endpoint's name does: an aliased wall carries
                    // the *same* line, so a handle that pinned an endpoint there still pins it
                    // here. ★ When the **name** folds (a tangent corner's `Pierce` onto
                    // the `ThreePlane` of its planes with this class) the pin is derived again
                    // from the representative and the line — [`combinatorics::pin_for`], the one
                    // rule — because a pin is a fact about the name beside it, not luggage.
                    let wall = aliases.canon_wall(wc, s.wall);
                    let mut end = [NodeId::three_planes(Canon3::three([0, 1, 2])); 2];
                    let mut end_h = s.end_h;
                    for k in 0..2 {
                        end[k] = aliases.canon_point(s.end[k]);
                        if end[k] != s.end[k] {
                            if let Some(pin) = combinatorics::pin_for(jd, wc, wall, end[k]) {
                                end_h[k] = pin;
                            }
                        }
                    }
                    MergedSeg {
                        wall,
                        end,
                        end_h,
                        merged: Vec::new(),
                        sense: None,
                    }
                }
            })
            .merged
            .push((s.solid, s.kind));
    }
    // Deterministic order: first appearance.
    order
        .into_iter()
        .map(|k| groups.remove(&k).unwrap())
        .collect()
}
