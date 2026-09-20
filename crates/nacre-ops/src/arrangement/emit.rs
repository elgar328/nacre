use super::*;
/// Emit the result faces on plane class `wc` for a boolean `kind`. A `+1` cell is a face of the
/// result iff its two chambers disagree under `keep` (material on one side of W, void on the
/// other); its `-1` holes (`nesting.holes`) ride along as inner rings. The DCEL cell ring is
/// already CCW about `n_out(wc)` (the walk stored `winding == +1`) and a hole cell is CW
/// (`winding == -1`) — exactly the `LocalFace.inner` contract ("kept material on the loop's left"),
/// so both are emitted verbatim; `flip` alone carries the chamber and `assemble_fuse_cut` reverses
/// outer and inner together, making the result normal point out of the kept solid:
/// `flip = keep_above == (orient_sign(wc) > 0)`.
///
/// Only `+1` cells are hosts — a `-1` cell is either a hole (emitted as some host's inner ring) or
/// the void root — so `-1` cells are skipped, never emitted as their own face.
///
/// **Output contract:** every ring vertex is a [`combinatorics::NodeId`] triple, including triples that coincide
/// with an original A/B vertex (a cap corner); the weld table canonicalizes such a W-triple onto
/// the same result vertex, or `assemble_fuse_cut`'s manifold guard rejects. This brick emits
/// all-`Seam`; the reconciliation and assembly are later.
#[allow(clippy::too_many_arguments)]
pub(super) fn emit_faces(
    kind: BoolKind,
    labels: &[Label],
    cells: &[Cell],
    edges: &ClassEdges<'_>,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    holes: &HashMap<usize, Vec<usize>>,
) -> EmitOut {
    let planes = jd.planes;
    // [`combinatorics::NodeId`] is the vertex-identity enum; the local `Node` (this module's
    // three-valued-scan struct, a *different* type that happens to share the word) shadows the
    // glob-imported name here, which is why the identity is always spelled out in full.
    // ★ The wall travels with the ring. A half-edge was *told* which plane its edge rides, and
    // that is the one fact a name cannot always give back (see `draft::Ring`). A circle cell
    // (pseudo-half-edge past the segment range) has no nodes — its boundary is the cylinder
    // class itself.
    let bound_of = |cell: &Cell| -> crate::draft::Bound {
        if let Some(&he) = cell.half_edges.first()
            && let HalfEdgeKind::Circle(i) = edges.kind(he)
        {
            return crate::draft::Bound::Circle {
                cyl: edges.circles[i].cyl,
            };
        }
        // ★★ `Ring.walls` is a carrier now (`combinatorics::Wall`), so an arc half-edge has something
        // true to put here at last: its cylinder class and which way this side travels. What is
        // still missing is downstream — `edge_for` refuses an arc carrier by the population's
        // name until the arc-casting cell teaches it the ordered circle key.
        crate::draft::Bound::Ring(crate::draft::Ring::new(
            cell.half_edges.iter().map(|&he| edges.origin(he)).collect(),
            cell.half_edges
                .iter()
                .map(|&he| match edges.kind(he) {
                    HalfEdgeKind::Seg(i) => crate::combinatorics::Wall::Plane(edges.segs[i].wall),
                    // `edge_at`'s own convention, carried not re-derived: `MergedArc::end` runs
                    // counter-clockwise about the axis, so the even half-edge travels that way
                    // and its twin the other.
                    HalfEdgeKind::Arc(i) => crate::combinatorics::Wall::Arc {
                        cyl: edges.arcs[i].cyl,
                        ccw: he % 2 == 0,
                    },
                    // Same shape for a ruling: `MergedRuling::end` ascends the axis, the even
                    // half-edge travels up. `edge_for` refuses this wall by the population's
                    // name (`RulingBoundNotYet`) until the panel cell teaches it the key.
                    HalfEdgeKind::Ruling(i) => {
                        let r = &edges.rulings[i];
                        crate::combinatorics::Wall::Ruling {
                            cyl: r.cyl,
                            side: r.side,
                            up: he % 2 == 0,
                        }
                    }
                    // Mirrors `origin`'s statement: the nodes map above already refused it.
                    HalfEdgeKind::Circle(_) => {
                        unreachable!("a circle's pseudo-half-edge has no vertex to leave")
                    }
                })
                .collect(),
        ))
    };
    let mut out = Vec::new();
    for (c, cell) in cells.iter().enumerate() {
        if cell.winding != 1 {
            continue; // a -1 cell is a hole or the void root, never a face of its own
        }
        let l = labels[c];
        let keep_above = keep(kind, l[0], l[2]);
        let keep_below = keep(kind, l[1], l[3]);
        if keep_above == keep_below {
            continue; // material the same on both sides ⇒ not a result face here
        }
        let flip = keep_above == (planes[wc].frame_sign > 0);
        let inner: Vec<crate::draft::Bound> = holes
            .get(&c)
            .map(|hs| hs.iter().map(|&h| bound_of(&cells[h])).collect())
            .unwrap_or_default();
        out.push(LocalFace {
            surf: crate::planes::ClassIx::Plane(wc),
            outer: bound_of(cell),
            inner,
            flip,
        });
    }
    // ★★ **Every circle's disk label, whether or not a face was kept**. The four bits
    // of a disk cell say which solid's material lies immediately above and below this plane
    // *inside the circle* — which is exactly the chamber of the cylinder slab that starts here,
    // so the band pass reads them instead of casting a witness ray.
    //
    // ★ Collected from the **cells**, deliberately outside the keep filter above. A through
    // hole's outermost bands end on the cylinder's own cap classes, and `Cut` drops those cap
    // faces — following the faces would leave those bands with no label to read, while the cell
    // (and its label) is right here.
    let disk_labels = edges
        .circles
        .iter()
        .enumerate()
        .filter_map(|(i, mc)| {
            let he = edges.he_count() + 2 * i;
            let c = cells.iter().position(|cell| cell.half_edges == [he])?;
            #[cfg(test)]
            assert!(
                mc.whole_marks()
                    .ok()
                    .and_then(|w| disk_side_agrees(&w, &labels[c]).ok())
                    .unwrap_or(true),
                "a whole circle's disk does not carry its own solid: {:?}",
                labels[c]
            );
            Some((mc.cyl, labels[c]))
        })
        .collect();
    // ★ **A cut circle's labels, per arc** ([`ArcLabel`]) — same discipline as the disk labels
    // above: collected from the cells, outside the keep filter.
    //
    // ★★★★★ **Which half-edge borders the disk side, derived — and the factor that was missing.**
    // `MergedArc::end` runs counter-clockwise about the *axis*, so the even half-edge travels
    // `+θ̂`. The walk keeps a cell on the **left of its travel in the root face's frame**, whose
    // outward is `n_out = frame_sign · n_P` ([`crate::planes::WorkingPlane::frame_sign`] — the
    // sentence's one home; the ∥ road's `ruling_interior_is_even` cites it too), so
    //
    // ```text
    //   left = n_out × θ̂ = frame_sign·sign(n_P·m̂)·(m̂ × θ̂) = −frame_sign·axis_up·r̂
    //   ⇒ the even half-edge borders the disk  ⟺  axis_up · frame_sign = +1
    // ```
    //
    // The same product is already written one module over — `combinatorics`'s
    // `smooth_extremum_winding` computes `winding = ccw · axis_up · frame_sign`, and a circle
    // cell's boundary winds `+1` on the disk side, which is this statement rearranged.
    //
    // ★★★★★ **The same rule missing `frame_sign` is refuted.** Reading `axis_up` alone is
    // wrong. Measured: a boss cutting a plate and
    // the **notch it leaves** put the same plane (z = 1), the same stored normal and the same two
    // arcs on **opposite** half-edges — the boss's top cap is `Forward` there and the notch's
    // ceiling is `Reversed`, so `frame_sign` is the one thing that differs. The label that comes
    // back without it is the annulus cell's, all four bits true, which no disk-side cell there
    // can be. Over the lib suite the factor decides **112** of some 4,300 arcs, and without it
    // every one of them is wrong, silently: the only consumer is the chart's `read_cell`,
    // which mostly refuses
    // before reading a cut end.
    //
    // ★★ **The geometry still watches it, on every arc** ([`disk_side_probe`]): a cell cannot
    // straddle the circle, so a rational corner's radial side names the side that whole cell is
    // on — and where a corner speaks it must agree with this rule. Measured 3,260 of 3,260 on the
    // census corpus. That check is what makes this a *derivation* rather than a third guess.
    let ns_arcs = 2 * edges.segs.len();
    let arc_labels = edges
        .arcs
        .iter()
        .enumerate()
        .filter_map(|(i, ma)| {
            // ★ `plus_t_is_above` and nothing spelled beside it: the inline `normal().dot(axis)`
            // that used to stand here was that function's second spelling, letter for letter.
            let axis_up = crate::planes::plus_t_is_above(&jd.planes[wc], &ma.def);
            let even = ns_arcs + 2 * i;
            let he = even + usize::from(axis_up != (jd.planes[wc].frame_sign > 0));
            #[cfg(test)]
            disk_side_probe::record(
                cells
                    .iter()
                    .position(|q| q.half_edges.contains(&even))
                    .and_then(|cx| cell_side(jd, edges, &cells[cx], &ma.def)),
                cells
                    .iter()
                    .position(|q| q.half_edges.contains(&(even + 1)))
                    .and_then(|cx| cell_side(jd, edges, &cells[cx], &ma.def)),
                he == even,
                jd.planes[wc].frame_sign,
            );
            let c = cells
                .iter()
                .position(|cell| cell.half_edges.contains(&he))?;
            #[cfg(test)]
            assert!(
                disk_side_agrees(&ma.merged, &labels[c]).unwrap_or(true),
                "{}: {:?} {:?}",
                disk_side_probe::NOT_OWN_SOLID,
                ma.merged,
                labels[c]
            );
            // ★ The label and the trace that made it, **from one visit to one arc** — see
            // [`ArcLabel`] for why they may not become two maps.
            Some((
                ma.cyl,
                ArcLabel {
                    ends: ma.end,
                    label: labels[c],
                    marks: ma.merged.clone(),
                },
            ))
        })
        .collect();
    (out, disk_labels, arc_labels)
}

/// [`emit_faces`]' product: the kept faces, the disk labels, and the cut circles' per-arc
/// disk-side labels.
type EmitOut = (Vec<LocalFace>, Vec<(usize, Label)>, Vec<(usize, ArcLabel)>);

/// A face's box, from its vertices' cached coordinates.
///
/// **`f64` is the right precision here and that is not a compromise.** Nothing this decides is
/// part of an answer — the culling spike below only *counts* what a box would skip.
#[cfg(test)]
pub(super) fn face_box(model: &Model, fh: Handle<Face>) -> [[f64; 2]; 3] {
    let mut b = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    for p in face_points(model, fh) {
        for (k, s) in b.iter_mut().enumerate() {
            s[0] = s[0].min(p[k]);
            s[1] = s[1].max(p[k]);
        }
    }
    b
}

/// A face's vertex coordinates, from the caches.
///
#[cfg(test)]
fn face_points(model: &Model, fh: Handle<Face>) -> Vec<[f64; 3]> {
    let f = model.face(fh);
    let mut out = Vec::new();
    for lp in std::iter::once(&f.outer).chain(f.inner.iter()) {
        for he in &lp.half_edges {
            for &vh in model.edge(he.edge).vertices.iter() {
                out.push(model.vertex_point(vh).as_array());
            }
        }
    }
    out
}
