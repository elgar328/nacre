use super::*;
/// Why a region's cycles could not be classified into a band with holes — every abstention by
/// name. These are the assembly's own limits on chains and holes, stated where the
/// cycles are made ([`classify_cycles`]); the emitter refuses them as `ArcBoundNotYet`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CurvedAbstain {
    /// A step's passage of θ = 0 could not be named ([`seam_step`] declined).
    SeamUnnamed,
    /// A cycle winding more than once about the axis — not a simple boundary.
    Winding,
    /// The wrapping cycles are not exactly one of each sense (`+1` and `−1`).
    WrappingNotOneEach,
    /// A chain rim meeting the seam at more than two contact vertices, or passing θ = 0 other
    /// than once — the slit walk is not written for it.
    ChainContacts,
    /// A hole meets the seam at other than zero or two contacts.
    HoleContacts,
    /// More than one hole meets the seam.
    Bridging,
}

/// **A region's cycles classified by their winding about the axis, and the face they bound** —
/// `Band { lo, hi }` with holes. The chart's region emitter hands it its walk's cycles.
/// `rims_lo` / `rims_hi` are the whole (uncut) rims by the walk's sense.
pub(crate) fn classify_cycles(
    k: usize,
    flip: bool,
    rings: Vec<Ring>,
    rims_lo: Vec<usize>,
    rims_hi: Vec<usize>,
    rims: &crate::draft::HeldRims,
) -> Result<LocalFace, CurvedAbstain> {
    // A **contact** is a rim node *on* the seam (`CutRim::seam_is_node`) or an arc
    // that wraps past it (the loop builder splits that one at the rim's seam vertex) — the two
    // spellings are exclusive by `wrapping_rim`'s own guard, so neither is counted twice. It is
    // the slit's question — *which boundary vertices lie on θ = 0* — and not the winding's: a
    // hole that touches the seam ruling from one side never passes θ = 0 (winding 0, no
    // crossing) yet has two contacts, and the slit must still be spliced through it, or it
    // would run on top of the hole's own ruling where no shell guard can see it.
    //
    // ★ **Counted the way the walk that reads it reads it — two counts.** A hole is spliced by
    // `band_loop` at two contact **positions** of its ring (`hits = [i, j]`), so a hole counts
    // positions. A chain is entered at one contact **vertex**, picked by its station
    // (`reconstruct`'s `walk_of`), so a chain counts vertices: a chain pinched at a contact —
    // the inward wedge with its apex on the seam, past the top cap, visits the top rim's seam node
    // twice — has two contact vertices and three positions, and reading positions there abstained
    // before the shell guard could name the result touching itself along the apex ruling.
    let contacts = |r: &Ring| -> (usize, usize) {
        let n = r.nodes.len();
        let (mut positions, mut wraps) = (0usize, 0usize);
        let mut seam_nodes: Vec<NodeId> = Vec::new();
        for t in 0..n {
            let on_seam = rims.iter().any(|(&(kk, _), cr)| {
                kk == k && cr.seam_is_node && cr.nodes.first() == Some(&r.nodes[t])
            });
            if on_seam {
                positions += 1;
                if !seam_nodes.contains(&r.nodes[t]) {
                    seam_nodes.push(r.nodes[t]);
                }
            } else if matches!(
                wrapping_rim(
                    ClassIx::Cyl(k),
                    r.nodes[t],
                    r.nodes[(t + 1) % n],
                    r.walls[t],
                    rims,
                ),
                Ok(Some(_))
            ) {
                positions += 1;
                wraps += 1;
            }
        }
        (positions, seam_nodes.len() + wraps)
    };

    // 5. Every threaded cycle, classified by its **winding about the axis** — Σ [`seam_step`]
    //    over its arc steps, read off the rim's order table and no coordinate. `+1` is a lower
    //    boundary walked forward, `−1` an upper one walked backward, `0` a hole; a simple cycle
    //    on the lateral cannot wind twice. So the region's two rims are the `+1` cycle and the
    //    `−1` cycle, each a whole circle or a **chain** — the whole-circle case is this rule
    //    with the chains left out.
    struct Cycle {
        ring: Ring,
        w: i32,
        /// Contact positions — what a hole's splice reads.
        contacts: usize,
        /// Contact vertices — what a chain's walk reads.
        contact_vertices: usize,
        crossings: usize,
    }
    let mut cycles: Vec<Cycle> = Vec::with_capacity(rings.len());
    for r in rings {
        let n = r.nodes.len();
        let (mut w, mut crossings) = (0i32, 0usize);
        for t in 0..n {
            match seam_step(
                ClassIx::Cyl(k),
                r.nodes[t],
                r.nodes[(t + 1) % n],
                r.walls[t],
                rims,
            ) {
                Ok(Some((_, _, sign))) => {
                    w += i32::from(sign);
                    crossings += 1;
                }
                Ok(None) => {}
                Err(_) => return Err(CurvedAbstain::SeamUnnamed),
            }
        }
        if w.abs() > 1 {
            return Err(CurvedAbstain::Winding);
        }
        let (contacts, contact_vertices) = contacts(&r);
        cycles.push(Cycle {
            ring: r,
            w,
            contacts,
            contact_vertices,
            crossings,
        });
    }

    // 6. ★★ **Only shapes the outer walk can bridge.** A hole that meets the seam generator is
    //    spliced into the band's outer boundary (`band_loop`), and that splice is written for
    //    **two** contacts on **at most one** hole; a chain is walked from one contact, and the
    //    slit's argument holds for a chain with at most two contact vertices and exactly one passage of
    //    θ = 0 (the population measured; a chain that also crosses elsewhere is not arranged,
    //    and no walk is written for it — the `band_loop` precedent). Anything else abstains
    //    here rather than reaching the honest reject there — this pass protects the pass, and
    //    the whole-result judgements keep their witnesses (the charter `merge_component`
    //    records).
    let mut lo: Vec<Rim> = rims_lo.into_iter().map(Rim::Circle).collect();
    let mut hi: Vec<Rim> = rims_hi.into_iter().map(Rim::Circle).collect();
    let mut holes: Vec<Ring> = Vec::new();
    let mut bridging = 0usize;
    for c in cycles {
        match c.w {
            0 => {
                match c.contacts {
                    0 => {}
                    2 => bridging += 1,
                    _ => return Err(CurvedAbstain::HoleContacts),
                }
                holes.push(c.ring);
            }
            _ => {
                if c.contact_vertices > 2 || c.crossings != 1 {
                    return Err(CurvedAbstain::ChainContacts);
                }
                if c.w == 1 {
                    lo.push(Rim::Chain(c.ring));
                } else {
                    hi.push(Rim::Chain(c.ring));
                }
            }
        }
    }
    if bridging > 1 {
        return Err(CurvedAbstain::Bridging);
    }
    let (Ok([lo]), Ok([hi])) = (<[Rim; 1]>::try_from(lo), <[Rim; 1]>::try_from(hi)) else {
        return Err(CurvedAbstain::WrappingNotOneEach);
    };

    Ok(LocalFace {
        surf: ClassIx::Cyl(k),
        outer: Bound::Band { lo, hi },
        inner: holes.into_iter().map(Bound::Ring).collect(),
        flip,
    })
}
