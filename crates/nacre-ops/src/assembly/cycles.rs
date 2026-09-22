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
    /// A chain rim meeting the seam at more than two contacts, or passing θ = 0 other than
    /// once — the slit walk is not written for it.
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
    cut_rims: &crate::draft::CutRims,
) -> Result<LocalFace, CurvedAbstain> {
    // A **contact** is a node the split put *on* the seam (`CutRim::seam_is_node`) or an arc
    // that wraps past it (the loop builder splits that one at the rim's seam vertex) — the two
    // spellings are exclusive by `wrapping_rim`'s own guard, so neither is counted twice. It is
    // the slit's question — *which boundary vertices lie on θ = 0* — and not the winding's: a
    // hole that touches the seam ruling from one side never passes θ = 0 (winding 0, no
    // crossing) yet has two contacts, and the slit must still be spliced through it, or it
    // would run on top of the hole's own ruling where no shell guard can see it.
    let contacts = |r: &Ring| -> usize {
        let n = r.nodes.len();
        (0..n)
            .filter(|&t| {
                let on_seam = cut_rims.iter().any(|(&(kk, _), cr)| {
                    kk == k && cr.seam_is_node && cr.nodes.first() == Some(&r.nodes[t])
                });
                on_seam
                    || matches!(
                        wrapping_rim(
                            ClassIx::Cyl(k),
                            r.nodes[t],
                            r.nodes[(t + 1) % n],
                            r.walls[t],
                            cut_rims,
                        ),
                        Ok(Some(_))
                    )
            })
            .count()
    };

    // 5. Every threaded cycle, classified by its **winding about the axis** — Σ [`seam_step`]
    //    over its arc steps, read off the split's order table and no coordinate. `+1` is a lower
    //    boundary walked forward, `−1` an upper one walked backward, `0` a hole; a simple cycle
    //    on the lateral cannot wind twice. So the region's two rims are the `+1` cycle and the
    //    `−1` cycle, each a whole circle or a **chain** — the whole-circle case is this rule
    //    with the chains left out.
    struct Cycle {
        ring: Ring,
        w: i32,
        contacts: usize,
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
                cut_rims,
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
        let contacts = contacts(&r);
        cycles.push(Cycle {
            ring: r,
            w,
            contacts,
            crossings,
        });
    }

    // 6. ★★ **Only shapes the outer walk can bridge.** A hole that meets the seam generator is
    //    spliced into the band's outer boundary (`band_loop`), and that splice is written for
    //    **two** contacts on **at most one** hole; a chain is walked from one contact, and the
    //    slit's argument holds for a chain with at most two contacts and exactly one passage of
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
                if c.contacts > 2 || c.crossings != 1 {
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
