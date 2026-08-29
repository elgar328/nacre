//! **A cylinder's own chart, as an arrangement's line set** (capability D, first rung).
//!
//! The plane side long ago stopped walking faces by hand: each plane class gets a **cell complex**
//! (`arrangement`'s `walk_cells → nest_cells → label_cells → emit_faces`) and the case-work went
//! with it. The cylinder side did not. It still has three hand-written walks — `bands.rs`' band and
//! panel roads, and `boolean.rs`' `band_loop` slit-leg and `merge_curved_group` — and the last four
//! cells were each a falsehood inside one of them.
//!
//! `docs/design.md` names the way out: the lateral has an **isometric chart** `(z, r·θ)`, so the
//! same engine can run there — the seam is only the chart's cut line, and notches, holes and
//! θ-panels all become cells. This module is that road's **first** piece: the arrangement's two
//! axes, and nothing else.
//!
//! ## Why the lines are orthogonal — a derivation, not a measurement
//!
//! The M6-2a population gate (`planes.rs`) sorts every `(plane class, cylinder class)` pair into
//! five, and only two leave a mark on the lateral:
//!
//! | plane vs axis | distance | on the lateral | ☑ suite |
//! |---|---|---|---|
//! | `n ∥ m` | — | a **circle** — the chart's *horizontal* | 1128 |
//! | `n ⊥ m` | `d = 0` | **two rulings** — the chart's *vertical* | 75 |
//! | `n ⊥ m` | `0 < d ≤ r` | refused, `WallMeetsLateral` (capability B) | 7 |
//! | `n ⊥ m` | `d > r` | nothing | 1494 |
//! | oblique | — | refused, `ObliqueCylinderCut` (M6-3) | 2 |
//!
//! So within this milestone the chart carries **axis-parallel segments only**. Not a grid, though:
//! a ruling exists over a finite axis interval and a cut circle is arcs, so the cells are what
//! rectilinear *segments* cut out — which is why building them is its own rung.
//!
//! ★ The gate is the reason, so the day the gate changes this does too. Capability B opening
//! `0 < d ≤ r` adds vertical lines at **irrational** θ — and the chart never reads a θ *value*,
//! only an order, so nothing here would break. That is plausible and **unmeasured**; it is not
//! a reason this module leans on.
//!
//! ## What is not here
//!
//! The **cells**, the labels, and any cutover. The two roads this will one day replace are not
//! touched, because the census below has to measure what they answer *today*.

use crate::arrangement::{Curved, RulingExtent};
use crate::boolean::{Bound, LocalFace};
use crate::planes::{ClassIx, WorkingCyl, WorkingPlane};
use crate::tolerant::Judge;
use crate::{BoolError, combinatorics};
use nacre_scalar::Rat;

/// One **horizontal** line: a ⊥ plane class, and where it crosses the axis.
///
/// ★ `t` is exact (`Rat`) — a ⊥ class's axis parameter is a rational, and `bands::param` is the
/// one spelling of it. A class whose parameter is too wide for `Rat` refuses rather than being
/// skipped, exactly as `bands_of` already does: a missing band boundary merges two regions whose
/// membership differs, which is a closed and wrong answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ZLine {
    pub(crate) class: usize,
    pub(crate) t: Rat,
}

/// One **vertical** segment: a ruling, and the axis interval it spans.
///
/// Carried straight from [`RulingExtent`] — the plane pass decides this and the chart must not
/// decide it again.
///
/// ★ `#[allow(dead_code)]` for the same reason [`RulingExtent`] carries it: this rung *builds* the
/// line set and measures it; the cells that read `wall`/`side`/`end`/`z` are the next one.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub(crate) struct ThetaSeg {
    pub(crate) wall: usize,
    pub(crate) side: i8,
    /// The piece's own two branch nodes: the θ **order** is asked of these
    /// (`circular_order_about_seam` via `arrangement::circular_order`), never of a coordinate.
    pub(crate) end: [combinatorics::NodeId; 2],
    /// `[lower, upper]` on the axis.
    pub(crate) z: [Rat; 2],
}

/// One cylinder class's chart: the two axes, and nothing else yet.
#[derive(Clone, Debug, Default)]
pub(crate) struct Chart {
    pub(crate) z_lines: Vec<ZLine>,
    pub(crate) theta: Vec<ThetaSeg>,
}

/// **Every ⊥ class that cuts this cylinder, and every ruling on it** — per cylinder *class*.
///
/// ★★★★★ **The sources are `bands_of`'s, and one thing it does is deliberately *not* done here.**
/// That function ends with `keyed.retain(|(t, _)| t ∈ row.span)`: it clips to the **face** it was
/// asked about, because its row is a face and not a class. A class has no such span — one lateral
/// surface can carry several faces with gaps between them, which is exactly why `CylRow` became
/// per-face. So the clip is dropped, and that is a *repair* rather than a loss: `CylRow`'s own doc
/// warns that merging the spans into one `min..max` would **invent a band where the solid has no
/// face at all**, and on a chart that gap is simply a cell that keeps nothing — the argument the
/// plane side already makes for its unbounded cells.
///
/// ★ The cylinder class itself is one solid's surface, not both operands': `planes.rs` interns
/// cylinders by `Handle<Surface>` (planes by geometry), and two solids share no handles. That is
/// sound here only because a **coincident pair is refused** — `cylinders_clear` turns coaxial
/// cylinders away as `CylinderPairContact`. What does see both operands is what the chart is made
/// *of*: the plane classes, which are interned geometrically.
pub(crate) fn chart_of(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    k: usize,
    plane_faces: &[LocalFace],
    curved: &Curved,
) -> Result<Chart, BoolError> {
    let def = &cyls[k].def;
    let mut classes: Vec<usize> = Vec::new();
    let push = |c: usize, v: &mut Vec<usize>| {
        if !v.contains(&c) {
            v.push(c);
        }
    };
    // (i) the circles the plane arrangement actually emitted, and (ii) the cut rims, which emit
    // arcs instead and so are invisible to (i). `bands_of` reads both for the same reason.
    for lf in plane_faces {
        let ClassIx::Plane(c) = lf.surf else { continue };
        if std::iter::once(&lf.outer)
            .chain(lf.inner.iter())
            .any(|b| matches!(b, Bound::Circle { cyl } if *cyl == k))
        {
            push(c, &mut classes);
        }
    }
    for &(kk, c) in curved.cut_rims.keys() {
        if kk == k {
            push(c, &mut classes);
        }
    }
    // (iii) every ⊥ class, full stop. `bands_of` only takes the ones that land on a face's rim
    // because it is answering about a face; a chart covers the whole axis.
    for c in 0..jd.planes.len() {
        let Some(coeffs) = combinatorics::class_coeffs_rat(jd, c) else {
            continue;
        };
        if nacre_scalar::parallel_rat(&[coeffs[0], coeffs[1], coeffs[2]], &def.dir()) {
            push(c, &mut classes);
        }
    }
    let mut z_lines: Vec<ZLine> = classes
        .into_iter()
        .map(|class| crate::bands::axis_param(jd, class, def).map(|t| ZLine { class, t }))
        .collect::<Result<_, BoolError>>()?;
    z_lines.sort_by_key(|l| l.t);
    z_lines.dedup_by_key(|l| l.t);

    let theta: Vec<ThetaSeg> = curved
        .rulings
        .get(&k)
        .map(|v| {
            v.iter()
                .map(|r: &RulingExtent| {
                    let (lo, hi) = if r.z[0] <= r.z[1] {
                        (r.z[0], r.z[1])
                    } else {
                        (r.z[1], r.z[0])
                    };
                    ThetaSeg {
                        wall: r.wall,
                        side: r.side,
                        end: r.end,
                        z: [lo, hi],
                    }
                })
                .collect()
        })
        .unwrap_or_default();
    Ok(Chart { z_lines, theta })
}

/// **What the chart's line set looks like, beside what the hand-written roads answer today.**
///
/// The first rung of capability D ships no capability: it builds the two axes and *measures* them.
/// The numbers decide the next rung's design — above all whether the cells are worth walking or
/// worth overlaying — and until they are in, that design would rest on an unmeasured premise.
///
/// ★ Called from the one place the plane arrangement's faces and [`Curved`] are both in hand. It
/// reads; nothing downstream reads it.
pub(crate) fn census(
    jd: &Judge<'_, WorkingPlane>,
    cyls: &[WorkingCyl],
    plane_faces: &[LocalFace],
    curved: &Curved,
    rows: &[crate::bands::CylRow],
) {
    for (k, _) in cyls.iter().enumerate() {
        // ★ A class whose axis parameter cannot be stated is skipped, not recorded: `bands_of`
        // refuses on the very same `axis_param` a few lines further down, so this is not an
        // outcome the chart introduces. ☑ Measured **0** across the suite either way.
        let Ok(chart) = chart_of(jd, cyls, k, plane_faces, curved) else {
            continue;
        };
        // How today's road sees the same cylinder: one row per lateral *face*, each with its own
        // clipped span. The chart has one line set for the whole class.
        let mine: Vec<&crate::bands::CylRow> = rows.iter().filter(|r| r.class == k).collect();
        // ★★★★★ **The claim is asserted here, where the fact is made.** A test that reads the
        // ledger afterwards sees only what ran before it — this session measured that exact hole
        // (a census reported 764 of 2174 because it sat mid-suite). Asserting at the record makes
        // the coverage total and names the offending test in the panic.
        //
        // ★ The relation is **⊇, not equality**: a class's z-lines include every ⊥ plane, whether
        // or not a face reaches it, while a row's band boundaries are clipped to that face. So
        // what is checked is *cover*, and the surplus is what the chart gains — the gaps a row
        // cannot describe at all.
        let covers_rows = mine.iter().all(|r| {
            [r.span[0], r.span[1]]
                .iter()
                .all(|t| chart.z_lines.iter().any(|l| l.t == *t))
        });
        assert!(
            covers_rows,
            "a chart lost a band boundary its own class's row had: cyl {k}, \
             z-lines {:?}, row spans {:?}",
            chart.z_lines.iter().map(|l| l.class).collect::<Vec<_>>(),
            mine.iter().map(|r| r.span).collect::<Vec<_>>()
        );
        probe::push(probe::Row {
            z_lines: chart.z_lines.len(),
            theta: chart.theta.len(),
            rows: mine.len(),
        });
    }
}

/// The census's ledger — the same shape as `ruling_probe`/`panel_probe`: filled where the fact is
/// made, read by one test that reports it.
pub(crate) mod probe {
    use std::sync::Mutex;

    /// One chart's shape: how many horizontal lines, how many vertical segments, and how many
    /// **rows** the hand-written road split the same class into.
    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Row {
        pub(crate) z_lines: usize,
        pub(crate) theta: usize,
        pub(crate) rows: usize,
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    pub(crate) fn push(r: Row) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
}

#[cfg(test)]
mod tests {
    use super::probe::{ROWS, Row};
    use crate::BoolKind;
    use nacre_math::{Point3, Vector3};
    use nacre_topo::Model;

    /// **The chart census is live** — the claim itself is asserted where the fact is made
    /// (`census`, "a chart covers its own class's rows"), because a test that reads the ledger
    /// sees only the booleans that ran before it.
    ///
    /// ★ So this builds its **own** cylinder boolean rather than leaning on the suite's order:
    /// a test that passes only when something else ran first is a test that fails under a filter,
    /// and this session has already left one of those behind.
    #[test]
    fn the_chart_census_is_running() {
        let before = ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .len();
        let mut m = Model::new();
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let drill = m.add_cylinder(
            Point3::from_array([2.0, 2.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            0.5,
            4.0,
        );
        m.rebuild_adjacency();
        crate::boolean(&mut m, BoolKind::Cut, plate, drill).expect("a through bore");
        let rows = ROWS
            .lock()
            .expect("the probe's lock is never held across a panic")
            .clone();
        assert!(
            rows.len() > before,
            "a boolean with a cylinder recorded no chart"
        );
        // A through bore: the plate's two caps and the drill's own two rims are the ⊥ classes it
        // meets, and **no wall holds its axis** — so the chart is all horizontal, which is the
        // shape 198 of the suite's 263 charts have (☑ measured). One row, because one lateral face.
        let Row {
            z_lines,
            theta,
            rows: n,
        } = *rows.last().expect("the census recorded a chart");
        assert_eq!((theta, n), (0, 1), "a plain bore's chart");
        assert!(z_lines >= 2, "a band needs two boundaries, got {z_lines}");
    }
}
