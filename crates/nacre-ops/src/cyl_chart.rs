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
    // ★ Deduped by **parameter**, and that is exact rather than lucky: two classes ⊥ to this axis
    // at one `t` are the *same plane*, which is `bands_of`'s own words for why it matches span ends
    // against `t` without a tolerance. So one of the two indices survives arbitrarily — which is
    // why everything downstream compares `t`, never the class index.
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
        //
        // ★★★★★ **The comparison asks `bands_of` itself, not the row's two ends.** The first
        // spelling checked `r.span[0]`/`r.span[1]` and the message still claimed "a band
        // boundary" — a weak check wearing a strong claim, which is the shape this session kept
        // catching. A row's boundaries come from three sources (emitted circle carriers, cut rims,
        // span ends), so only that function can say what they are.
        //
        // ★★ **Compared by `t`, never by class index.** The chart dedups its z-lines by axis
        // parameter, which is sound because two ⊥ classes at one `t` *are* one plane (`bands_of`
        // says so in its own words) — but it means the surviving line may carry the other class's
        // index. Matching on the index would fail for a reason that is not a lost boundary.
        //
        // ★★★★ **A refusal here is skipped, never a panic.** `bands_of` runs for real a few lines
        // below and may honestly refuse (a class whose axis parameter is too wide for `Rat`);
        // this census sits *before* it, so panicking would turn a legitimate reject into a crash
        // — an instrument changing the answer it came to watch. ☑ Measured 0 either way.
        for r in &mine {
            let Ok(bands) = crate::bands::bands_of(r, plane_faces, jd, &curved.cut_rims) else {
                continue;
            };
            for (lo, hi) in bands {
                for c in [lo, hi] {
                    let Ok(t) = crate::bands::axis_param(jd, c, &r.def) else {
                        continue;
                    };
                    assert!(
                        chart.z_lines.iter().any(|l| l.t == t),
                        "a chart lost a band boundary its own class's row had: cyl {k}, \
                         class {c}, z-lines {:?}",
                        chart.z_lines.iter().map(|l| l.class).collect::<Vec<_>>()
                    );
                }
            }
        }
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
    use nacre_scalar::Rat;
    use nacre_topo::Model;

    /// **The chart census is live, and every chart it records has the shape a chart must have.**
    ///
    /// The real claim — "a chart covers its own class's rows" — is asserted in `census`, where the
    /// fact is made, because a test that reads the ledger sees only what ran before it.
    ///
    /// ★★★★★ **And it may not read its *own* row out of that ledger.** The first spelling did
    /// (`rows.last()`), which is only its own under `--test-threads=1`; the workspace gate runs
    /// **parallel**, and a filtered run caught it immediately — the last row belonged to a
    /// `bands::tests` fixture (θ = 6 where a bore has 0). So what is asserted here is **universal
    /// over every recorded chart**, which no interleaving can break:
    ///
    /// * `rows >= 1` — a cylinder class exists because a face made it (`cyl_rows` refuses
    ///   otherwise), so the hand-written road always has at least one row for it.
    /// * `z_lines >= 2` — a band needs two boundaries; a lateral face's own two rims are ⊥ classes
    ///   and both are in the set. ☑ Measured minimum across the suite: exactly 2.
    /// * `theta` is **even** — a plane holding the axis cuts the lateral in *two* rulings
    ///   (`QuadRoot::{Lo, Hi}`), and nothing else contributes a vertical line. ☑ Measured: 0, 4,
    ///   6, 8.
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
        for r in &rows {
            let Row {
                z_lines,
                theta,
                rows: n,
            } = *r;
            assert!(n >= 1, "a class with no row: {r:?}");
            assert!(z_lines >= 2, "a band needs two boundaries: {r:?}");
            assert!(
                theta % 2 == 0,
                "a wall cuts two rulings, not {theta}: {r:?}"
            );
        }
    }

    /// **The derivation, on an axis that is not `+z`** — both roads, on a genuinely tilted one.
    ///
    /// ★★★★★ **Why this fixture had to be built.** Everything above is derived from the gate's
    /// branch table and depends on no axis direction, but every number beside it was measured on
    /// a corpus that is effectively axis-aligned: the kit builder raises cylinders on world XY
    /// (`+z` only), the ops corpus is `+z`/`+x`/`+y`, and all three *tilted* cylinders in it are
    /// **reject** fixtures — so no tilted axis had ever reached the chart. On `+z` against an
    /// axis-aligned box a class's axis parameter is just a `z` coordinate and the division in
    /// `bands::axis_param` is trivial; tilted, it is a real rational one.
    ///
    /// ★★★★ **Both halves are stated on the exact road, and that is the whole trick.**
    /// A *Pythagorean frame* — `u = (0.6, 0.8, 0)`, `v = (−0.48, 0.36, 0.8)` — lifts to axes that
    /// are exactly orthonormal in the rationals, so `SketchPlane::from_axes` takes the
    /// world-rational path and every face of the prism raised on it states narrow coefficients
    /// (asserted below, small integers). Its normal `(0.64, −0.48, 0.6)` is the cylinder's axis,
    /// and the frame's own `u` is a **unit rational perpendicular** to it — exactly the `ref_dir`
    /// `Model::add_cylinder_exact` requires.
    ///
    /// ★★★★★ **`Model::add_cylinder` cannot state this, and that is not a kernel limit.** That
    /// entry normalizes an `f64` axis and lifts its cap points back out of `any_perpendicular`'s
    /// computed floats, so on a tilted axis the caps land outside the narrow window and the
    /// population gate honestly declines (`CylinderGateUndecided` — measured while building this).
    /// Its own doc calls it *the test entry*; the production road is `add_cylinder_exact`, and on
    /// that road the tilted case is not the irrational case.
    ///
    /// ☑ What this establishes: the chart's two axes are built, and D1a's universal claims hold,
    /// on an axis with no zero component — through-bore (circles only) **and** a wall holding the
    /// axis (rulings).
    #[test]
    fn the_chart_stands_on_a_tilted_axis() {
        // Two runs of one fixture: the bore centred inside the prism (⊥ classes only), then
        // centred **on** a wall, which puts that wall through the axis and cuts rulings.
        // The prism is 2 x 2 x 1 and the bore has r = 1/2, so the removed volume is `pi/4`
        // whole and half of that when the axis lies in a wall.
        let bore = std::f64::consts::FRAC_PI_4;
        for (name, base, want_rulings, want_vol) in [
            ("through bore", [-0.52, 1.64, 0.2], false, 4.0 - bore),
            (
                "half bore on a wall",
                [0.08, 2.44, 0.2],
                true,
                4.0 - bore / 2.0,
            ),
        ] {
            let before = ROWS
                .lock()
                .expect("the probe's lock is never held across a panic")
                .len();
            let mut m = Model::new();
            let plane = crate::SketchPlane::from_axes(
                Point3::from_array([0.0; 3]),
                Vector3::from_array([0.6, 0.8, 0.0]),
                Vector3::from_array([-0.48, 0.36, 0.8]),
            );
            let frame = match crate::apply(
                &mut m,
                &crate::Operation::DatumPlane {
                    def: crate::DatumDef::Stated(plane),
                },
            ) {
                Ok(crate::OpOutput::DatumPlane { frame, .. }) => frame,
                other => panic!("stating the tilted plane: {other:?}"),
            };
            let p = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
            let crate::OpOutput::Extrude { solid, faces, .. } = crate::apply(
                &mut m,
                &crate::Operation::Extrude {
                    frame,
                    profile: crate::Profile2d::polygon(vec![
                        p(0.0, 0.0),
                        p(2.0, 0.0),
                        p(2.0, 2.0),
                        p(0.0, 2.0),
                    ])
                    .expect("a square"),
                    dist: 1.0,
                },
            )
            .expect("the tilted prism") else {
                unreachable!()
            };
            // ★★★★ **The fixture qualifies its own population** (the S2 census lesson, which this
            // file's neighbours keep citing): every face must state narrow world coefficients, and
            // each must be either ⊥ or ∥ to the axis **exactly** — otherwise the gate would be
            // answering the oblique branch and the chart would never be reached.
            let axis = [0.64, -0.48, 0.6].map(|x| Rat::from_decimal(x).expect("a short decimal"));
            for &f in &faces {
                let n = m
                    .surface_name
                    .get(&m.faces.get(f).surface)
                    .and_then(|nm| nm.narrow())
                    .copied()
                    .expect("a tilted prism's face states itself");
                let n = [n[0], n[1], n[2]];
                let dot = nacre_scalar::dot_sign_rat(&n, &axis);
                assert!(
                    dot == nacre_scalar::Orient::Zero || nacre_scalar::parallel_rat(&n, &axis),
                    "a face neither ⊥ nor ∥ to the axis would be the oblique branch"
                );
            }
            let q = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            let (cyl, _) = m
                .add_cylinder_exact(
                    base.map(q),
                    axis,
                    // ★ The frame's own `u`: unit, and `u · axis = 0.384 − 0.384 = 0` exactly.
                    [q(0.6), q(0.8), q(0.0)],
                    q(0.5),
                    q(3.0),
                    None,
                )
                .expect("the exact road states a tilted cylinder");
            m.rebuild_adjacency();
            let out = crate::boolean(&mut m, BoolKind::Cut, solid, cyl)
                .unwrap_or_else(|e| panic!("{name}: {e:?}"));
            assert!(
                nacre_validate::validate(&m).is_empty(),
                "{name}: the tilted result must be sound"
            );
            // ★★★★★ **The oracle is on this model, derived from the inputs.** Green plus a clean
            // `validate` would say only that *some* solid was built; the volume says the tilted
            // boolean removed the right material, and it is computed from the prism's and bore's
            // own dimensions rather than from anything the kernel answered.
            let [s] = out[..] else {
                panic!("{name}: one solid")
            };
            let got = nacre_props::mass_props(&m, s).expect("mass props").volume;
            assert!(
                (got - want_vol).abs() < 1e-9,
                "{name}: volume {got} vs {want_vol}"
            );
            let rows = ROWS
                .lock()
                .expect("the probe's lock is never held across a panic")
                .clone();
            assert!(rows.len() > before, "{name}: no chart was recorded");
            // ★★★ Read as a **set difference**, not `last()` (a `rows.last()` spelling had to be
            // corrected in D1a — the gate runs this suite in parallel). ★ And it stays an
            // *existential* claim over a **shared** ledger, so it can be diluted by a concurrent
            // test but never carries the weight alone: what proves this fixture right is the
            // volume above, on this model.
            let want = if want_rulings { 6 } else { 0 };
            assert!(
                rows[before..]
                    .iter()
                    .any(|r| r.z_lines == 4 && r.theta == want),
                "{name}: expected a chart with {want} rulings, got {:?}",
                &rows[before..]
            );
        }
    }
}
