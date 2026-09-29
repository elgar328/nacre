//! Curved cells: lateral cycles, crossings, rulings, disk-in-disk and rim witnesses.

use super::*;

/// A result's lateral faces' boundary cycles by kind — (rims, chains, panels, holes) summed
/// over the laterals; an unnamed cycle panics (the census locks that none is). The tracer's
/// inputs are set up as `pinned_ends_ordered` does.
fn lateral_cycle_census(m: &Model, r: Handle<Solid>) -> [usize; 4] {
    let faces_tab = collect_planes(m, r).expect("the result's face table");
    let mut surf_ix: HashMap<Handle<Face>, usize> = HashMap::new();
    for (i, pi) in faces_tab.iter().enumerate() {
        if let Some(fh) = pi.face() {
            surf_ix.insert(fh, i);
        }
    }
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, plane_ix, cyl_surfs) = dense_planes(&faces_tab, &canon);
    let inc = combinatorics::edge_faces(m, r, &surf_ix).expect("edge incidence");
    let cyls: Vec<crate::planes::WorkingCyl> = cyl_surfs
        .iter()
        .map(|&surf| {
            let def = m.world_cylinder_def(surf).expect("a world cylinder");
            let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(surf) else {
                unreachable!()
            };
            crate::planes::WorkingCyl {
                surf,
                def,
                realized: *cache,
                owner: crate::planes::SolidSide::A,
                shared: Vec::new(),
            }
        })
        .collect();
    let jd = crate::planes::test_judge(&planes);
    let mut tally = [0usize; 4];
    for &fh in &m.shell(m.solid(r).outer).faces {
        let fp = surf_ix[&fh];
        if !matches!(plane_ix[fp], crate::planes::ClassIx::Cyl(_)) {
            continue;
        }
        let cycles = combinatorics::lateral_cycles(m, fh, fp, &inc, &jd, &plane_ix, &cyls)
            .unwrap_or_else(|e| panic!("a lateral's cycles could not be named: {e:?}"));
        for (kind, _) in cycles {
            tally[match kind {
                combinatorics::CycleKind::Rim => 0,
                combinatorics::CycleKind::Chain => 1,
                combinatorics::CycleKind::Panel => 2,
                combinatorics::CycleKind::Hole => 3,
            }] += 1;
        }
    }
    tally
}

/// **The re-operation census**: can a boolean's *result* be an operand again? For every
/// family × kind on a 4×4×2 plate with an r = 0.5 boss, the result is cut with a cube that lies
/// far outside it (`[20, 21]³`) — the geometry cannot change, so the second boolean exercises only
/// what an operand needs: every face named as rings, every lateral row stated. What is locked is
/// the **outcome by name** (`Ok` with the volume unchanged, or the refusal's name — a
/// `TraceDeclined` by its `kind`, the face handle being a witness rather than a lock).
///
/// The table in the body is the statement — one row per [`BOSS_FAMILIES`] entry, each kind's
/// outcome beside the boundary cycles its result's laterals carry — and the distribution asserted
/// at the end reads it as a whole, so a moved cell is seen beside its neighbours.
///
/// ★ The probe beside it counts the loops carrying an `OnSeam` joint at any vertex — outer loops
/// of plane faces and holes of every face (a lateral's outer loop is joined at the seam by
/// construction, so it is not counted). Its hole count is the zero that says a hole never carries
/// one: a hole touching the seam is spliced into the outer walk.
#[test]
fn reop_census_families_reoperate_or_decline_by_name() {
    #[derive(Debug, PartialEq, Eq, Clone, Copy)]
    enum Reop {
        Ok,
        Empty,
        First(RejectReason),
        Rejected(RejectReason),
        Declined(DeclineKind),
    }
    use Reop::*;
    // A lateral's boundary cycles per kind, summed over the result's lateral faces:
    // (rims, chains, panels, holes).
    type Cyc = [usize; 4];
    const RIMS: Cyc = [2, 0, 0, 0];
    const RIMS_HOLE: Cyc = [2, 0, 0, 1];
    const CHAIN: Cyc = [1, 1, 0, 0];
    const PANEL: Cyc = [0, 0, 1, 0];
    const NONE: Cyc = [0; 4];
    // (name, [Fuse, Cut, Common], the cycles each result's laterals carry) — one row per
    // `BOSS_FAMILIES` entry, in its order.
    type Family = (&'static str, [Reop; 3], [Cyc; 3]);
    let families: [Family; 17] = [
        ("through", [Ok, Ok, Ok], [[4, 0, 0, 0], RIMS, RIMS]),
        // A boss standing on the plate: the cut removes nothing and leaves no lateral at all.
        ("on top", [Ok, Ok, Empty], [RIMS, NONE, NONE]),
        ("flush", [Ok, Ok, Ok], [RIMS, RIMS, RIMS]),
        ("wall -y", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall +y", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall -x", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("wall +x", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("corner", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("corner-lo", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("offmid", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("half wall", [Ok, Ok, Ok], [CHAIN, PANEL, PANEL]),
        ("half wall, cap below", [Ok, Ok, Ok], [CHAIN, PANEL, PANEL]),
        ("half +x", [Ok, Ok, Ok], [CHAIN, PANEL, PANEL]),
        ("half +x, cap below", [Ok, Ok, Ok], [CHAIN, PANEL, PANEL]),
        // The offset wall: a notch hole in the band, a panel for the cut and the common
        // — the wall families' shape, with the chord off the diameter.
        ("offset-out", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("offset-in", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
        ("offset-irr", [Ok, Ok, Ok], [RIMS_HOLE, PANEL, PANEL]),
    ];
    let mut seam_joints = (0usize, 0usize);
    let mut table: Vec<String> = Vec::new();
    let mut mismatches = 0usize;
    let mut cycle_mismatches: Vec<String> = Vec::new();
    for ((name, want, want_cyc), &(fam, base, h)) in families.into_iter().zip(&BOSS_FAMILIES) {
        assert_eq!(name, fam, "the table's rows follow the corpus");
        for ((kind, want), want_cyc) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
            .into_iter()
            .zip(want)
            .zip(want_cyc)
        {
            let (mut m, plate, boss) = boss_family(base, h);
            let got = match boolean(&mut m, kind, plate, boss) {
                Err(BoolError::Rejected { reason, .. }) => First(reason),
                Err(e) => panic!("{name} {kind:?}: first op {e:?}"),
                Result::Ok(out) if out.is_empty() => Empty,
                Result::Ok(out) => {
                    m.rebuild_adjacency();
                    // The seam-joint probe, over the result's faces.
                    let solid = m.solid(out[0]);
                    for sh in std::iter::once(solid.outer).chain(solid.cavities.iter().copied()) {
                        for &fh in &m.shell(sh).faces {
                            let face = m.face(fh);
                            let lateral = matches!(
                                m.surface(face.surface),
                                nacre_topo::Surface::Cylinder { .. }
                            );
                            let joint = |lp: &nacre_topo::Loop| {
                                lp.half_edges.len() >= 2
                                    && lp.half_edges.iter().any(|&he| {
                                        matches!(
                                            *m.vertex(m.he_start(he)),
                                            nacre_topo::Vertex::OnSeam(_)
                                        )
                                    })
                            };
                            if !lateral && joint(&face.outer) {
                                seam_joints.0 += 1;
                            }
                            seam_joints.1 += face.inner.iter().filter(|lp| joint(lp)).count();
                        }
                    }
                    // The lateral faces' boundary cycles, named as the tracer will read
                    // them — the outer loop cut at its slits, the spliced hole recovered.
                    let got_cyc = lateral_cycle_census(&m, out[0]);
                    if got_cyc != want_cyc {
                        cycle_mismatches
                            .push(format!("{name} {kind:?}: {got_cyc:?} ← want {want_cyc:?}"));
                    }
                    let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                    let far =
                        m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
                    m.rebuild_adjacency();
                    match boolean(&mut m, BoolKind::Cut, out[0], far) {
                        Result::Ok(r) => {
                            assert_eq!(r.len(), 1, "{name} {kind:?}: the far cut keeps one solid");
                            m.rebuild_adjacency();
                            let issues = nacre_validate::validate(&m);
                            assert!(issues.is_empty(), "{name} {kind:?}: {issues:?}");
                            let v = nacre_props::mass_props(&m, r[0]).expect("props").volume;
                            assert!(
                                (v - v0).abs() < 1e-9,
                                "{name} {kind:?}: the far cut changed the volume {v0} → {v}"
                            );
                            Ok
                        }
                        Err(BoolError::Rejected {
                            reason: RejectReason::TraceDeclined { kind, .. },
                            ..
                        }) => Declined(kind),
                        Err(BoolError::Rejected { reason, .. }) => Rejected(reason),
                        Err(e) => panic!("{name} {kind:?}: {e:?}"),
                    }
                }
            };
            table.push(format!(
                "{name} {kind:?}: {got:?}{}",
                if got == want { "" } else { "  ← want " }
            ));
            if got != want {
                let last = table.len() - 1;
                table[last].push_str(&format!("{want:?}"));
                mismatches += 1;
            }
        }
    }
    // The whole table at once, so a moved cell is read beside its neighbours.
    assert_eq!(mismatches, 0, "re-operation table:\n{}", table.join("\n"));
    assert!(
        cycle_mismatches.is_empty(),
        "lateral cycles:\n{}",
        cycle_mismatches.join("\n")
    );
    // The distribution, so a drift in the table is read as a whole.
    let count = |p: fn(&Reop) -> bool| -> usize {
        families
            .iter()
            .flat_map(|(_, w, _)| w.iter())
            .filter(|r| p(r))
            .count()
    };
    // Every kind of every family re-operates but one: `on top`'s Common, which is empty.
    assert_eq!(count(|r| *r == Ok), 50, "{table:?}");
    // ★ No row waits on the chart's reader. A new population must restate this zero.
    assert_eq!(
        count(|r| *r == Rejected(RejectReason::CylinderGateUndecided)),
        0,
        "no row waits on the chart's reader"
    );
    assert_eq!(
        count(|r| *r == Rejected(RejectReason::DegenerateFace)),
        0,
        "a circle-bounded face is never degenerate"
    );
    assert_eq!(count(|r| *r == Declined(DeclineKind::CylSpan)), 0);
    assert_eq!(
        count(|r| *r == Declined(DeclineKind::OuterRing)),
        0,
        "every ring of a result face has a name"
    );
    // Six of the outer loops are the offset families' Fuse rows — two whole-circle rims each
    // that pass the seam vertex. The wall families' seam sits on the wall itself (a chord end, a
    // pierce vertex), which is why their Fuse rows add none of these.
    assert_eq!(seam_joints, (14, 0), "seam-joint loops (outer, holes)");
}

/// **A radius-`0.5` disk centred on `base`'s `xy`, clipped by an axis-aligned rectangle — the
/// closed form.** The chord length integrated over the rectangle's `x` span:
/// `∫ [min(c₁, s) + min(c₀, s)] dt` with `s(t) = √(r² − t²)`, `c₁ = y₁ − cy`, `c₀ = cy − y₀`, over
/// `t ∈ [x₀ − cx, x₁ − cx]` clipped to where the chord meets `[y₀, y₁]` (`|t| ≤ √(r² − d²)`, `d`
/// the centre's distance to that span, 0 when the centre is inside it — there `min(c, s)` is
/// never negative). `∫ min(c, s) dt` splits at `±√(r² − c²)`: `c` between, `s` outside, with
/// `∫ s = (t·s + r²·asin(t/r))/2`. Exact as a closed form in `π`, `asin`, `sqrt`, and it reads
/// every rectangle — the edge through the centre (a halving),
/// the edge off it (a segment), two edges, a rectangle inside the disk.
///
/// ★ One spelling, read by [`removed_by`] and [`first_volume`] alike — they ask the same question
/// about the same disk, and two copies would be free to drift. `disk_in_rect_is_the_closed_form`
/// locks it against the known values and a numeric integration.
fn disk_in_rect(base: [f64; 3], rect: [[f64; 2]; 2]) -> f64 {
    let r = 0.5;
    let (cx, cy) = (base[0], base[1]);
    let (x0, x1) = (rect[0][0], rect[0][1]);
    let (y0, y1) = (rect[1][0], rect[1][1]);
    if x1 <= x0 || y1 <= y0 {
        return 0.0;
    }
    // The chord at `t` meets the y span only where `|t| ≤ t_max`.
    let d_out = (y0 - cy).max(cy - y1).max(0.0);
    if d_out >= r {
        return 0.0;
    }
    let t_max = (r * r - d_out * d_out).sqrt();
    let a = (x0 - cx).max(-t_max);
    let b = (x1 - cx).min(t_max);
    if b <= a {
        return 0.0;
    }
    // ∫ s(t) dt, the antiderivative.
    let big_f = |t: f64| {
        let t = t.clamp(-r, r);
        (t * (r * r - t * t).max(0.0).sqrt() + r * r * (t / r).asin()) / 2.0
    };
    // ∫_a^b min(c, s(t)) dt.
    let int_min = |c: f64| -> f64 {
        if c >= r {
            return big_f(b) - big_f(a);
        }
        if c <= 0.0 {
            // The edge is on the far side of the centre: within the domain `s ≥ |c|`, so the
            // minimum is the (non-positive) constant.
            return c * (b - a);
        }
        let tc = (r * r - c * c).sqrt();
        let (lo, hi) = (a.max(-tc), b.min(tc));
        let flat = if hi > lo { c * (hi - lo) } else { 0.0 };
        let left = if a < -tc {
            big_f((-tc).min(b)) - big_f(a)
        } else {
            0.0
        };
        let right = if b > tc {
            big_f(b) - big_f(tc.max(a))
        } else {
            0.0
        };
        flat + left + right
    };
    int_min(y1 - cy) + int_min(cy - y0)
}

/// [`disk_in_rect`] against the values the corpus has always used and the offset wall's
/// measured segments, then against a midpoint integration on rectangles of every shape (an oracle
/// that shares nothing with the closed form).
#[test]
fn disk_in_rect_is_the_closed_form() {
    use std::f64::consts::PI;
    let (r, c) = (0.5, [1.0, 1.0, 0.0]);
    let seg = |d: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let near = |a: f64, b: f64| (a - b).abs() < 1e-12;
    let full = PI * r * r;
    assert!(
        near(disk_in_rect(c, [[0.0, 2.0], [0.0, 2.0]]), full),
        "the whole disk"
    );
    assert!(
        near(disk_in_rect(c, [[1.0, 2.0], [0.0, 2.0]]), full / 2.0),
        "a halving"
    );
    assert!(
        near(disk_in_rect(c, [[1.0, 2.0], [1.0, 2.0]]), full / 4.0),
        "two halvings"
    );
    assert!(
        near(disk_in_rect(c, [[1.3, 2.0], [0.0, 2.0]]), seg(0.3)),
        "the segment beyond x = 1.3"
    );
    assert!(
        near(disk_in_rect(c, [[0.0, 1.3], [0.0, 2.0]]), full - seg(0.3)),
        "the disk less it"
    );
    assert!(
        near(disk_in_rect(c, [[1.3, 2.0], [1.0, 2.0]]), seg(0.3) / 2.0),
        "half a segment"
    );
    assert!(
        near(disk_in_rect(c, [[0.0, 0.7], [0.0, 2.0]]), seg(0.3)),
        "the far side's segment"
    );
    assert!(
        near(disk_in_rect(c, [[0.9, 1.1], [0.9, 1.1]]), 0.04),
        "a rectangle inside the disk"
    );
    assert!(
        near(disk_in_rect(c, [[1.6, 2.0], [0.0, 2.0]]), 0.0),
        "a rectangle beyond it"
    );
    assert!(
        near(disk_in_rect(c, [[0.0, 2.0], [1.0, 1.0]]), 0.0),
        "an empty span"
    );
    // The numeric oracle: the chord length sampled at midpoints.
    let numeric = |rect: [[f64; 2]; 2]| -> f64 {
        let n = 400_000;
        let (x0, x1) = (rect[0][0].max(c[0] - r), rect[0][1].min(c[0] + r));
        if x1 <= x0 {
            return 0.0;
        }
        let h = (x1 - x0) / n as f64;
        (0..n)
            .map(|i| {
                let t = x0 + (i as f64 + 0.5) * h - c[0];
                let s = (r * r - t * t).max(0.0).sqrt();
                let (lo, hi) = (rect[1][0].max(c[1] - s), rect[1][1].min(c[1] + s));
                (hi - lo).max(0.0) * h
            })
            .sum()
    };
    for rect in [
        [[0.0, 2.0], [0.8, 1.2]],
        [[0.7, 1.3], [0.6, 1.1]],
        [[1.1, 1.45], [0.55, 0.9]],
        [[0.55, 1.2], [1.15, 2.0]],
        [[1.2, 2.0], [1.3, 2.0]],
        [[0.6, 1.4], [0.6, 1.4]],
    ] {
        let (a, b) = (disk_in_rect(c, rect), numeric(rect));
        assert!(
            (a - b).abs() < 1e-7,
            "{rect:?}: closed form {a} vs numeric {b}"
        );
    }
}

/// **The volume of a family's first result**, from the same closed form — the plate `4×4×2`, the
/// boss `πr²h`, and their overlap (the disk clipped by the plate's footprint, times the axial
/// overlap): `Fuse = plate + boss − both`, `Cut = plate − both`, `Common = both`.
///
/// ★★★★★ **`removed_by` never checked this.** It locks what a *tool* removes from the first
/// result, taking that result's own volume as the baseline — so a first operation could be wrong
/// by a whole feature and every row would still agree with itself. The day `corner-lo` opened is
/// the day that mattered: a family that had only ever been refused arrived with no independent
/// number of its own. ☑ Its three kinds land on 34.748893572 / 31.607300918 / 0.392699082, to
/// 7e-15 or exactly.
fn first_volume(kind: BoolKind, base: [f64; 3], h: f64) -> f64 {
    let seg = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] { [a[0].max(b[0]), a[1].min(b[1])] };
    let len = |a: [f64; 2]| (a[1] - a[0]).max(0.0);
    let r = 0.5;
    let plate = 4.0 * 4.0 * 2.0;
    let boss = std::f64::consts::PI * r * r * h;
    let both =
        disk_in_rect(base, [[0.0, 4.0], [0.0, 4.0]]) * len(seg([base[2], base[2] + h], [0.0, 2.0]));
    match kind {
        BoolKind::Fuse => plate + boss - both,
        BoolKind::Cut => plate - both,
        BoolKind::Common => both,
    }
}

/// The volume a tool box removes from a family's **first** result, summed from parts: the plate's
/// box ∩ tool, the boss's disk ∩ the tool's footprint times the axial overlap, and the doubly
/// counted disk ∩ plate ∩ tool — `Fuse = plate + boss − both`, `Cut = plate − both`,
/// `Common = both`. A disk clipped by a rectangle is [`disk_in_rect`]'s closed form — `πr²/2ᵏ`
/// for `k` edges through the centre, a segment for an edge off it, whatever the corpus
/// and panics rather than approximating.
fn removed_by(kind: BoolKind, base: [f64; 3], h: f64, tool: [[f64; 3]; 2]) -> f64 {
    let seg = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] { [a[0].max(b[0]), a[1].min(b[1])] };
    let len = |a: [f64; 2]| (a[1] - a[0]).max(0.0);
    let plate = [[0.0, 0.0, 0.0], [4.0, 4.0, 2.0]];
    let axis = |b: [[f64; 3]; 2], i: usize| [b[0][i], b[1][i]];
    let box_vol = |a: [[f64; 3]; 2], b: [[f64; 3]; 2]| -> f64 {
        (0..3).map(|i| len(seg(axis(a, i), axis(b, i)))).product()
    };
    let disk_in_rect = |rect: [[f64; 2]; 2]| -> f64 { disk_in_rect(base, rect) };
    let tool_xy = [axis(tool, 0), axis(tool, 1)];
    let plate_xy = [seg(tool_xy[0], [0.0, 4.0]), seg(tool_xy[1], [0.0, 4.0])];
    let boss_z = [base[2], base[2] + h];
    let l_bt = len(seg(boss_z, axis(tool, 2)));
    let l_pbt = len(seg(seg(boss_z, [0.0, 2.0]), axis(tool, 2)));
    let plate_tool = box_vol(plate, tool);
    let boss_tool = l_bt * disk_in_rect(tool_xy);
    let both = l_pbt * disk_in_rect(plate_xy);
    match kind {
        BoolKind::Fuse => plate_tool + boss_tool - both,
        BoolKind::Cut => plate_tool - both,
        BoolKind::Common => both,
    }
}

/// **The crossing census**: a boolean's result cut by a tool whose faces actually **cross**
/// the result's rings — where the far cube of the re-operation census could see only whether the
/// rows are stated. Five tools per family × kind: a slab through the plate's middle
/// (`z ∈ [0.5, 1.5]`, its ⊥ caps crossing the wall faces that carry a boss's rulings, and parting
/// the result in two), a slab above it (`[2.5, 3.5]`, crossing the bands standing on the plate), a
/// slab below (`[−1.5, −0.5]`), a box whose wall passes **through the boss's axis** (so it crosses
/// the caps' arcs and cuts the lateral along two rulings), and a **wall slab** — the mid slab's
/// height, past the line one unit inside the wall the boss stands on, so its caps cross that
/// wall's rulings while the plate stays one solid. What is locked is the outcome by name —
/// `Ok(n)` with the volume equal to the first result's minus [`removed_by`] — and the counts.
///
/// The `want` table in the body is the statement, one row per [`BOSS_FAMILIES`] entry; the
/// constants name the shapes several families share, and the distribution asserted at the end
/// reads the table as a whole.
#[test]
fn crossing_census_slabs_and_through_axis_walls_by_name() {
    #[derive(Debug, PartialEq, Clone, Copy)]
    enum Cross {
        Ok(usize),
        First,
        Empty,
        Rejected(RejectReason),
        Declined(DeclineKind),
    }
    use Cross::*;
    let mid = [[-1.0, -1.0, 0.5], [5.0, 5.0, 1.5]];
    let top = [[-1.0, -1.0, 2.5], [5.0, 5.0, 3.5]];
    let bottom = [[-1.0, -1.0, -1.5], [5.0, 5.0, -0.5]];
    // The through-axis wall: `y = by` for a boss on an x-wall or in the interior, `x = bx` for one
    // on a y-wall — the wall that is not the plate's own.
    let axis_wall = |b: [f64; 3]| -> [[f64; 3]; 2] {
        if b[1] == 0.0 || b[1] == 4.0 {
            [[b[0], b[1] - 3.0, -2.0], [b[0] + 3.0, b[1] + 3.0, 6.0]]
        } else {
            [[b[0] - 3.0, b[1], -2.0], [b[0] + 3.0, b[1] + 3.0, 6.0]]
        }
    };
    // The wall slab: the mid slab's height, but only past the line one unit inside the wall the
    // boss stands on — its cap classes cross that wall face on the boss's rulings while the plate
    // stays one solid (one question, apart from the mid slab's second one: splitting the
    // result). For an interior boss it clears the boss and is a control.
    let wall_slab = |b: [f64; 3]| -> [[f64; 3]; 2] {
        if b[0] == 4.0 {
            [[3.0, -1.0, 0.5], [6.0, 5.0, 1.5]]
        } else if b[0] == 0.0 {
            [[-2.0, -1.0, 0.5], [1.0, 5.0, 1.5]]
        } else if b[1] == 4.0 {
            [[-1.0, 3.0, 0.5], [5.0, 6.0, 1.5]]
        } else if b[1] == 0.0 {
            [[-1.0, -2.0, 0.5], [5.0, 1.0, 1.5]]
        } else {
            [[3.0, -1.0, 0.5], [6.0, 5.0, 1.5]]
        }
    };
    const TOOLS: [&str; 5] = ["mid", "top", "bottom", "axis wall", "wall slab"];
    // Per family, per kind: the five tools' outcomes, in `TOOLS` order.
    use RejectReason::{NoClearRay, PierceVertexUnnamed, RingHasNoWitness, RulingBoundNotYet};
    // A band result: the mid slab parts the plate into two solids, the top and bottom slabs cut
    // the standing boss, the wall slab clears the boss. The through-axis wall builds: the chord's
    // cell reads its cylinder from the class table, and a wall **panel** ring hands over the
    // midpoint of an edge whose ends are one solve's two roots (`combinatorics::conjugate_midpoint`).
    // The Common column reads the same.
    const BAND: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    // A wall boss fused: the wall slab's caps cross the plate's wall face on the boss's
    // **rulings** — named as pierce nodes — and the cut builds with its exact volume.
    // The mid slab does the same and then splits the result in two, and the label decides:
    // every deciding probe's rays run clear of the boss, and a miss counts 0 against a
    // panel/chain lateral, so both halves classify and build with their exact volumes. The
    // through-axis wall crosses the bitten cap's **arc** and builds.
    const WALL_FUSE: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    // A wall boss cut (its lateral a panel — the notch's wall): the top and bottom slabs miss
    // the result (nothing of the boss stands outside the plate), the wall slab crosses the
    // panel's rulings and builds with its exact volume, the mid slab splits the result and both
    // halves classify by rays that miss the boss (as `WALL_FUSE`), and the through-axis wall
    // crosses the bite's arcs and reads each sector as the **run** of the rim's arcs it spans,
    // so it builds with its exact volume too.
    const PANEL_CUT: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    // A half boss fused (its lateral a chain, its top cap a half-disk): the cap's
    // chord is a graze on the wall class and the result builds — the top and bottom slabs miss
    // the boss, the wall slab crosses the chain's rulings with its exact volume, the mid slab
    // splits the result and both halves classify by rays that miss the boss (as `WALL_FUSE`),
    // the through-axis wall crosses the cap's arc. The **cap below** variants, whose seated
    // half-disk is the boss's *bottom* cap, read the same.
    const HALF_FUSE: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    // A half boss's Common builds under every tool: a ring bounded by one circle names that
    // circle's centre (`nesting::ring_own_circle`), and a wall panel's chord names its midpoint.
    const HALF_COMMON: [Cross; 5] = [Ok(1); 5];
    // A wall boss's Common — the inner half-cylinder alone, two half-disk caps and a flat side.
    // A **slab** parts it into two components whose faces are half-disk caps, a panel and a
    // lateral — not one vertex among them — and the depth classification takes its witness from
    // a cut cap read as the disk it is (`combinatorics::face_circle`, offered by
    // `combinatorics::coord_probes`). The **through-axis wall** cuts the other way and leaves one
    // solid, each sector read as the run of the rim's arcs it spans. ☑ In every family that
    // takes this constant that tool removes real volume, so each `Ok(1)` is a real cut.
    const COMMON: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(2)];
    // A half boss's Cut (a notch with a half-disk ceiling): the arc label carries its own side,
    // so the ceiling's two ends agree and the notch re-cuts — the mid slab parts it in two, the
    // others take one solid, all with their exact volumes. The through-axis wall's sector is not
    // an adjacent node pair of the rim; read as the **run** between the nearest rim nodes, it
    // answers as well.
    const HALF_CUT: [Cross; 5] = [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)];
    let want: [[[Cross; 5]; 3]; 17] = [
        [BAND, BAND, BAND], // through
        // on top: the cut leaves the plate alone.
        [BAND, [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)], [Empty; 5]],
        [BAND, BAND, BAND],             // flush
        [WALL_FUSE, PANEL_CUT, COMMON], // wall -y
        [WALL_FUSE, PANEL_CUT, COMMON], // wall +y
        [WALL_FUSE, PANEL_CUT, COMMON], // wall -x
        [WALL_FUSE, PANEL_CUT, COMMON], // wall +x
        // corner: the through-axis tool's wall `x = 4` is the plate's own wall (coplanar
        // contact), and the box lies wholly outside the plate beyond it, so `removed_by` is 0
        // for both kinds and what the volume oracle locks there is that the result is
        // **unchanged**. The **mid** slab parts the Fuse and the Cut in two, and the severed
        // halves' depth classification casts its rays from the plate's corner — every ray meets a
        // ring corner there, and a corner on the ray is a decision, not a tie. The Common is a
        // quarter cylinder alone, and either slab (mid, or the wall slab — both z ∈ [½, 1½])
        // parts it in two; its halves' rays from the axis vertex cross the *other* half's
        // **lateral**, a panel whose boundary loops answer the crossing question on the
        // cylinder's own chart: two solids, decided on the first attempt. The through-axis wall
        // cuts the other way and leaves one.
        [
            [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)],
            [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)],
            [Ok(2), Ok(1), Ok(1), Ok(1), Ok(2)],
        ],
        // corner-lo: its axis sits on **both** plate walls, so the chords through the circle are
        // radii and the sector outside the plate is reflex at the centre — the ring's true
        // extremum is in an arc, not at a node, and the winding is read there. The tools land on
        // the roads the other corner family sits on: the mid slab's Fuse and Cut decide at the
        // corner, the Common's halves at the other half's lateral.
        [
            [Ok(2), Ok(1), Ok(1), Ok(2), Ok(1)],
            [Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)],
            [Ok(2), Ok(1), Ok(1), Ok(0), Ok(2)],
        ],
        // offmid: the top slab's cap at z = 2.5 meets the notch's rulings above the plate — the
        // severed top piece has no vertex, its coordinate probe forks the other component's
        // mixed ring to the rational walk, and rays that miss the boss decide (the one cell
        // outside the mid column that differs from `WALL_FUSE`).
        [[Ok(2), Ok(2), Ok(1), Ok(1), Ok(1)], PANEL_CUT, COMMON],
        [HALF_FUSE, HALF_CUT, HALF_COMMON], // half wall
        [HALF_FUSE, HALF_CUT, HALF_COMMON], // half wall, cap below
        [HALF_FUSE, HALF_CUT, HALF_COMMON], // half +x
        [HALF_FUSE, HALF_CUT, HALF_COMMON], // half +x, cap below
        // ★ The offset wall. Every tool builds as it does for the wall families: the
        // mid slab parts the result, the through-axis wall (`y ≥ 2`, ⊥ to the plate's wall)
        // halves the segment. `offset-out`'s Common is a 0.2-deep segment prism: its halves'
        // corners are pierce names (the vertex probe drops them) and the cut cap's centre lies
        // outside a segment thinner than `r/2`, so the witness is one of the two points the cap
        // offers per chord (`combinatorics::ring_interior_candidates`): two solids.
        [[Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)], PANEL_CUT, COMMON],
        [[Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)], PANEL_CUT, COMMON], // offset-in
        [[Ok(2), Ok(1), Ok(1), Ok(1), Ok(1)], PANEL_CUT, COMMON], // offset-irr
    ];
    let mut table: Vec<String> = Vec::new();
    let mut mismatches = 0usize;
    let mut tally: Vec<(Cross, usize)> = Vec::new();
    // Owned, so the tie and deciding rows a refused case prints below are this corpus's own.
    crate::ledger::owned(|| {
        for (&(name, base, h), want) in BOSS_FAMILIES.iter().zip(want) {
            for (kind, want) in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common]
                .into_iter()
                .zip(want)
            {
                let tools = [mid, top, bottom, axis_wall(base), wall_slab(base)];
                for ((tool_name, tool), want) in TOOLS.into_iter().zip(tools).zip(want) {
                    // The mixed road's abstentions, attributed to this cell.
                    let tie0 = crate::combinatorics::tie_probe::len();
                    let dec0 = crate::assembly::probe::deciding::ROWS.len();
                    let (mut m, plate, boss) = boss_family(base, h);
                    let got = match boolean(&mut m, kind, plate, boss) {
                        Err(BoolError::Rejected { .. }) => First,
                        Err(e) => panic!("{name} {kind:?}: first op {e:?}"),
                        Result::Ok(out) if out.is_empty() => Empty,
                        Result::Ok(out) => {
                            m.rebuild_adjacency();
                            let v0 = nacre_props::mass_props(&m, out[0]).expect("props").volume;
                            // ★ The first result's own volume against the closed form, so the
                            // baseline `removed_by` subtracts from is itself checked. Summed over the
                            // whole result rather than guarded on one solid: a guard is a branch that
                            // can go quietly untaken, and the sum is right however many pieces a first
                            // operation leaves.
                            let v_first: f64 = out
                                .iter()
                                .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
                                .sum();
                            let w = first_volume(kind, base, h);
                            assert!(
                                (v_first - w).abs() < 1e-9,
                                "{name} {kind:?}: first volume {v_first} ≠ {w}"
                            );
                            let t = m.add_cuboid(
                                Point3::from_array(tool[0]),
                                Point3::from_array(tool[1]),
                            );
                            m.rebuild_adjacency();
                            match boolean(&mut m, BoolKind::Cut, out[0], t) {
                                Result::Ok(r) => {
                                    m.rebuild_adjacency();
                                    let issues = nacre_validate::validate(&m);
                                    assert!(
                                        issues.is_empty(),
                                        "{name} {kind:?} × {tool_name}: {issues:?}"
                                    );
                                    let v: f64 = r
                                        .iter()
                                        .map(|&s| {
                                            nacre_props::mass_props(&m, s).expect("props").volume
                                        })
                                        .sum();
                                    let expect = v0 - removed_by(kind, base, h, tool);
                                    assert!(
                                        (v - expect).abs() < 1e-9,
                                        "{name} {kind:?} × {tool_name}: volume {v} ≠ {v0} − removed = {expect}"
                                    );
                                    Ok(r.len())
                                }
                                Err(BoolError::Rejected {
                                    reason: RejectReason::TraceDeclined { kind, .. },
                                    ..
                                }) => Declined(kind),
                                Err(BoolError::Rejected { reason, .. }) => Rejected(reason),
                                Err(e) => panic!("{name} {kind:?} × {tool_name}: {e:?}"),
                            }
                        }
                    };
                    if got == Rejected(NoClearRay) {
                        let hist = crate::combinatorics::tie_probe::since(tie0);
                        eprintln!("tie {name} {kind:?} x {tool_name}: {hist:?}");
                        let dec = crate::assembly::probe::deciding::ROWS.mine_since(dec0);
                        for r in dec.iter().filter(|r| !r.2) {
                            eprintln!(
                                "deciding {name} {kind:?} x {tool_name}: exhausted offered {} ties {:?}",
                                r.1, r.3
                            );
                        }
                        let decided: Vec<usize> = dec.iter().filter(|r| r.2).map(|r| r.0).collect();
                        eprintln!(
                            "deciding {name} {kind:?} x {tool_name}: decided tried {decided:?}"
                        );
                    }
                    match tally.iter_mut().find(|(c, _)| *c == got) {
                        Some((_, n)) => *n += 1,
                        None => tally.push((got, 1)),
                    }
                    let ok = want == got;
                    table.push(format!(
                        "{name} {kind:?} × {tool_name}: {got:?}{}",
                        if ok {
                            String::new()
                        } else {
                            format!("  ← want {want:?}")
                        }
                    ));
                    if !ok {
                        mismatches += 1;
                    }
                }
            }
        }
    });
    assert_eq!(
        mismatches,
        0,
        "crossing table:\n{}\n\ntally: {tally:?}",
        table.join("\n")
    );
    // The distribution, so a moved cell is read as a whole.
    let count = |p: fn(&Cross) -> bool| -> usize {
        want.iter().flatten().flatten().filter(|c| p(c)).count()
    };
    assert_eq!(count(|c| matches!(c, Ok(_))), 250, "{tally:?}");
    // ★ Each refusal name below is empty in this corpus, and each zero is an emptiness of *this*
    // corpus, not of the road that names it — a new population must restate it.
    //
    // `NoClearRay`: every probe of the 3D depth classification abstained. Rays from a corner on
    // the ray decide there, and a ray crossing a panel or chain lateral is answered by the
    // lateral's own boundary loops.
    assert_eq!(count(|c| *c == Rejected(NoClearRay)), 0);
    // `RingHasNoWitness`: a ring with no witness at all (not «every witness blocked») — a thin
    // segment's cut cap offers two points per chord.
    assert_eq!(count(|c| *c == Rejected(RingHasNoWitness)), 0);
    // `PierceVertexUnnamed`: the names road builds carriers from its walls, so nothing dies at
    // ring construction.
    assert_eq!(count(|c| *c == Rejected(PierceVertexUnnamed)), 0);
    // `RulingBoundNotYet`: a run whose boundary ruling has no node on the interval's z-line — a
    // region crossing that line transversally in one sector while ending on it in another, which
    // the region emitter walks.
    assert_eq!(count(|c| *c == Rejected(RulingBoundNotYet)), 0);
    // `CylinderGateUndecided`, the chart's own refusal: a sector it cannot name as one arc is
    // read as the run of arcs it spans.
    assert_eq!(
        count(|c| *c == Rejected(RejectReason::CylinderGateUndecided)),
        0
    );
    // `First`: no family's first operation is refused — `corner-lo`'s reflex sector has its
    // winding read at the region's extremum, not at a node the ring bulges past.
    assert_eq!(count(|c| *c == First), 0);
    assert_eq!(count(|c| *c == Empty), 5);
}

/// **The ⊥ road runs on a panel and on a chain, and says what the walk-through predicts.**
/// Locked through the road's own ledger (`cycle_probe`) — the carved extents and spans are an
/// intermediate no result shows — and every class is traced by hand (`trace_every_class`), since
/// the production driver stops at the first decline and, in serial mode, would never reach the
/// classes this reads. Stations are axis parameters from the
/// boss's base (z = −1), so `t = z + 1`.
///
/// Panel (wall +x Cut, cut by the wall slab): at `z = 0.5` the class runs strictly inside the
/// panel's range — outer `Crosses`, the two rulings crossed carve **one** extent (the outer half)
/// away, one span left;
/// at `z = 0` the class is the panel's own lower arc — no outer answer, one graze carved, one
/// span. Chain (half wall Fuse, cut by the far cube): `z = −1` is the whole rim — a graze, nothing
/// carved; `z = 0` is inside — `Crosses`, the plate-bottom arc carves a graze on the inner half,
/// two spans; `z = 1` is the chain's top — no outer answer, the boss-top arc grazes, one span.
#[test]
fn a_cycle_is_carved_on_its_classes() {
    crate::ledger::owned(|| {
        {
            let (mut m, plate, boss) = boss_family([4.0, 2.0, -1.0], 4.0);
            let out = boolean(&mut m, BoolKind::Cut, plate, boss).expect("the notch builds");
            m.rebuild_adjacency();
            let t = m.add_cuboid(
                Point3::from_array([3.0, -1.0, 0.5]),
                Point3::from_array([6.0, 5.0, 1.5]),
            );
            m.rebuild_adjacency();
            crate::arrangement::trace_every_class(&m, out[0], t)
                .expect("the panel's classes trace");
        }
        {
            let (mut m, plate, boss) = boss_family([2.0, 0.0, -1.0], 2.0);
            let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the half boss builds");
            m.rebuild_adjacency();
            let far = m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
            m.rebuild_adjacency();
            crate::arrangement::trace_every_class(&m, out[0], far)
                .expect("the chain's classes trace");
        }
    });
    use crate::arrangement::CylOnClass;
    // ★ This test's own records: the two traces above and nothing else, so a station's rows
    // are its fixture's and can be counted rather than merely found.
    let hits = crate::arrangement::cycle_probe::HITS.mine();
    // `kinds = [rims, chains, panels, holes]`; every recorded entry of the shape at the station
    // must read the same, and the station must have been reached. A graze's `body_above` is
    // written in the class's stored frame, so the rim check asks only for the kind of answer.
    let check = |kinds: [usize; 4],
                 t: f64,
                 outer: fn(Option<CylOnClass>) -> bool,
                 carved: usize,
                 spans: usize| {
        let rows: Vec<_> = hits
            .iter()
            .filter(|h| h.kinds == kinds && h.t == t)
            .collect();
        assert!(!rows.is_empty(), "no record for {kinds:?} at t = {t}");
        for h in rows {
            assert!(outer(h.outer), "{kinds:?} at t = {t}: {h:?}");
            assert_eq!(
                (h.carved, h.spans),
                (carved, spans),
                "{kinds:?} at t = {t}: {h:?}"
            );
        }
    };
    let crosses = |o: Option<CylOnClass>| o == Some(CylOnClass::Crosses);
    let grazes = |o: Option<CylOnClass>| matches!(o, Some(CylOnClass::Grazes { .. }));
    let absent = |o: Option<CylOnClass>| o.is_none();
    let panel = [0, 0, 1, 0];
    check(panel, 1.5, crosses, 1, 1);
    check(panel, 1.0, absent, 1, 1);
    let chain = [1, 1, 0, 0];
    check(chain, 0.0, grazes, 0, 1);
    check(chain, 1.0, crosses, 1, 2);
    check(chain, 2.0, absent, 1, 1);
}

/// **The rulings road sweeps a chain rim**: on the half wall boss's wall class `y = 0`,
/// each ruling is `Transversal` from the whole rim at `z = −1` up to the plate's bottom `z = 0`
/// (both halves of the boss are there), then a `Graze` along the chain's own ruling edge up to
/// the boss's top `z = 1` — the collinear run whose arcs arrive from the outer half and leave into
/// the inner one toggles the face off above it. A lock of the **rule**, not of a result: every
/// re-operation of a half boss is still refused at its cap's chord, so
/// the classes are traced by hand and the road's ledger is read for this boss (`origin
/// (2, 0, −1)`, rulings spanning `t ∈ [0, 2]`).
#[test]
fn a_chain_sweeps_its_rulings() {
    let (mut m, plate, boss) = boss_family([2.0, 0.0, -1.0], 2.0);
    let out = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the half boss builds");
    m.rebuild_adjacency();
    let far = m.add_cuboid(Point3::from_array([20.0; 3]), Point3::from_array([21.0; 3]));
    m.rebuild_adjacency();
    crate::ledger::owned(|| {
        crate::arrangement::trace_every_class(&m, out[0], far).expect("the chain's classes trace")
    });
    // ★ This trace's own carvings. The origin and span still name the boss, because the trace
    // carves other classes of the same model too.
    let rows: Vec<Vec<crate::arrangement::SegKind>> = crate::arrangement::ruling_probe::CARVED
        .mine()
        .iter()
        .filter(|c| c.origin == [2.0, 0.0, -1.0] && c.span == [0.0, 2.0])
        .map(|c| c.kinds.clone())
        .collect();
    // Both rulings of the chain are swept, and only those two.
    assert_eq!(rows.len(), 2, "both rulings sweep: {rows:?}");
    for kinds in &rows {
        assert!(
            matches!(
                kinds[..],
                [
                    crate::arrangement::SegKind::Transversal { .. },
                    crate::arrangement::SegKind::Graze { .. }
                ]
            ),
            "a chain's ruling swept as {kinds:?}"
        );
    }
}

/// **The ring walk cuts a run where an arc departs, and the departure's side flanks the pieces.**
/// A hand-built ring on the unit cube's classes: `(0,0,0) → (1,0,0) → (1,1,0) → (1,1,1) →
/// (1,0,1) → (0,0,1)`, read against the class `x = 1` — four nodes on the line between two off it
/// (both on the `x < 1` side). With the edge `(1,1,0) → (1,1,1)` declared a departure to side σ,
/// the line is met in two runs, `[(1,0,0),(1,1,0)]` and `[(1,1,1),(1,0,1)]`: the first is flanked
/// by the off-line node and σ, the second by σ and the off-line node — so both pieces cross
/// exactly when σ is the *other* side, and negating σ swaps the answer. With the edge on the line
/// the four nodes are one run that touches and turns back. The closure supplies σ, so no
/// cylinder is needed: this locks the walk's rule, and `arc_departure_side`'s sign is locked by
/// the re-operation census (a wall boss's plate face is such a run with the arc's own σ).
#[test]
fn the_walk_cuts_a_run_at_a_departure() {
    use crate::combinatorics::{EdgeMeet, Feature, RingWalk, ring_against_plane, side_of};
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(Point3::from_array([5.0; 3]), Point3::from_array([6.0; 3]));
    m.rebuild_adjacency();
    let setup = crate::arrangement::plane_index_setup(&m, a, b).expect("setup");
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    // A class by the coordinate all three of its triangle's points share.
    let class = |axis: usize, at: f64| {
        (0..setup.geom.len())
            .find(|&c| {
                setup.geom[c]
                    .witness_coords()
                    .iter()
                    .all(|p| (p.as_array()[axis] - at).abs() < 1e-12)
            })
            .expect("a class of the unit cube")
    };
    let (x0, x1, y0, y1, z0, z1) = (
        class(0, 0.0),
        class(0, 1.0),
        class(1, 0.0),
        class(1, 1.0),
        class(2, 0.0),
        class(2, 1.0),
    );
    let ring = [
        NodeId::three_planes(Canon3::three([x0, y0, z0])),
        NodeId::three_planes(Canon3::three([x1, y0, z0])),
        NodeId::three_planes(Canon3::three([x1, y1, z0])),
        NodeId::three_planes(Canon3::three([x1, y1, z1])),
        NodeId::three_planes(Canon3::three([x1, y0, z1])),
        NodeId::three_planes(Canon3::three([x0, y0, z1])),
    ];
    let off = side_of(&jd, &[], ring[0], x1).expect("a side");
    assert_ne!(off, 0);
    let runs = |walk: RingWalk| -> Vec<(usize, usize, bool, i8)> {
        let RingWalk::Met(f) = walk else {
            panic!("the ring meets the line")
        };
        f.into_iter()
            .map(|f| match f {
                Feature::Run {
                    first,
                    len,
                    flanks_differ,
                    flank,
                } => (first, len, flanks_differ, flank),
                Feature::Crossing { .. } => panic!("no edge crosses x = 1 strictly"),
            })
            .collect()
    };
    for sigma in [off, -off] {
        let got = runs(ring_against_plane(&jd, &[], &ring, x1, |i| {
            Some(if i == 2 {
                EdgeMeet::Departs(sigma)
            } else {
                EdgeMeet::On
            })
        }));
        let crosses = sigma != off;
        assert_eq!(
            got,
            vec![(1, 2, crosses, off), (3, 2, crosses, sigma)],
            "σ = {sigma}, off-line side {off}"
        );
    }
    let got = runs(ring_against_plane(&jd, &[], &ring, x1, |_| {
        Some(EdgeMeet::On)
    }));
    assert_eq!(got, vec![(1, 4, false, off)]);
}

/// **The disk-side rule is derived, and the cells watch it** — the counterpart of the assertions
/// in `arrangement::disk_side_probe`, which fire where the fact is made and so see every chart in
/// the binary. Those say the geometry never contradicts the rule; what a test adds is the other
/// half: that the population is large, that the geometry actually speaks for most of it, and —
/// the one that would have caught the original defect — that the **`frame_sign` factor is
/// exercised**. A rule whose deciding factor is constant over the corpus is a rule nothing has
/// tested, which is exactly how it shipped wrong the first time.
///
/// ★ Read as floors, not literals. The rows are this test's own, but some of its operations are
/// **refused**, and the two builds do not do the same amount of work on a refused input: the
/// parallel one evaluates the remaining classes, the serial one stops at the first error. So the
/// shape of the population is held, not its size.
#[test]
fn the_disk_side_rule_is_derived_and_the_cells_watch_it() {
    // ★ **The test builds the population it measures**, and owns it: the very row this needs is a
    // `Reversed` root face — the **notch** a half boss leaves, whose ceiling class carries
    // `frame_sign = -1` — and a proportion read off the binary's rows would be diluted by
    // whatever else was running. So it runs the boss corpus itself (each family's first op, then
    // a slab cut of the result) and reads its own rows.
    crate::ledger::owned(|| {
        for (_, base, h) in BOSS_FAMILIES {
            for kind in [BoolKind::Fuse, BoolKind::Cut] {
                let (mut m, plate, boss) = boss_family(base, h);
                let Result::Ok(out) = boolean(&mut m, kind, plate, boss) else {
                    continue;
                };
                // The second operation is what puts a `Reversed` root face on a cut circle's
                // class (a notch's ceiling), which is where the frame factor decides.
                m.rebuild_adjacency();
                let t = m.add_cuboid(
                    Point3::from_array([-1.0, -1.0, 0.5]),
                    Point3::from_array([5.0, 5.0, 1.5]),
                );
                m.rebuild_adjacency();
                let _ = boolean(&mut m, BoolKind::Cut, out[0], t);
            }
        }
    });
    let rows = crate::arrangement::disk_side_probe::ROWS.mine();
    assert!(
        !rows.is_empty(),
        "the fixture produced no arc labels at all"
    );
    let sh = *crate::arrangement::cyl_chart::probe::shadow::COUNTS
        .lock()
        .unwrap();
    eprintln!(
        "SHADOW agree {} disagree {} vert_only {} horiz_only {} both_silent {} one_ruling {}",
        sh.0, sh.1, sh.2, sh.3, sh.4, sh.5
    );
    let checked = rows.iter().filter(|r| r.checked).count();
    let both = rows.iter().filter(|r| r.both_spoke).count();
    let frame_neg = rows.iter().filter(|r| r.frame_negative).count();
    // The geometry watches most of the population, so the agreement assertion is not a formality.
    assert!(
        checked * 2 > rows.len(),
        "the cells named the side for only {checked} of {} arcs",
        rows.len()
    );
    // An arc whose two cells both name a side is where "never the same side" is really tested.
    assert!(
        both > 0,
        "no arc had a witness on both sides — that premise check never ran"
    );
    // ★ The factor must decide some arcs, or this corpus cannot tell the rule from one without
    // it.
    assert!(
        frame_neg > 0,
        "no class with frame_sign = -1 reached the rule — the factor is untested here"
    );
}

/// **The corner boss builds, and only the extremum road makes it.** Its axis sits exactly on
/// **both** plate walls, so the chords through its circle are radii and the sector outside the
/// plate is reflex at the centre — the ring's smallest **node**, with the arc reaching further.
/// Reading the turn there inverted the sign, the void seed landed on a bounded sector, and one
/// solid's material was flipped across the whole component; the emitter then refused every one of
/// this family's booleans.
///
/// ★★★★ **The lock is the outcome, not a counter.** A process-global counter snapshotted around
/// this body would not measure it: every other test's booleans write to the same counter in
/// parallel, so the delta is not this population's. The behaviour is the honest lock, and it is a real one:
/// reverting `arc_extremum_winding`'s early return puts all three kinds back to
/// `CylinderGateUndecided` here and in three censuses.
#[test]
fn a_ring_whose_arc_bulges_past_its_nodes_reads_the_winding_there() {
    for (kind, want) in [
        (
            BoolKind::Fuse,
            first_volume(BoolKind::Fuse, [0.0, 0.0, -1.0], 4.0),
        ),
        (
            BoolKind::Cut,
            first_volume(BoolKind::Cut, [0.0, 0.0, -1.0], 4.0),
        ),
        (
            BoolKind::Common,
            first_volume(BoolKind::Common, [0.0, 0.0, -1.0], 4.0),
        ),
    ] {
        let (mut m, plate, boss) = boss_family([0.0, 0.0, -1.0], 4.0);
        let out = boolean(&mut m, kind, plate, boss).expect("the corner boss builds");
        m.rebuild_adjacency();
        let issues = nacre_validate::validate(&m);
        assert!(issues.is_empty(), "{kind:?}: {issues:?}");
        let v: f64 = out
            .iter()
            .map(|&sh| nacre_props::mass_props(&m, sh).expect("props").volume)
            .sum();
        assert!((v - want).abs() < 1e-9, "{kind:?}: volume {v} ≠ {want}");
    }
}

/// ★★★★★ **A circle names witnesses on its own rim, so a disk has something to offer the
/// converse.**
///
/// `inside`'s converse pass skips `Witness::In` without setting a flag, so a `Cell::Disk` in the
/// reverse position offers **only** its rim — without it every such descent would end in
/// [`crate::RejectReason::RingHasNoWitness`]. A nested pair of disks is the smallest shape that
/// walks that road: the small disk's centre lands inside the big one and the engine asks the
/// converse.
///
/// Both adapters answer disk↔disk with `disk_in_disk` before the engine sees it
/// (`nesting::cell_in_cell`, `boolean`'s merge road), so this configuration reaches
/// `cell_inside` **only from a test** — which is exactly why the road needs one.
///
/// With the rim in the supply the question never gets that far: a boundary witness that lands
/// strictly inside settles it at `Said::In` with no descent at all.
#[test]
fn a_disk_inside_a_disk_is_decided_by_its_rim() {
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([-10.0, -10.0, 0.0]),
        Point3::from_array([10.0, 10.0, 5.0]),
    );
    m.rebuild_adjacency();
    let faces_tab = collect_planes(&m, s).unwrap();
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, _plane_ix, _cyls) = dense_planes(&faces_tab, &canon);
    let jd = crate::planes::test_judge(&planes);
    let zero = nacre_exact::Rat::from_int(0);
    // A class the z axis actually meets — a cap, not a wall.
    let wc = (0..planes.len())
        .find(|&c| {
            combinatorics::class_coeffs_rat(&jd, c)
                .is_some_and(|k| k[0] == zero && k[1] == zero && k[2] != zero)
        })
        .expect("a cap class");
    let disk = |r: i128| {
        let q = nacre_exact::Rat::from_int;
        nacre_topo::CylinderDef::new(
            [zero, zero, zero],
            [zero, zero, q(1)],
            [q(1), zero, zero],
            nacre_exact::BigRat::from(q(r * r)), // the radius, stated as its square
        )
        .expect("a coaxial bore")
    };
    let (small, big) = (disk(2), disk(7));
    assert!(
        crate::nesting::cell_inside(
            &jd,
            NO_CYLS,
            wc,
            crate::nesting::Cell::Disk(&small),
            crate::nesting::Cell::Disk(&big),
        )
        .expect("the rim decides where the centre could only ask the converse"),
        "the small disk is inside the big one"
    );
    // And the other way round is a plain «out», by the same rim.
    assert!(
        !crate::nesting::cell_inside(
            &jd,
            NO_CYLS,
            wc,
            crate::nesting::Cell::Disk(&big),
            crate::nesting::Cell::Disk(&small),
        )
        .expect("decided"),
        "the big disk is not inside the small one"
    );
}

/// ★★★★★ **A circle cell is a circle only where the class is ⊥ to the axis.**
///
/// `Cell::Disk` has four rim witnesses, and their being *on the cell's boundary* rests on
/// the class plane being perpendicular to the cylinder's axis — a premise nothing enforced. Off it
/// the section is an ellipse: `û ⊥ axis` does not give `û ⊥ n`, so two of the four rim points leave
/// the plane, and `inside` answers a boundary `In` **without the converse**. That is a silent wrong
/// answer, which is why the two consumers ask.
///
/// The fixture is the smallest oblique one that still forms everything: a cuboid's `x̂` class with
/// `dir = (3,4,0)`, `ref_dir = (0,0,1)`. It is oblique (`n × m = (0,0,4) ≠ 0`), the centre still
/// forms (`n·m = 3 ≠ 0`), and the frame still stands (`‖dir‖ = 5`, `‖ref_dir‖ = 1`, both rational)
/// — so `û₁ = (0,0,1)` lies in the plane while `û₂ = (4,−3,0)/5` does not. ★ An axis like `(1,1,0)`
/// measures nothing here: `‖dir‖² = 2` is irrational, no frame forms, and the rim is absent with or
/// without the guard.
#[test]
fn an_oblique_class_carries_no_circle_and_so_offers_no_rim() {
    let q = nacre_exact::Rat::from_int;
    let zero = q(0);
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([-10.0, -10.0, 0.0]),
        Point3::from_array([10.0, 10.0, 5.0]),
    );
    m.rebuild_adjacency();
    let faces_tab = collect_planes(&m, s).unwrap();
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, _plane_ix, _cyls) = dense_planes(&faces_tab, &canon);
    let jd = crate::planes::test_judge(&planes);
    let class = |pick: &dyn Fn(&[nacre_exact::Rat; 4]) -> bool| {
        (0..planes.len())
            .find(|&c| combinatorics::class_coeffs_rat(&jd, c).is_some_and(|k| pick(&k)))
            .expect("the cuboid has all three")
    };
    let wall = class(&|k: &[nacre_exact::Rat; 4]| k[1] == zero && k[2] == zero && k[0] != zero);
    let cap = class(&|k: &[nacre_exact::Rat; 4]| k[0] == zero && k[1] == zero && k[2] != zero);
    let def =
        |o: [nacre_exact::Rat; 3], d: [nacre_exact::Rat; 3], e: [nacre_exact::Rat; 3], r: i128| {
            nacre_topo::CylinderDef::new(o, d, e, nacre_exact::BigRat::from(q(r * r)))
                .expect("a statable cylinder")
        };

    // ⊥: the four rim witnesses stand.
    let upright = def([zero; 3], [zero, zero, q(1)], [q(1), zero, zero], 2);
    assert_eq!(crate::nesting::rim_witness_count(&jd, cap, &upright), 4);
    // Oblique: the same cylinder statement, asked on a class its axis is not perpendicular to.
    // The frame forms and the centre forms — only the premise fails, and only the rim goes.
    let slanted = def([zero; 3], [q(3), q(4), zero], [zero, zero, q(1)], 2);
    assert!(
        nacre_exact::cyl_unit_frame(&slanted.dir(), &slanted.ref_dir()).is_some(),
        "the fixture must be one where the frame stands, or it measures nothing"
    );
    assert!(combinatorics::circle_centre_rat(&jd, wall, &slanted).is_some());
    assert_eq!(
        crate::nesting::rim_witness_count(&jd, wall, &slanted),
        0,
        "an ellipse's boundary is not `centre ± r·û`"
    );
}

/// ★★★★★ **And radii do not decide an ellipse.** `disk_in_disk` compares `(r_b − r_a)²`
/// with the centre distance, which is the right question only when both sections are circles.
///
/// This fixture is the half where it is **silently wrong**, not merely conservative: both axes are
/// `(3,4,0)`, so each section is an ellipse elongated along `ŷ` (semi-minor `r`, semi-major `5r/3`),
/// and the two centres are separated by `4/3` **along that major axis**. Shrinking `ŷ` by `3/5`
/// turns both into circles and the offset into `4/5`, so `4/5 + 1 ≤ 2` — the small disk really is
/// inside the big one. Comparing radii would answer `Ok(false)` — `(2−1)² = 1` is not
/// `> (4/3)² = 16/9` — so an oblique class is refused by name (`ObliqueCircleClass`).
///
/// ⚠ Concentric or minor-axis offsets are answered **correctly** even on an oblique class, so a
/// fixture that did not pin the offset's direction would show the guard "answering" rather than
/// the radius comparison being wrong.
#[test]
fn an_oblique_class_refuses_two_disks_rather_than_comparing_radii() {
    let q = nacre_exact::Rat::from_int;
    let zero = q(0);
    let mut m = Model::new();
    let s = m.add_cuboid(
        Point3::from_array([-10.0, -10.0, 0.0]),
        Point3::from_array([10.0, 10.0, 5.0]),
    );
    m.rebuild_adjacency();
    let faces_tab = collect_planes(&m, s).unwrap();
    let canon = plane_classes(&crate::planes::test_judge(&faces_tab));
    let (planes, _plane_ix, _cyls) = dense_planes(&faces_tab, &canon);
    let jd = crate::planes::test_judge(&planes);
    let wall = (0..planes.len())
        .find(|&c| {
            combinatorics::class_coeffs_rat(&jd, c)
                .is_some_and(|k| k[1] == q(0) && k[2] == q(0) && k[0] != q(0))
        })
        .expect("an x-normal wall class");
    // Same axis for both, so the two section centres differ by exactly the origins' difference.
    let axis = [q(3), q(4), zero];
    let seam = [zero, zero, q(1)];
    let small =
        nacre_topo::CylinderDef::new([zero; 3], axis, seam, nacre_exact::BigRat::from(q(1)))
            .expect("small");
    let big = nacre_topo::CylinderDef::new(
        [zero, nacre_exact::Rat::new(4, 3).unwrap(), zero],
        axis,
        seam,
        nacre_exact::BigRat::from(q(4)), // radius 2, as r²
    )
    .expect("big");
    // The offset really is `4/3` along `ŷ`, in the class — stated so the fixture cannot drift.
    let (ca, cb) = (
        combinatorics::circle_centre_rat(&jd, wall, &small).expect("centre a"),
        combinatorics::circle_centre_rat(&jd, wall, &big).expect("centre b"),
    );
    let delta: Vec<_> = (0..3)
        .map(|k| cb[k].checked_sub(ca[k]).expect("a rational offset"))
        .collect();
    assert_eq!(
        delta,
        vec![zero, nacre_exact::Rat::new(4, 3).unwrap(), zero]
    );
    let err = crate::nesting::disk_in_disk(&jd, wall, &small, &big)
        .expect_err("radii cannot decide an ellipse");
    assert!(
        matches!(
            err,
            BoolError::Rejected {
                reason: crate::RejectReason::ObliqueCircleClass,
                ..
            }
        ),
        "{err:?}"
    );
}

/// ★★★★★ **One rule for "the point this edge names", with every arm — the plainest included.**
///
/// Spelled as producers per supply — chained verbatim in one, counted again in its instrument,
/// absent on the component road one dimension up — a planar component whose every corner grazes
/// runs out of witnesses and refuses `NoClearRay` where the shape's own name is
/// `SelfTouchingResult` (`contact_separates.rs`'s ⑤c). "A chord's midpoint" and "an edge's
/// interior point" are one rule — *a point inside a straight edge* — and it must cover **the
/// plainest edge there is**, two three-plane corners joined by a straight step, which an arm
/// that starts by asking for a pierce name misses.
///
/// The oracle is computed **from the corners**, not by calling the arm a second time.
#[test]
fn a_straight_edge_between_rational_corners_names_its_own_midpoint() {
    let (planes, p, outer, _hole) = holed_face_rings("dimple");
    let jd = crate::planes::test_judge(&planes);
    let ring = combinatorics::ring_from_names(p, &outer).unwrap();
    let half = nacre_exact::Rat::new(1, 2).unwrap();
    let mut seen = 0usize;
    for e in &ring {
        let got: Vec<_> = combinatorics::edge_interior_points(&jd, NO_CYLS, e).collect();
        // Exactly one arm answers a rational-cornered straight edge: a careless fourth arm, or one
        // whose guards overlap, shows up here.
        assert_eq!(got.len(), 1, "one point per plain edge: {got:?}");
        let (a, b) = (
            combinatorics::node_coords_rat(&jd, e.node).expect("a rational corner"),
            combinatorics::node_coords_rat(&jd, e.to).expect("a rational corner"),
        );
        let want: Vec<_> = (0..3)
            .map(|k| a[k].checked_add(b[k]).unwrap().checked_mul(half).unwrap())
            .collect();
        assert_eq!(got[0].to_vec(), want, "the midpoint of the edge's two ends");
        // And it is *on* the edge, which is the invariant the component road's `Probe` rests on:
        // strictly between the ends in every coordinate that separates them.
        for k in 0..3 {
            if a[k] != b[k] {
                let (lo, hi) = if a[k] < b[k] {
                    (a[k], b[k])
                } else {
                    (b[k], a[k])
                };
                assert!(
                    lo < got[0][k] && got[0][k] < hi,
                    "strictly inside in axis {k}"
                );
                seen += 1;
            }
        }
    }
    assert!(
        seen > 0,
        "the fixture must have edges that separate coordinates"
    );
}

/// ★★★★★ **A plain bored cube builds at any size, and the answer is right.**
///
/// A bore through a cube is the least exotic input this kernel has, and it must not be
/// **refused** (`WitnessNotRational`) because a coordinate needs a long decimal at a small scale.
/// With a chart normal the four-coefficient canonicalisation leaves non-primitive, the plane
/// offset's denominator passes `sqrt(i128)` between `1/3·1e-3` and `1/3·1e-4` (measured), with
/// nothing about the geometry changing across that line; the chart reads a primitive normal
/// ([`crate::combinatorics`]'s `primitive_normal`).
///
/// The sizes below span **seven orders of magnitude**, each with a full 17-digit mantissa so the
/// denominators are as bad as f64 can make them; the last two are the regression half. The oracle
/// is `s³ − πr²s`.
///
/// The **radius** is `s/4` where the sketch can state it and a short decimal below that: a sketched
/// circle carries `r²` as a `Rat`, and a full-digit radius under `1e-4` squares past `i128`
/// (`SketchError::Undecidable` — the application's road refuses the circle before any solid
/// exists). The widths this test is about are the cube's and the centre's, which stay full-digit.
///
/// ⚠ **The tolerance is relative and it belongs to the oracle, not to the kernel.** At `s = 2e-7`
/// the volume is ~6.4e-21, where any absolute epsilon is meaningless; and `π` does not cancel
/// here, so the oracle is f64 arithmetic while the kernel's answer is not. Precedent for the
/// relative form: `pocket_on_a_random_slanted_face_is_valid_or_rejects`' `1e-9 * 8.0`.
#[test]
fn a_bored_cube_builds_at_any_size_and_its_volume_is_right() {
    for (s, r) in [
        (2.0000000000000003e-7, 5e-8),
        ((1.0 / 3.0) * 1e-4, 8e-6),
        ((1.0 / 7.0) * 1e-3, 3.5e-5),
        (1.0 / 3.0, 1.0 / 12.0),
        (2.0000000000000003, 2.0000000000000003 / 4.0),
    ] {
        let h = s / 2.0;
        let mut m = Model::new();
        let block = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([s; 3]));
        let bore = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([h, h, -h]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            r,
            s * 2.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Cut, block, bore)
            .unwrap_or_else(|e| panic!("s = {s:e}: {e:?}"));
        assert_eq!(out.len(), 1, "s = {s:e}: a bored cube is one body");
        assert!(
            nacre_validate::validate(&m).is_empty(),
            "s = {s:e}: {:?}",
            nacre_validate::validate(&m)
        );
        let vol = nacre_props::mass_props(&m, out[0])
            .expect("mass props")
            .volume;
        let want = s * s * s - std::f64::consts::PI * r * r * s;
        assert!(
            (vol - want).abs() <= 1e-12 * want,
            "s = {s:e}: volume {vol:e} against {want:e}"
        );
    }
}

/// ★ **The rim witness is the point the cylinder's own statement names.** The extrude
/// road states a unit, perpendicular frame (a sketch frame is checked orthonormal exactly, and a
/// whole circle's `ref_dir` is a rim chord divided by its own radius), so the unit frame comes
/// back as the statement itself and `centre + r·û₁` is that circle's seam point — not new
/// geometry, just the point the arrangement's uncut circle drops the node for.
#[test]
fn a_rim_witness_is_the_statements_own_seam_point() {
    let mut m = Model::new();
    let op = cylinder_op(&m, [3.0, -4.0], 2.5, 6.0);
    apply(&mut m, &op).expect("the cylinder extrudes");
    let def = lone_cylinder_def(&m);
    let (u1, u2) = nacre_exact::cyl_unit_frame(&def.dir(), &def.ref_dir())
        .expect("the extrude road's frame is unit and perpendicular");
    assert_eq!(u1, def.ref_dir(), "û₁ is the statement's own ref_dir");
    assert_eq!(
        u2,
        nacre_exact::cross3_rat(&def.dir(), &def.ref_dir()).expect("m × ê"),
        "û₂ is m × ê, already unit"
    );
    // The post-condition the supply asserts where it mints: a rim point is on the rim, exactly.
    let r = def
        .radius_exact()
        .expect("a stated radius has a rational spelling");
    for u in [u1, u2] {
        let mut p = def.origin();
        for k in 0..3 {
            p[k] = p[k].checked_add(r.checked_mul(u[k]).unwrap()).unwrap();
        }
        assert_eq!(
            nacre_exact::quad::cylinder_radial_side(&p, &def.origin(), &def.dir(), def.r2()),
            nacre_exact::Orient::Zero,
        );
    }
}

/// ★ **A ring with no three-plane corner is answered by the witnesses it does have.**
/// A nesting question is answered by a witness of the source cell, and the supply that names
/// those witnesses is one spelling for every arm. With every corner filleted the plate's cap ring
/// offers eight pierce corners, four edge-interior points and **no** three-plane name at all — an
/// arm that read only three-plane names would find nothing, and the plate could not enter any
/// boolean, whatever the other solid was (this one is a hundred units away).
///
/// The audit is also the negative control: it must *see* that shape, or its zeros mean nothing.
#[test]
fn a_ring_with_no_three_plane_corner_is_answered_by_the_witnesses_it_has() {
    let mut m = Model::new();
    let plate = rounded_plate(&mut m, 4, 4);
    let far = m.add_cuboid(
        Point3::from_array([200.0, -5.0, -5.0]),
        Point3::from_array([210.0, 5.0, 5.0]),
    );
    m.rebuild_adjacency();
    // The probe records only work a test owns, and `mine` then reads this fixture's questions
    // alone — no other test's boolean can put a row in this list.
    let out = crate::ledger::owned(|| boolean(&mut m, BoolKind::Fuse, plate, far));
    let rows = nesting::nesting_probe::ROWS.mine();

    // Two bodies: the plate and the far box, each untouched.
    let out = out.expect("the rounded plate enters a boolean");
    assert_eq!(out.len(), 2, "the plate and the box, apart");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    let want = (90.0 * 50.0 - 4.0 * (25.0 - std::f64::consts::PI * 25.0 / 4.0)) * 12.0
        - 4.0 * std::f64::consts::PI * 3.5 * 3.5 * 12.0
        + 1000.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
    // ★ The negative control: the instrument must *see* the failing shape, or its zeros mean
    // nothing. A ring with no three-plane name, asked against a disk.
    // ★ `route` first. A row from the shared-node or two-radii road
    // describes an offer nobody read, and 85% of the rows are of that kind.
    let nameless_vs_disk: Vec<_> = rows
        .iter()
        .filter(|r| r.route == nesting::Route::Engine && !r.a_disk && r.named == 0 && r.b_disk)
        .collect();
    // ★ Counted on this boolean's own questions (224 of them): 24 ask the cap ring — which has
    // no three-plane corner — against a disk. That the shape is reached at all is the negative
    // control; how the 24 divide among the plate's four bores is not a claim this makes.
    assert_eq!(
        nameless_vs_disk.len(),
        24,
        "the audit must see the ring that has no three-plane corner asked against a disk: {rows:?}"
    );
    // And those rings do have exact witnesses — of a kind this arm does not ask for.
    for r in &nameless_vs_disk {
        assert_eq!(r.coords, 0, "no three-plane coordinate: {r:?}");
        assert!(
            r.pierce > 0,
            "but every corner is a rational pierce point: {r:?}"
        );
        // The rest of the supply, measured: a rounded outline's corners are tangencies, so no
        // edge is a whole chord (`chord` 0) and the ring is not a circle — but its straight edges
        // *do* name interior points, which is one more exact witness this arm never asks for.
        assert_eq!(r.chord, 0, "no whole chord: {r:?}");
        assert!(r.edge > 0, "but its edges name interior points: {r:?}");
        assert!(!r.circle, "and the ring is not a circle: {r:?}");
        assert!(!r.b_mixed, "a disk target is not a mixed ring: {r:?}");
    }
    // ★ A disk source describes its own supply. Before the rim
    // it read all zeros in a ring's vocabulary, so no row ever said what that arm had to offer —
    // which is how a population counted off these rows came out wrong.
    let disks: Vec<_> = rows
        .iter()
        .filter(|r| r.route == nesting::Route::Engine && r.a_disk)
        .collect();
    assert_eq!(disks.len(), 40, "the plate's bores ask as disks: {rows:?}");
    for r in &disks {
        assert_eq!(r.rim, 4, "a disk offers its four rim witnesses: {r:?}");
    }
}

/// ★ **The premise, measured**: where both roads can answer one question, they agree.
/// The engine picks whichever witness comes first, so «any witness decides» has to
/// be a fact and not a hope. A single disagreement would mean the answer depends on the witness,
/// which is the one thing that would make the unification wrong.
#[test]
fn the_two_roads_never_disagree() {
    crate::ledger::owned(|| {
        for (fillets, bores) in [(3, 4), (4, 0), (3, 0)] {
            for k in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
                let mut m = Model::new();
                let plate = rounded_plate(&mut m, fillets, bores);
                let boss = m.add_cuboid(
                    Point3::from_array([10.0, -20.0, 12.0]),
                    Point3::from_array([25.0, 20.0, 62.0]),
                );
                m.rebuild_adjacency();
                let _ = boolean(&mut m, k, plate, boss);
            }
        }
    });
    let rows = nesting::nesting_probe::ROWS.mine();
    // ★ `route` first here too. A shared-node pair is «not comparable» — the engine
    // discards that question — so agreeing about it measures nothing the kernel acts on.
    let both: Vec<_> = rows
        .iter()
        .filter(|r| r.route == nesting::Route::Engine)
        .filter_map(|r| r.roads)
        .collect();
    // ★ A floor, not a count: some of these operations are refused, and the parallel build
    // evaluates the classes after the first error while the serial one stops there — so how many
    // questions a refused case asks is a property of the build, not of the kernel.
    assert!(
        !both.is_empty(),
        "the fixtures must reach questions both roads can answer"
    );
    let disagree = both.iter().filter(|(a, b)| a != b).count();
    assert_eq!(
        disagree,
        0,
        "the name road and the coordinate road answered differently on {disagree} of {} questions",
        both.len()
    );
}

/// ★ **An island standing inside a round hole**, a shape the corpus did not hold: two
/// bodies whose top faces are coplanar, one of them a ring that lies inside the other's circular
/// hole and covers its centre.
///
/// ★★ **It does not walk the merge road's converse clause** — measured, and recorded rather than
/// dressed up. That road has a converse clause («a centre inside a ring does not mean the disk is
/// inside it»), and this is the shape that would make a spelling without it claim two owners for
/// one circle. It passes without the clause too:
/// `unify_coplanar_faces` looks for a circle's owner **within its own edge-connected group**, and
/// an island that does not touch the rim is not in that group. For a competing outer to be in it,
/// something must reach the disk's face without cutting its rim — which no producer here makes.
/// So the converse clause in that road guards a population that cannot yet arise, and this test
/// locks the shape rather than the clause.
#[test]
fn an_island_inside_a_round_hole_does_not_claim_the_hole() {
    let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
    let mut m = Model::new();
    let profile = stated(vec![
        line(p2(-20.0, -20.0), p2(20.0, -20.0)),
        line(p2(20.0, -20.0), p2(20.0, 20.0)),
        line(p2(20.0, 20.0), p2(-20.0, 20.0)),
        line(p2(-20.0, 20.0), p2(-20.0, -20.0)),
        stated::circle(p2(0.0, 0.0), 10.0),
    ])
    .unwrap()
    .remove(0);
    let frame = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: plate, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile,
            dist: 5.0,
        },
    )
    .expect("the holed plate extrudes") else {
        unreachable!()
    };
    let island = m.add_cuboid(
        Point3::from_array([-3.0, -3.0, 0.0]),
        Point3::from_array([3.0, 3.0, 5.0]),
    );
    m.rebuild_adjacency();
    let out = boolean(&mut m, BoolKind::Fuse, plate, island).expect("the two bodies fuse");
    assert_eq!(out.len(), 2, "they do not touch, so they stay two");
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "{:?}",
        nacre_validate::validate(&m)
    );
    let v: f64 = out
        .iter()
        .map(|&s| nacre_props::mass_props(&m, s).expect("props").volume)
        .sum();
    let want = (40.0 * 40.0 - std::f64::consts::PI * 100.0) * 5.0 + 6.0 * 6.0 * 5.0;
    assert!((v - want).abs() < 1e-9, "{v} vs {want}");
}
