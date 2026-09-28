//! Cylinders in the arrangement: the ruling trace, lateral and mixed parity, cut circles, caps,
//! seated circles.

use super::*;

/// ★★ **The production path, run on a cylinder-bearing model.**
///
/// A `.plane()` call on the production road with cylinder rows flowing into it —
/// the `work` class list, the decline witness, the report's `class_of` — aborts the
/// kernel on the first drill, and reading the ~30 call sites does not find one. **Running
/// the road is what finds them**, which is why this test exists.
///
/// What it asserts is deliberately weak on geometry and strong on survival: the tracer
/// completes, every face it emits is still a plane face, and the
/// drilled cap carries its circular hole.
#[test]
fn the_production_road_survives_a_cylinder_past_the_stopper() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // a through-hole: circles on both box caps
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        class_owner,
        n_a,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let trace_in = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &cyls,
        Default::default(),
    );
    let (faces, _, _) = trace_result_faces(
        &m,
        BoolKind::Cut,
        a,
        b,
        &jd,
        &faces_tab,
        &plane_ix,
        &cyls,
        n_a,
        &class_owner,
        // `Proved` on purpose: the reuse guard must be what turns the shortcut off, not the
        // caller. Without it `pass_through` would carry planes across and drop the lateral.
        crate::reuse::ClassReuse::Proved,
        &trace_in,
    )
    .expect("the drill population traces");
    assert!(
        faces.iter().all(|f| matches!(f.surf, ClassIx::Plane(_))),
        "the plane arrangement emits plane faces only; bands are C4b-2"
    );
    // The box: 4 walls + 2 caps, and each cap carries the drill's circular hole.
    assert_eq!(faces.len(), 6, "{faces:?}");
    let holed = faces
        .iter()
        .filter(|f| {
            f.inner
                .iter()
                .any(|b| matches!(b, crate::draft::Bound::Circle { .. }))
        })
        .count();
    assert_eq!(holed, 2, "both caps are drilled: {faces:?}");
}

/// The rulings-road harness: the through-boss geometry
/// — a plate and a cylinder overlapping its full height, axis exactly
/// on the `x = 40` wall plane — assembled **past the standing gate** from production parts:
/// `plane_index_setup_inner` (the gate-free half) plus the cylinder table built the way the
/// gate builds it. Returns everything the two locks below read, and the record the
/// gate produces: `(wall class, cylinder)` marked not-proven-clear.
#[allow(clippy::type_complexity)]
fn armed_through_boss() -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::arrangement::PlaneSetup,
    usize,
    std::collections::HashSet<(usize, usize)>,
) {
    armed_through_boss_z(-10.0, 50.0)
}

/// [`armed_through_boss`] with the boss's axial extent chosen: `z_lo` and height `h`. The
/// default runs through the plate (`−10`, `50`); a boss whose lower cap sits **inside** the
/// plate (`10`, `30`) has a lateral whose lower boundary is a **chain** — arcs at z = 10
/// (outside the plate) and z = 20 (inside) joined by rulings — the staircase the corner
/// rule is watched on.
#[allow(clippy::type_complexity)]
fn armed_through_boss_z(
    z_lo: f64,
    h: f64,
) -> (
    Model,
    Handle<Solid>,
    Handle<Solid>,
    crate::arrangement::PlaneSetup,
    usize,
    std::collections::HashSet<(usize, usize)>,
) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 20.0]),
    );
    let boss = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([40.0, 20.0, z_lo]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        5.0,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    let (mut setup, cyl_surfs) =
        crate::arrangement::plane_index_setup_inner(&m, plate, boss).unwrap();
    for &surf in &cyl_surfs {
        let nacre_topo::Surface::Cylinder { def, .. } = m.surface(surf) else {
            unreachable!("a cylinder row carries a cylinder truth")
        };
        let nacre_geom::Surface::Cylinder(cache) = m.surface_cache(surf) else {
            unreachable!("push_cylinder_raw pairs them, so a cylinder truth has a cylinder cache")
        };
        setup.cyls.push(crate::planes::WorkingCyl {
            surf,
            def: def.clone(),
            realized: *cache,
            owner: crate::planes::SolidSide::A,
        });
    }
    let wc = setup
        .geom
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[0] - 40.0).abs() < 1e-12)
        })
        .expect("the x = 40 wall class");
    let crossings: std::collections::HashSet<(usize, usize)> = [(wc, 0)].into_iter().collect();
    (m, plate, boss, setup, wc, crossings)
}

/// ★ **The gate's record arms the rulings road — and only the record.** With the
/// through-boss pair listed, the boss's trace on the wall class is the rectangle: two
/// rulings (one per side, Pierce-named ends, no plain segments) and two cap chords
/// (`Lo`/`Hi` roots on the cap classes). With the record empty — every production call
/// today — the same trace is empty: the negative control that pins "an empty record
/// changes nothing", which is the very population an unconditionally-firing arm broke
/// (14 arc-family tests red, measured).
#[test]
fn the_gates_record_arms_the_ruling_trace() {
    let (m, _plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut tr = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        crossings.clone(),
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(
        tr.segs.len(),
        2,
        "one chord per cap — a segment of the sweep since cell ⑩"
    );
    let mut sides: Vec<i8> = tr.rulings.iter().map(|r| r.side).collect();
    sides.sort_unstable();
    assert_eq!(sides, [-1, 1], "one ruling per side");
    for r in &tr.rulings {
        for nd in r.end {
            let (_, cyl, _) = combinatorics::pierce_name(nd).expect("a Pierce end");
            assert_eq!(cyl, 0);
        }
        assert!(matches!(r.kind, SegKind::Transversal { .. }));
    }
    for c in &tr.segs {
        assert_eq!(
            c.end_h,
            [combinatorics::EndPin::Cylinder; 2],
            "a chord's ends are the two roots"
        );
        // The chord rides the cap's own class — a ⊥ plane at z = −10 or z = 40.
        let z = setup.geom[c.wall].witness_coords()[0].as_array()[2];
        assert!(
            setup.geom[c.wall]
                .witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
                && ((z + 10.0).abs() < 1e-12 || (z - 40.0).abs() < 1e-12),
            "cap class at z = {z}"
        );
    }
    // The negative control: today's record.
    let mut tr0 = Trace::default();
    trace_one_of(
        &m,
        boss,
        SolidSide::B,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_b,
        &setup.plane_ix,
        Default::default(),
        &mut tr0,
    );
    assert!(
        tr0.rulings.is_empty() && tr0.segs.is_empty(),
        "empty record, empty road"
    );
}

/// ★ **The armed arrangement digests the rectangle** — production bricks end to end on the
/// wall class. Each ruling is cut at its T-junctions with the plate's `z = 0` and `z = 20`
/// lines (three pieces each), those lines split at the same Pierce nodes, the chords close
/// the far ends, and the walk closes the subdivision: one unbounded contour and five
/// bounded cells — plate-left, plate-right, the overlap band, and the rectangle's two
/// overhangs.
#[test]
fn the_armed_wall_class_walks_to_closed_cells() {
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let tr = trace_on_class_of(
        &m,
        plate,
        boss,
        wc,
        &jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_a,
        &setup.inc_b,
        &setup.plane_ix,
        crossings,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(&jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
    let split = drop_newsless(split).unwrap();
    let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
    // This fixture's class is parallel to the only cylinder in it, so it traces rulings and
    // no circle. Not a general law: a class may carry both, when the circle
    // belongs to a *different*, perpendicular cylinder.
    assert!(circles.is_empty(), "this class meets no ⊥ cylinder");
    let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
    assert_eq!(rulings.len(), 2);
    let edges = ClassEdges::of(
        &jd,
        &setup.cyls,
        wc,
        &split,
        &circles,
        &rulings,
        &Aliases::default(),
    )
    .unwrap();
    assert_eq!(
        edges.rulings.len(),
        6,
        "each ruling cut at z = 0 and z = 20"
    );
    let (cells, _face_of) = walk_cells(&jd, &setup.cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 6, "five bounded cells and the outer");
    assert_eq!(
        cells.iter().filter(|c| c.winding == -1).count(),
        1,
        "one connected skeleton, one outer contour"
    );
}

/// One armed class's arrangement, through the production bricks (the chain
/// `the_armed_wall_class_walks_to_closed_cells` spells out) — for the locks that need
/// several classes' products at once.
fn armed_class_edges<'a>(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &'a crate::arrangement::PlaneSetup,
    jd: &Judge<'a, WorkingPlane>,
    wc: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
) -> ClassEdges<'static> {
    let tr = trace_on_class_of(
        m,
        plate,
        boss,
        wc,
        jd,
        &setup.cyls,
        &setup.planes,
        &setup.surf_ix,
        &setup.inc_a,
        &setup.inc_b,
        &setup.plane_ix,
        crossings.clone(),
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    let merged = merge_coincident(jd, &tr.segs, wc, &Aliases::default());
    let split = split_at_crossings(jd, &setup.cyls, wc, &merged, &mut Aliases::default()).unwrap();
    let split = drop_newsless(split).unwrap();
    let circles = merge_circles(&tr.circles, &setup.cyls, &Aliases::default()).unwrap();
    let rulings = merge_rulings(&tr.rulings, &setup.cyls, &Aliases::default());
    let e = ClassEdges::of(
        jd,
        &setup.cyls,
        wc,
        &split,
        &circles,
        &rulings,
        &Aliases::default(),
    )
    .unwrap();
    // Owned copies so the borrows above may end — a test convenience, not a production shape.
    ClassEdges {
        segs: std::borrow::Cow::Owned(e.segs.into_owned()),
        arcs: std::borrow::Cow::Owned(e.arcs.into_owned()),
        rulings: std::borrow::Cow::Owned(e.rulings.into_owned()),
        circles: std::borrow::Cow::Owned(e.circles.into_owned()),
        cut_rims: e.cut_rims,
    }
}

/// **A lateral face answers the crossing question by the same parity, on its own chart**
/// — the digon oracle's lateral twin, on the through-boss Fuse's lateral. The
/// truth is known in world coordinates, so every crossing a rational ray makes with the
/// cylinder is checked: rays `{y = y₀, z = z₀}` over a half-step lattice (two roots each,
/// at x = 40 ± √(25 − (y₀ − 20)²) — irrational θ) and the **station column**
/// `{x = 40, z = z₀}`, whose roots are exactly the rulings' points (40, 15) and (40, 25):
/// on the ruling within the notch's z (a corner at its ends), on the face beyond it.
///
/// Two builds. **Through** (caps at −10 and 40): a band between the caps' whole circles
/// with the plate's notch as its one hole (the plate-side half-circle × z ∈ (0, 20), two
/// arcs on the cut rims and two rulings on the wall x = 40). **Staircase** (caps at 10 and
/// 40, the lower cap inside the plate): the lower boundary is a chain — the arc at z = 10
/// outside the plate, the arc at z = 20 inside, the two rulings between — so at the
/// station (40, 25) one arc ends `hi` and the other starts `lo`: the corner rule's one
/// discriminating population. ☑ In a notch both arcs share their ends, and «both ends
/// count» is invisible there — measured; the staircase is why the second build exists.
/// The other station, (40, 15), is the seam (the fixture puts it on −y): a
/// seam-incident root with an arc above it is the one tie the loops road keeps.
///
/// ☑ Measured on the through build's rays: on the two whole-circle rims alone, «opposite
/// sides of the two rim planes» and the loops road's «exactly one rim above» agree on all
/// 2,180 rays — every root, the graze on a rim included — which is why no banded arm
/// exists.
#[test]
fn the_lateral_parity_agrees_with_the_notch_it_bounds() {
    for staircase in [false, true] {
        lateral_lattice(staircase);
    }
}

fn lateral_lattice(staircase: bool) {
    use nacre_exact::quad::{CylinderMeet, plane_plane_cylinder, plane_side};
    use nacre_exact::{Orient, Rat};
    let caps = if staircase {
        [10.0, 40.0]
    } else {
        [-10.0, 40.0]
    };
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss_z(caps[0], caps[1] - caps[0]);
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, plane_faces, rows, _) =
        armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, caps);
    let fuse = crate::arrangement::cyl_chart::emit_lateral(
        BoolKind::Fuse,
        &jd,
        &setup.cyls,
        &plane_faces,
        &curved,
        &rows,
    )
    .unwrap();
    assert_eq!(fuse.len(), 1);
    let cf = crate::assembly::comp_face(&jd, &setup.cyls, &fuse[0]).unwrap();
    let combinatorics::CompSurf::Cylinder(def) = &cf.surf else {
        panic!("a lateral")
    };
    let combinatorics::BoundEdges::Lateral(loops) = &cf.outer else {
        panic!("loops, got {:?}", cf.outer)
    };
    assert!(cf.inner.is_empty(), "a lateral's holes are among its loops");
    let rims: Vec<usize> = loops
        .iter()
        .filter_map(|l| match l {
            combinatorics::LateralLoop::Circle(c) => Some(*c),
            combinatorics::LateralLoop::Ring(_) => None,
        })
        .collect();
    if staircase {
        assert_eq!(rims.len(), 1, "the upper cap's whole circle");
        assert_eq!(loops.len(), 2, "and the chain below");
    } else {
        assert_eq!(rims.len(), 2, "the caps' two whole circles");
        assert_eq!(loops.len(), 3, "and the notch as one ring");
    }
    let (o, mm, r2) = (def.origin(), def.dir(), def.r2());
    let rat = |k: i128| Rat::new(k, 2).unwrap();
    let ri = Rat::from_int;
    let zero = ri(0);
    let neg = |v: Rat| zero.checked_sub(v).unwrap();
    let x40 = [ri(1), zero, zero, neg(ri(40))];
    let y20 = [zero, ri(1), zero, neg(ri(20))];
    let (cap_lo, cap_hi) = (ri(caps[0] as i128), ri(caps[1] as i128));
    // The notch's z range: the plate's own (0, 20) through the boss, or — with the lower
    // cap inside the plate — the chain's step from the cap (10) to the plate's top (20).
    let (n_lo, n_hi) = (if staircase { cap_lo } else { zero }, ri(20));
    // The truth on the cylinder, by the root's side of the wall and the ray's z. `seam`:
    // the root is the station (40, 15), the seam generator; a seam-incident root with an
    // arc above it is the tie the loops road keeps (`arc_span`'s `SeamRoot`; z is asked
    // first, so with no arc above the rims still decide).
    let truth = |x_side: Orient, z0: Rat, seam: bool| -> Option<bool> {
        if seam && z0 < n_hi {
            return None;
        }
        if z0 >= cap_hi {
            return if z0 == cap_hi { None } else { Some(false) };
        }
        if z0 < cap_lo {
            return Some(false);
        }
        if z0 == cap_lo {
            // The lower cap: through the plate a whole rim (a tie everywhere); in the
            // staircase an arc outside the plate only — a tie there, a corner on the wall,
            // and nothing inside, where the face starts at the plate's top.
            return match (staircase, x_side) {
                (true, Orient::Negative) => Some(false),
                _ => None,
            };
        }
        match x_side {
            // Outside the plate: the band, whole.
            Orient::Positive => Some(true),
            // Inside the plate's footprint: the notch (a hole, or the chain's step), its
            // arcs the boundary — both arcs through the plate, the upper one alone in the
            // staircase (its lower boundary there is the cap, handled above).
            Orient::Negative => {
                if (!staircase && z0 == n_lo) || z0 == n_hi {
                    None
                } else {
                    Some(!(z0 > n_lo && z0 < n_hi))
                }
            }
            // On the wall: a ruling for z within the notch (corners at its ends), the
            // face beyond.
            Orient::Zero => {
                if z0 >= n_lo && z0 <= n_hi {
                    None
                } else {
                    Some(true)
                }
            }
        }
    };
    // on face, off (hole), boundary, tangent, station on-ruling, station on-face, seam ties
    let mut n = [0usize; 7];
    let mut ask = |pa: [Rat; 4], pb: [Rat; 4], z0: Rat, station: bool| {
        let roots = match plane_plane_cylinder(&pa, &pb, &o, &mm, r2).unwrap() {
            CylinderMeet::Pair { line, s } => (line, s),
            CylinderMeet::Tangent { .. } => {
                n[3] += 1;
                return;
            }
            other => panic!("an ordinary ray: {other:?}"),
        };
        let (line, s) = roots;
        for root in &s {
            let x_side = plane_side(&x40, &line, root);
            let seam = station && plane_side(&y20, &line, root) == Orient::Negative;
            let want = truth(x_side, z0, seam);
            let got = combinatorics::loop_parity(&jd, def, loops, &line, root);
            assert_eq!(
                got, want,
                "staircase {staircase} z0 {z0:?} side {x_side:?} station {station} seam {seam}"
            );
            match (got, station, seam) {
                (Some(true), false, _) => n[0] += 1,
                (Some(false), false, _) => n[1] += 1,
                (None, false, _) => n[2] += 1,
                (None, true, true) if z0 < n_lo || z0 > n_hi => n[6] += 1,
                (None, true, _) => n[4] += 1,
                (Some(_), true, _) => n[5] += 1,
            }
        }
    };
    for k in 0..=108 {
        let z0 = ri(-12).checked_add(rat(k)).unwrap();
        let pz = [zero, zero, ri(1), neg(z0)];
        for j in 0..=20 {
            let y0 = ri(15).checked_add(rat(j)).unwrap();
            ask([zero, ri(1), zero, neg(y0)], pz, z0, false);
        }
        ask(x40, pz, z0, true);
    }
    eprintln!(
        "lateral lattice (staircase {staircase}): on {} hole {} boundary {} tangent {} \
             station on-ruling {} station on-face {} seam ties {}",
        n[0], n[1], n[2], n[3], n[4], n[5], n[6]
    );
    assert!(
        n.iter().all(|&c| c > 0),
        "every arm has a population: {n:?}"
    );
}

/// ★ **A cut circle is a band boundary**. On the armed
/// through-boss, the per-class products of the four ⊥ classes are collected the production
/// way (`per_class` on each), and the chart must break the lateral at the two **cut**
/// circles (z = 0, z = 20) as well as the rims: three intervals, not one full-height band
/// (`boundary_lines` and `chart_of`'s z-lines agree on the four). The middle interval is
/// both-cut and is emitted as panel rings; blinding **both** of the chart's axes (the arc
/// labels and the rulings') leaves the reader no answer there and `emit_lateral` refuses by
/// name (`CylinderGateUndecided`) — blinding only the arcs does not, because the vertical
/// lines answer in their place.
/// The through-boss's curved carriers and plane faces, collected the production way
/// (`per_class` on each of the four ⊥ classes and the wall class), with the cylinder's rows
/// and the four ⊥ classes `[z0, z20, cap_lo, cap_hi]`. ★ This mirrors production's fold by
/// hand (`trace_result_faces`' accumulation). A drift between them is not caught by
/// anything: a test would simply start measuring a map the boolean never builds.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn armed_curved(
    m: &Model,
    plate: Handle<Solid>,
    boss: Handle<Solid>,
    setup: &crate::arrangement::PlaneSetup,
    jd: &Judge<'_, WorkingPlane>,
    wc: usize,
    crossings: &std::collections::HashSet<(usize, usize)>,
    caps: [f64; 2],
) -> (
    Curved,
    Vec<LocalFace>,
    Vec<crate::bands::CylRow>,
    [usize; 4],
) {
    let z_class = |z: f64| -> usize {
        setup
            .geom
            .iter()
            .position(|p| {
                p.witness_coords()
                    .iter()
                    .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
            })
            .unwrap_or_else(|| panic!("a class at z = {z}"))
    };
    let (z0, z20, cap_lo, cap_hi) = (
        z_class(0.0),
        z_class(20.0),
        z_class(caps[0]),
        z_class(caps[1]),
    );
    // The production carriers: disk labels from every ⊥ class, cut rims from the cut ones.
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut plane_faces: Vec<LocalFace> = Vec::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    // ★ The wall class too: the chart's rulings are the wall's pieces, and without them the
    // both-cut middle is one whole-circle cell the emitter cannot read per sector.
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in [z0, z20, cap_lo, cap_hi, wc] {
        let edges = armed_class_edges(m, plate, boss, setup, jd, c, crossings);
        let staged = per_class(jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in &edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        plane_faces.extend(staged.faces);
    }
    let curved = Curved {
        aliases: Aliases::default(),
        disk_labels,
        arc_labels,
        cut_rims,
        rulings,
    };
    let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
    (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi])
}

#[test]
fn a_cut_circle_bounds_the_bands() {
    let (m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let (curved, plane_faces, rows, [z0, z20, cap_lo, cap_hi]) =
        armed_curved(&m, plate, boss, &setup, &jd, wc, &crossings, [-10.0, 40.0]);
    let (disk_labels, arc_labels, cut_rims) =
        (&curved.disk_labels, &curved.arc_labels, &curved.cut_rims);
    assert!(cut_rims.contains_key(&(0, z0)) && cut_rims.contains_key(&(0, z20)));
    assert!(disk_labels.contains_key(&(0, cap_lo)) && disk_labels.contains_key(&(0, cap_hi)));
    assert_eq!(arc_labels[&(0, z0)].len(), 2, "two arcs, two sector labels");
    assert_eq!(arc_labels[&(0, z20)].len(), 2);
    assert_eq!(rows.len(), 1);
    // ★ A cut circle is a band boundary (the rulings ladder): the chart's boundary rule names
    // the two cut rims beside the caps, and the chart's own lines are exactly those four.
    let def = &setup.cyls[0].def;
    let t_of = |c: usize| crate::bands::axis_param(&jd, c, def).unwrap();
    let expect: Vec<Rat> = [cap_lo, z0, z20, cap_hi].map(t_of).to_vec();
    assert_eq!(
        crate::arrangement::cyl_chart::boundary_lines(&jd, 0, def, &plane_faces, &curved, &rows)
            .unwrap(),
        expect,
        "three intervals, cut circles included"
    );
    let chart = crate::arrangement::cyl_chart::chart_of(&jd, &setup.cyls, 0, &plane_faces, &curved)
        .unwrap();
    assert_eq!(
        chart.z_lines.iter().map(|l| l.t).collect::<Vec<_>>(),
        expect,
        "the chart's lines are the boundaries and nothing else here"
    );
    // ★ **The emitter answers with regions.** Under `keep` the wall splits the middle
    // interval into two sectors, one kept and one not, and the kept one joins the whole
    // bands above and below it: **fuse** emits one face — a `Band` between the two cap
    // circles with the unkept sector as its one **hole** (a 4-node ring: two arcs on the
    // cut rims' own nodes, two rulings); **cut** keeps the other sector alone — one 4-node
    // `Ring` face. Which sector is the outer one is the assembly's and the volume oracle's to
    // measure (the gate-opening cell), not this harness's to re-derive.
    let faces_for = |kind: BoolKind| -> Vec<LocalFace> {
        crate::arrangement::cyl_chart::emit_lateral(
            kind,
            &jd,
            &setup.cyls,
            &plane_faces,
            &curved,
            &rows,
        )
        .unwrap()
    };
    let (fuse, cut) = (faces_for(BoolKind::Fuse), faces_for(BoolKind::Cut));
    assert_eq!(fuse.len(), 1, "fuse: one lateral face, a band with a hole");
    assert_eq!(cut.len(), 1, "cut: one lateral face, the kept sector");
    let rim = &cut_rims[&(0, z0)];
    // The ccw arc a ring carries on the lower cut rim, as an ordered node pair.
    let arc_on_z0 = |r: &crate::draft::Ring| -> [NodeId; 2] {
        let n = r.nodes.len();
        assert_eq!(n, 4, "two arcs and two rulings");
        let arcs = r
            .walls
            .iter()
            .filter(|w| matches!(w, crate::combinatorics::Wall::Arc { .. }))
            .count();
        assert_eq!(arcs, 2);
        for i in 0..n {
            let (a, b) = (r.nodes[i], r.nodes[(i + 1) % n]);
            if let crate::combinatorics::Wall::Arc { ccw, .. } = r.walls[i]
                && rim.nodes.contains(&a)
                && rim.nodes.contains(&b)
            {
                return if ccw { [a, b] } else { [b, a] };
            }
        }
        panic!("no arc on the lower cut rim");
    };
    let crate::draft::Bound::Band { lo, hi } = &fuse[0].outer else {
        panic!("fuse: a band between the caps, got {:?}", fuse[0].outer);
    };
    assert!(
        matches!(
            (lo, hi),
            (crate::draft::Rim::Circle(_), crate::draft::Rim::Circle(_))
        ),
        "the caps' whole circles are the band's rims"
    );
    assert_eq!(
        fuse[0].inner.len(),
        1,
        "the unkept sector is the band's one hole"
    );
    let crate::draft::Bound::Ring(hole) = &fuse[0].inner[0] else {
        panic!("a hole is a ring");
    };
    let crate::draft::Bound::Ring(panel) = &cut[0].outer else {
        panic!("cut: the kept sector is a ring, got {:?}", cut[0].outer);
    };
    assert!(cut[0].inner.is_empty());
    // ★ **Fuse's hole is cut's panel**: the sector fuse drops (inside the plate) is exactly
    // the sector cut keeps (the groove's wall), so the two rings carry the same ccw arc on
    // the lower cut rim — «complementary panels», stated for regions.
    let (pf, pc) = (arc_on_z0(hole), arc_on_z0(panel));
    assert_eq!(pf, pc, "fuse's hole is cut's panel");
    // Negative control: without the sector labels the both-cut interval's cells have no
    // speaking end, and the emitter refuses the class rather than guess a chamber. Called
    // directly, not through the census — this hand-broken input violates the very premise
    // (`src0_present == 0`) the census asserts where the fact is made.
    let mut blind = Curved {
        aliases: Aliases::default(),
        disk_labels: curved.disk_labels.clone(),
        arc_labels: HashMap::new(),
        cut_rims: curved.cut_rims.clone(),
        rulings: curved.rulings.clone(),
    };
    blind.arc_labels.clear();
    // ★ **Blinding one axis is not blinding the reader** — the chart has two, and the rulings
    // answer where the rims are silent (the vertical read's own cell). So this says the weaker,
    // truer thing: with the arcs gone the emitter still gets an answer, and
    // it is only when **both** axes are blinded that it refuses by name.
    assert!(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &plane_faces,
            &blind,
            &rows
        )
        .is_ok(),
        "the vertical lines answer what the blinded rims cannot"
    );
    for v in blind.rulings.values_mut() {
        for r in v.iter_mut() {
            r.label = None;
        }
    }
    assert!(matches!(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &plane_faces,
            &blind,
            &rows
        ),
        Err(BoolError::Rejected {
            reason: RejectReason::CylinderGateUndecided,
            ..
        })
    ));
}

/// The armed through-boss fuse **assembled to the end of the road** — the pipeline the
/// cell-3 lock walks, extracted so the consumer locks (props' θ-range integrals, tess's
/// open-rim merge) measure the very same solid through one spelling: every class's
/// arrangement, the coplanar unify, the band/panel pass, the seam table, `reconstruct`.
/// Returns the model with the one welded solid live, plus the setup and the wall class.
fn armed_assembled_through_boss() -> (Model, crate::arrangement::PlaneSetup, usize) {
    let (mut m, plate, boss, setup, wc, crossings) = armed_through_boss();
    let jd = Judge::new(&setup.geom, setup.standard, &setup.notes);
    let mut faces: Vec<LocalFace> = Vec::new();
    let mut disk_labels: crate::arrangement::DiskLabels = HashMap::new();
    let mut arc_labels: crate::arrangement::ArcLabels = HashMap::new();
    let mut cut_rims: CutRims = HashMap::new();
    let mut rulings: HashMap<usize, Vec<RulingExtent>> = HashMap::new();
    for c in 0..setup.geom.len() {
        let edges = armed_class_edges(&m, plate, boss, &setup, &jd, c, &crossings);
        let staged = per_class(&jd, &setup.cyls, BoolKind::Fuse, c, &edges).unwrap();
        for (cyl, label) in &staged.disk_labels {
            disk_labels.insert((*cyl, c), *label);
        }
        // ★ This mirrors production's fold by hand (`trace_result_faces`' accumulation). A
        // drift between them is not caught by anything: the test would simply start measuring
        // a map the boolean never builds.
        for (cyl, al) in &staged.arc_labels {
            arc_labels.entry((*cyl, c)).or_default().push(al.clone());
        }
        for (cyl, rim) in &edges.cut_rims {
            cut_rims.insert((*cyl, c), rim.clone());
        }
        for (cyl, r) in staged.ruling_extents {
            rulings.entry(cyl).or_default().push(r);
        }
        faces.extend(staged.faces);
    }
    let curved = Curved {
        aliases: Aliases::default(),
        disk_labels,
        arc_labels,
        cut_rims,
        rulings,
    };
    let faces = crate::assembly::unify_coplanar_faces(faces, &jd, &setup.cyls).unwrap();
    let rows = crate::bands::cyl_rows(&setup.planes, &setup.plane_ix, setup.n_a).unwrap();
    let mut faces = faces;
    faces.extend(
        crate::arrangement::cyl_chart::emit_lateral(
            BoolKind::Fuse,
            &jd,
            &setup.cyls,
            &faces,
            &curved,
            &rows,
        )
        .unwrap(),
    );
    let seam = seam_table(&m, &faces, &setup.cyls, &jd).unwrap();
    let out = crate::assembly::reconstruct(
        &mut m,
        &jd,
        &seam,
        &faces,
        &setup.cyls,
        &curved.cut_rims,
        None,
        crate::assembly::Tangencies::none(),
    )
    .unwrap();
    assert_eq!(out.len(), 1, "one welded solid");
    m.restore_live(out);
    m.rebuild_adjacency();
    (m, setup, wc)
}

/// ★ **The armed assembly welds the panels** (the lock): the whole
/// through-boss fuse, from production parts past the standing gate — every class's
/// arrangement, the coplanar unify, the band/panel pass, the seam table, and
/// `reconstruct` — comes back one solid with a clean `validate`. The ruling edges are
/// pinned structurally: exactly two straight lateral edges (the kept outer panel's), each
/// used exactly twice, carriers stated as **the edge's own fact** — the cylinder and the
/// wall plane — and the outer panel's seam-holding arc is split at an `OnSeam` vertex
/// (the panel road runs the wrap-arc split the Band road already had).
#[test]
fn the_armed_assembly_welds_the_panels() {
    let (m, setup, wc) = armed_assembled_through_boss();
    let issues = nacre_validate::validate(&m);
    assert!(issues.is_empty(), "{issues:?}");
    // The ruling edges, structurally: straight lateral edges = carriers {cylinder, a plane}
    // with a Line curve (the seam's [cyl, cyl] spelling is excluded by the mixed pair).
    let cyl_surf = setup.cyls[0].surf;
    let wall_surf = setup.geom[wc].surf;
    let reach = m.reachable();
    let mut rulings = 0;
    let mut i = 0u32;
    while let Some(eh) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(eh);
        if !reach.edges.contains(&eh) {
            continue;
        }
        let mixed = (e.surfaces[0] == cyl_surf) != (e.surfaces[1] == cyl_surf);
        if !mixed {
            continue;
        }
        let curve = m.derive_edge_curve(e.surfaces, e.vertices).unwrap();
        if matches!(curve, nacre_geom::Curve::Line(_)) {
            rulings += 1;
            assert!(
                e.surfaces.contains(&wall_surf),
                "a ruling's plane carrier is the wall: {:?}",
                e.surfaces
            );
        }
    }
    assert_eq!(rulings, 2, "the kept outer panel's two rulings");
    // The seam split ran on the panel: an OnSeam vertex is reachable.
    let on_seam = reach
        .vertices
        .iter()
        .filter(|&&v| matches!(*m.vertex(v), nacre_topo::Vertex::OnSeam(_)))
        .count();
    assert!(on_seam >= 1, "the outer panel's wrap arc split at the seam");
}

/// ★ **The armed solid's mass properties are exact** (the props
/// instrument): `mass_props` on the assembled through-boss fuse, with the θ-range lateral
/// integrals live, answers the derived closed forms — volume `32000 + 1000π` (plate plus
/// the boss outside it), area `6200 + 425π` (walls 3000 + split x = 40 wall 600, caps
/// bitten `3200 − 25π`, boss disks `50π`, full bands `300π`, the outer half-panel `100π`).
/// The lifted probe's full-2π mis-answer (7849.34 = truth + the panel counted whole,
/// `+100π`) cross-checks the area derivation. This is the volume oracle standing while the
/// gate still refuses production input.
#[test]
fn the_armed_solids_mass_properties_are_exact() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let props = nacre_props::mass_props(&m, m.live_solids()[0]).unwrap();
    let volume = 32000.0 + 1000.0 * std::f64::consts::PI;
    let area = 6200.0 + 425.0 * std::f64::consts::PI;
    assert!(
        (props.volume - volume).abs() <= 1e-9 * volume,
        "volume {} != {volume}",
        props.volume
    );
    assert!(
        (props.area - area).abs() <= 1e-9 * area,
        "area {} != {area}",
        props.area
    );
}

/// ★ **The armed solid tessellates watertight** (the tess
/// instrument): the open-rim merge triangulates the θ-panels of the assembled through-boss
/// fuse — the outer panel's rims are two arcs split at the seam vertex, so this one solid
/// exercises the multi-polyline open chain *and* the seam-crossing unwrap — and every
/// undirected triangle edge is used exactly twice.
#[test]
fn the_armed_solid_tessellates_watertight() {
    let (m, _setup, _wc) = armed_assembled_through_boss();
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).unwrap();
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    let open = uses.values().filter(|&&n| n != 2).count();
    assert_eq!(open, 0, "the mesh is watertight");
}

/// **The mixed parity reads a bitten ring — exactly, on the production pieces.**
///
/// The straddling boss cuts the plate-top ring at `(40, 15)` and `(40, 25)`, so that ring
/// carries two pierce corners and one arc step; the overhang digon is chord + outer arc.
/// Probes are derived, not read back: the chart for `n = +z` picks `e1 = [0, −1, 0]`, so
/// the ray runs toward −y at fixed x. `[37, 30]` is inside and its ray crosses the **arc
/// twice** (`(37−40)² + (y−20)² = 25` → y = 16, 24) before the bottom edge — the arc arm is
/// what that probe measures, and a chord-minded parity would answer it wrong. `[37, 18]`
/// sits inside the bite (outside the face), `[50, 20]` outside everything; `[43, 20]` is
/// inside the overhang (one arc crossing at y = 16). Cells are identified structurally,
/// not by size: the outside cell (winding −1) owns the outer arc piece's twin, so that
/// piece's forward cell is the overhang digon; the other digon is the bite.
#[test]
fn the_mixed_parity_reads_a_bitten_ring() {
    with_bitten_rings(|jd, cyls, wc, big, overhang, bite| {
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
        let rat = |x: i128, y: i128| {
            [
                nacre_exact::Rat::from_int(x),
                nacre_exact::Rat::from_int(y),
                nacre_exact::Rat::from_int(20),
            ]
        };
        let ask = |ring: &[combinatorics::RingEdge], p: [nacre_exact::Rat; 3]| {
            combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring)
        };
        assert_eq!(ask(big, rat(12, 12)), Some(true), "plain interior");
        assert_eq!(
            ask(big, rat(37, 30)),
            Some(true),
            "interior whose ray crosses the inner arc twice"
        );
        assert_eq!(ask(big, rat(37, 18)), Some(false), "inside the bite");
        assert_eq!(ask(big, rat(50, 20)), Some(false), "outside everything");
        assert_eq!(
            ask(overhang, rat(43, 20)),
            Some(true),
            "inside the overhang"
        );
        assert_eq!(
            ask(overhang, rat(37, 18)),
            Some(false),
            "the bite is not the overhang"
        );
        assert_eq!(ask(bite, rat(37, 18)), Some(true), "inside the bite digon");
        assert_eq!(
            ask(bite, rat(43, 20)),
            Some(false),
            "the overhang is not the bite"
        );
    });
}

/// ★★★★★ **The digon's truth is known independently, so the parity can be swept rather than
/// spot-checked.** A digon of a chord and an arc is `disk ∩ half-space`, and both halves are
/// exact rational predicates the arrangement already owns —
/// [`nacre_exact::quad::cylinder_radial_side`] and the sign of the wall's plane equation. A
/// grid over the circle's neighbourhood therefore checks **every** answer, and it is the only
/// control here that crosses the straight arm, the arc arm **and their junction**: the
/// non-mixed road's oracle cannot see an arc, and a whole circle's arm collapses to
/// `Ordering::Equal` on a single step.
///
/// ☑ **What it does and does not reach, measured.** On this fixture the chord's ends are the
/// seam and the tangent columns, and a ray along the chord meets **both** arc ends at once —
/// so this sweep is green under a flipped arc-end departure sign (two ends on one ray flip
/// together and keep the parity) and green under the straight arm's old corner abstention
/// (negative controls, both measured). It is the oracle for the mixed road as a
/// whole — the straight arm, the arc arm and their junction on ordinary roots — and the
/// corner-bitten plate (`a_root_at_an_arc_end_is_a_corner_on_the_ray`) is the watch on the
/// arc-end arm, where a ray meets one end alone.
#[test]
fn the_mixed_parity_agrees_with_the_digon_it_bounds() {
    use nacre_exact::{Orient, Rat};
    with_bitten_rings(|jd, cyls, wc, _big, overhang, bite| {
        let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
        let def = cyls
            .iter()
            .map(|c| &c.def)
            .find(|d| *d.r2() == nacre_exact::BigRat::from(Rat::from_int(25)))
            .expect("the bitten circle");
        let (mut swept, mut on_boundary) = (0usize, 0usize);
        let (mut abstained, mut inside_seen) = (0usize, 0usize);
        let mut on_chord = 0usize;
        for i in 60..=100 {
            for j in 20..=60 {
                let p = [
                    Rat::new(i.into(), 2).unwrap(),
                    Rat::new(j.into(), 2).unwrap(),
                    Rat::from_int(20),
                ];
                let radial = nacre_exact::quad::cylinder_radial_side(
                    &p,
                    &def.origin(),
                    &def.dir(),
                    def.r2(),
                );
                // The chord is the plate's wall `x = 40`; the bite keeps `x < 40`.
                let wall = p[0].checked_sub(Rat::from_int(40)).unwrap();
                if radial == Orient::Zero || wall == Rat::from_int(0) {
                    on_boundary += 1;
                    // ★★★★★ **A point strictly inside the circle and *on* the chord is on
                    // both digons' boundary, and "inside" has no answer there.** These are the
                    // grid's sharpest points: the ray along the chart's first axis leaves such
                    // a point straight through **both arc endpoints** — the chord is a step
                    // lying along the ray, and the straight arm names the probe between its
                    // ends as the boundary before any arc is asked.
                    if radial == Orient::Negative && wall == Rat::from_int(0) {
                        for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                            assert_eq!(
                                combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                None,
                                "{who} must abstain on its own chord at {p:?}"
                            );
                        }
                        on_chord += 1;
                    }
                    continue;
                }
                let inside_bite = radial == Orient::Negative && wall < Rat::from_int(0);
                let inside_over = radial == Orient::Negative && wall > Rat::from_int(0);
                for (ring, want, who) in [
                    (bite, inside_bite, "bite"),
                    (overhang, inside_over, "overhang"),
                ] {
                    // ★ **An abstention is allowed and a wrong answer is not.** The ray runs
                    // along the chart's first axis, so a grid point whose ray meets a ring
                    // corner has no parity to report — the caller's remedy is another point,
                    // which is exactly what `coord_probes` does with its candidate list. What
                    // the sweep locks is that every answer it *does* give is the truth.
                    match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                        Some(got) => {
                            assert_eq!(got, want, "{who} at {p:?}");
                            swept += 1;
                            inside_seen += usize::from(want);
                        }
                        None => abstained += 1,
                    }
                }
            }
        }
        // Neither vacuous nor all-outside: the grid straddles the circle and the chord, and
        // the two digons between them own a real interior.
        assert!(swept > 1_000, "answers swept: {swept}");
        assert!(inside_seen > 100, "points inside a digon: {inside_seen}");
        assert!(
            on_boundary > 0,
            "the grid meets the boundary: {on_boundary}"
        );
        // Recorded rather than bounded: every abstention here is the tangent ray (the rows
        // `y = 15, 25` graze the circle — measured: 160 of 160, no corner among
        // them), and its count is a property of this grid, not of the rule.
        assert!(abstained > 0, "abstentions: {abstained}");
        // What those abstentions were, by kind.
        {
            let rows = combinatorics::tie_probe::ROWS.mine();
            let mut hist: Vec<(combinatorics::tie_probe::Tie, usize)> = Vec::new();
            for t in &rows {
                match hist.iter_mut().find(|(k, _)| k == t) {
                    Some((_, c)) => *c += 1,
                    None => hist.push((*t, 1)),
                }
            }
            eprintln!("tie digon abstained {abstained} on_chord {on_chord} kinds {hist:?}");
        }
        assert!(
            on_chord > 0,
            "points on the chord inside the circle: {on_chord}"
        );
    });
}

/// The bitten-plate fixture, handed to `f` as `(judge, cyls, class, big, overhang, bite)`.
/// The bitten fixture: one wall (`x = 40`) cuts the circle into two digons.
fn with_bitten_rings(
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    with_cut_rings([40.0, 40.0, 20.0], false, f);
}

/// A plate `[0, plate]` with a z-cylinder (r = 5, h = 10) standing on its top at
/// `(40, 20)`: the class carrying the cut circle, its judge, and the three bounded rings —
/// the bitten top, the overhang (disk ∖ plate) and the bite (disk ∩ plate). The plate's
/// extent picks the cut: `[40, 40, 20]` bites with one wall (two digons); `[40, 24, 20]`
/// with two — the corner `(40, 24)` sits inside the circle and both rings are trigons
/// whose arc ends at the rational `(37, 24)`, off the seam and off the tangent columns.
/// `seam_off` puts the seam on `+x`, at `(45, 20)` — on no ring corner — instead of on `−y`,
/// which puts it on the chord end `(40, 15)`.
fn with_cut_rings(
    plate: [f64; 3],
    seam_off: bool,
    f: impl FnOnce(
        &Judge<'_, WorkingPlane>,
        &[crate::planes::WorkingCyl],
        usize,
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
        &[combinatorics::RingEdge],
    ),
) {
    let mut m = Model::new();
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array(plate));
    // The seam off every ring corner (`+x`, at `(45, 20)`) or on the chord end (`−y`, `(40, 15)`).
    let seam = if seam_off {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, -1.0, 0.0]
    };
    let b = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([40.0, 20.0, 20.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array(seam),
        5.0,
        10.0,
    )
    .solid;
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    // The one class whose circle is cut: find it by running the split everywhere.
    let mut found = None;
    for wc in 0..planes.len() {
        if !matches!(plane_ix.get(wc), Some(ClassIx::Plane(_)) | None) && wc < plane_ix.len() {
            continue;
        }
        if combinatorics::class_coeffs_rat(&jd, wc).is_none() {
            continue;
        }
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()) else {
            continue;
        };
        let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
            continue;
        };
        let Ok(edges) = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
        else {
            continue;
        };
        if edges.arcs.is_empty() {
            continue;
        }
        let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
        let ns = edges.segs.len();
        let na = edges.arcs.len();
        // Four cells share this class: the bitten plate-top face (+1, many steps), the
        // outside (−1 — it also borders the arc, via the outer bulge), and two digons
        // (chord + arc): the overhang (disk minus plate) and the bite (disk ∩ plate).
        // Tell them apart structurally: the outside owns the outer piece's twin, so the
        // outer piece's forward cell is the overhang; the remaining digon is the bite.
        let outside = cells
            .iter()
            .position(|c| c.winding == -1)
            .expect("the outside cell");
        let bounded: Vec<usize> = (0..cells.len())
            .filter(|&i| cells[i].winding == 1)
            .collect();
        assert_eq!(
            bounded.len(),
            3,
            "the bitten top, the overhang and the bite"
        );
        let big_ix = *bounded
            .iter()
            .max_by_key(|&&i| cells[i].half_edges.len())
            .expect("the bitten top ring");
        let outer_ai = (0..na)
            .find(|ai| face_of[&(2 * (ns + ai) + 1)] == outside)
            .expect("the outer piece borders the outside");
        let overhang_ix = face_of[&(2 * (ns + outer_ai))];
        let bite_ix = *bounded
            .iter()
            .find(|&&i| i != big_ix && i != overhang_ix)
            .expect("the bite");
        assert_eq!(
            cells[overhang_ix].half_edges.len(),
            cells[bite_ix].half_edges.len(),
            "the same walls cut the overhang and the bite"
        );
        let ring = |ix: usize| -> Vec<combinatorics::RingEdge> {
            cells[ix]
                .half_edges
                .iter()
                .map(|&he| edges.edge_at(he))
                .collect()
        };
        found = Some((wc, ring(big_ix), ring(overhang_ix), ring(bite_ix)));
        break;
    }
    let (wc, big, overhang, bite) = found.expect("one class carries the cut circle");
    f(&jd, &cyls, wc, &big, &overhang, &bite);
}

/// **A root at an arc's own end is a corner on the ray, and the arc's tangent there says
/// whether the arc counts it** — the half-open rule in the arc arm, on the one
/// fixture that reaches it.
///
/// The bitten fixture cannot: its chord's ends are the seam and the tangent columns, and
/// each abstains first (measured over the whole suite: the arc arm's end rule decides nothing
/// there, while the straight arm's corner tie speaks for the same corners). The plate's second wall
/// `y = 24` puts the corner `(40, 24)` inside the circle, so the bite's arc ends at
/// `(37, 24)` — a 3-4-5 point, rational, off the seam `(40, 15)` and off the tangent
/// columns `x = 35, 45` — and the lattice column `x = 37` sends its rays through that end
/// **alone**, beside an ordinary second root at `(37, 16)`. That is the single-end crossing
/// a global sign error cannot hide: two ends on one ray flip together and keep the parity
/// (the bitten fixture's chord), one end on a ray does not. Every answer the sweep gives
/// is the truth, the boundary abstains, and the arm decided the end exactly once per point
/// on the shooting side of that column, for both rings — under both seam placements, so
/// the seam-incident arms and the `(false, false)` arm each decide an end.
#[test]
fn a_root_at_an_arc_end_is_a_corner_on_the_ray() {
    use nacre_exact::{Orient, Rat};
    for seam_off in [false, true] {
        with_cut_rings(
            [40.0, 24.0, 20.0],
            seam_off,
            |jd, cyls, wc, _big, overhang, bite| {
                let coeffs = combinatorics::class_coeffs_rat(jd, wc).unwrap();
                let def = cyls
                    .iter()
                    .map(|c| &c.def)
                    .find(|d| *d.r2() == nacre_exact::BigRat::from(Rat::from_int(25)))
                    .expect("the cut circle");
                let decided0 = combinatorics::tie_probe::arc_end_decisions_mine();
                let (mut swept, mut abstained, mut inside_seen, mut boundary) = (0usize, 0, 0, 0);
                let mut column = [0usize; 2];
                let rat = |k: i128| Rat::new(k, 2).unwrap();
                for i in 60..=100 {
                    for j in 20..=60 {
                        let p = [rat(i), rat(j), Rat::from_int(20)];
                        let radial = nacre_exact::quad::cylinder_radial_side(
                            &p,
                            &def.origin(),
                            &def.dir(),
                            def.r2(),
                        );
                        let (x, y) = (p[0], p[1]);
                        let (wx, wy) = (x == Rat::from_int(40), y == Rat::from_int(24));
                        // Both rings' boundary: the circle and the two chords, `x = 40` for
                        // `15 ≤ y ≤ 24` and `y = 24` for `37 ≤ x ≤ 40`, ends included.
                        let on_chord = (wx && y >= Rat::from_int(15) && y <= Rat::from_int(24))
                            || (wy && x >= Rat::from_int(37) && x <= Rat::from_int(40));
                        if radial == Orient::Zero || on_chord {
                            boundary += 1;
                            for (ring, who) in [(bite, "bite"), (overhang, "overhang")] {
                                assert_eq!(
                                    combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring),
                                    None,
                                    "{who} must abstain on its boundary at {p:?}"
                                );
                            }
                            continue;
                        }
                        let in_disk = radial == Orient::Negative;
                        let on_plate = x < Rat::from_int(40) && y < Rat::from_int(24);
                        for (k, (ring, want, who)) in [
                            (bite, in_disk && on_plate, "bite"),
                            (overhang, in_disk && !on_plate, "overhang"),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            match combinatorics::point_in_mixed_ring(jd, cyls, &coeffs, &p, ring) {
                                Some(got) => {
                                    assert_eq!(got, want, "{who} at {p:?}");
                                    swept += 1;
                                    inside_seen += usize::from(want);
                                    if x == Rat::from_int(37) && radial == Orient::Positive {
                                        column[k] += 1;
                                    }
                                }
                                None => abstained += 1,
                            }
                        }
                    }
                }
                assert!(swept > 1_000, "answers swept: {swept}");
                assert!(inside_seen > 100, "points inside a ring: {inside_seen}");
                assert!(boundary > 0, "the grid meets the boundary: {boundary}");
                // The `x = 37` column outside the circle: 24 lattice points, every one answered by
                // both rings — the shooting side through the end, the other side by a miss — and
                // the end decided exactly once per point on the shooting side: 12 × 2 rings.
                assert_eq!(column, [24, 24], "the single-end column is answered");
                let decided = combinatorics::tie_probe::arc_end_decisions_mine() - decided0;
                eprintln!(
                    "arc-end sweep (seam_off {seam_off}): swept {swept} abstained {abstained} \
                 boundary {boundary} arc-end decisions {decided}"
                );
                // With the seam on the chord end `(40, 15)` (`seam_off = false`) the
                // `(true, false)` / `(false, true)` arms decide the `x = 37` column's shooting side:
                // 12 points × 2 rings. With the seam off every corner (`ref_dir = x`) both ends are
                // ordinary and the `(false, false)` arm decides — and the `x = 40` column's rays now
                // meet `(40, 15)` alone as well (its other root `(40, 25)` is the overhang's
                // interior), 12 more per ring, one of them a call that then abstains at the corner
                // `(40, 24)`; the column's seam-root abstentions are gone (182 → 160).
                assert_eq!(
                    decided,
                    if seam_off { 48 } else { 24 },
                    "the arc-end arm decided the single-end rays (seam_off {seam_off})"
                );
            },
        );
    }
}

/// **On a ring whose every corner is rational, the mixed road and the chart road are one
/// rule**: `point_in_ring_2d_rat` says Inside/Outside/OnBoundary and
/// `point_in_mixed_ring` `Some(true)`/`Some(false)`/`None`, and the pairs match at every
/// lattice point of every bounded ring of every class of two overlapping plates — squares
/// and L-shapes with reflex corners. A lattice at half steps shares a row with every corner
/// and every horizontal step, so every corner-on-the-ray configuration the half-open rule
/// distinguishes (both steps up, both down, one each, a step along the ray, the probe at
/// the corner) is on the sweep. A mixed road that abstained at such a corner would do so where
/// the chart road answers (measured over the whole suite's rational rings: 5 such corners, 0
/// disagreements); a whole-suite shadow of the two roads costs 61 % of the serial sweep and is
/// not kept — this lattice is the lock.
#[test]
fn the_two_roads_agree_on_every_rational_ring() {
    use nacre_exact::Rat;
    use nacre_geom::intersect::{RingSide, point_in_ring_2d_rat};
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([4.0, 4.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([2.0, 2.0, 0.0]),
        Point3::from_array([6.0, 6.0, 2.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let (mut rings_swept, mut points, mut inside, mut boundary, mut reflex) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    for wc in 0..planes.len() {
        let Some(coeffs) = combinatorics::class_coeffs_rat(&jd, wc) else {
            continue;
        };
        let tr = trace_on_class_of(
            &m,
            a,
            b,
            wc,
            &jd,
            &[],
            &faces_tab,
            &surf_ix,
            &inc_a,
            &inc_b,
            &plane_ix,
            Default::default(),
        );
        let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
        let Ok(split) = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()) else {
            continue;
        };
        let Ok(circles) = merge_circles(&tr.circles, &cyls, &Aliases::default()) else {
            continue;
        };
        let Ok(edges) = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default())
        else {
            continue;
        };
        let Ok((cells, _)) = walk_cells(&jd, &cyls, wc, &edges) else {
            continue;
        };
        let n = [coeffs[0], coeffs[1], coeffs[2]];
        let chart = combinatorics::Chart2dRat::of_normal(&n).unwrap();
        let (e1, e2) = chart.axes();
        // Axis-aligned classes: unit chart axes and a unit normal, so a chart point
        // `(X, Y)` lifts to `X·e1 + Y·e2 − d·n`.
        for e in [e1, e2, &n] {
            assert_eq!(nacre_exact::dot3_rat(e, e), Some(Rat::from_int(1)));
        }
        let lift = |xy: [Rat; 2]| -> [Rat; 3] {
            let mut p = [Rat::from_int(0); 3];
            for k in 0..3 {
                p[k] = xy[0]
                    .checked_mul(e1[k])
                    .unwrap()
                    .checked_add(xy[1].checked_mul(e2[k]).unwrap())
                    .unwrap()
                    .checked_sub(coeffs[3].checked_mul(n[k]).unwrap())
                    .unwrap();
            }
            p
        };
        for c in cells.iter().filter(|c| c.winding == 1) {
            let ring: Vec<combinatorics::RingEdge> =
                c.half_edges.iter().map(|&he| edges.edge_at(he)).collect();
            assert!(!combinatorics::ring_is_mixed(&ring), "a planar fixture");
            let nodes: Vec<NodeId> = ring.iter().map(|e| e.node).collect();
            let ring2 = chart.ring(&jd, &nodes).expect("rational corners");
            if ring2.len() > 4 {
                reflex += 1;
            }
            let bound = |k: usize, lo: bool| -> i128 {
                let it = ring2.iter().map(|q| {
                    assert_eq!(q[k].denom(), 1, "integer corners");
                    q[k].numer()
                });
                if lo {
                    it.min().unwrap() - 1
                } else {
                    it.max().unwrap() + 1
                }
            };
            for xi in 2 * bound(0, true)..=2 * bound(0, false) {
                for yi in 2 * bound(1, true)..=2 * bound(1, false) {
                    let q = [Rat::new(xi, 2).unwrap(), Rat::new(yi, 2).unwrap()];
                    let p = lift(q);
                    assert_eq!(
                        chart.project(&p),
                        Some(q),
                        "the lift is the chart's inverse"
                    );
                    let plane = point_in_ring_2d_rat(q, &ring2);
                    let mixed = combinatorics::point_in_mixed_ring(&jd, &cyls, &coeffs, &p, &ring);
                    match (plane, mixed) {
                        (RingSide::Inside, Some(true)) => inside += 1,
                        (RingSide::Outside, Some(false)) => {}
                        (RingSide::OnBoundary, None) => boundary += 1,
                        other => panic!("class {wc} ring {ring2:?} at {q:?}: {other:?}"),
                    }
                    points += 1;
                }
            }
            rings_swept += 1;
        }
    }
    assert!(rings_swept >= 3, "rings swept: {rings_swept}");
    assert!(reflex > 0, "an L-shaped ring: {reflex}");
    assert!(
        inside > 0 && boundary > 0,
        "points {points} inside {inside} boundary {boundary}"
    );
    eprintln!("two roads: rings {rings_swept} points {points} inside {inside} boundary {boundary}");
}

/// The gated drill population's fixture: a `[0,2]³` box and an axis-aligned cylinder at
/// `(1,1)`, r=0.5 — every wall is a full unit from the axis, so the population gate passes
/// and [`plane_index_setup`] hands the arrangement bricks a cylinder-bearing table.
fn drilled(m: &mut Model, z0: f64, h: f64) -> (Handle<Solid>, Handle<Solid>) {
    let a = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([2.0; 3]));
    let b = crate::fixtures::cylinder_with_seam(
        m,
        Point3::from_array([1.0, 1.0, z0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        0.5,
        h,
    )
    .solid;
    m.rebuild_adjacency();
    (a, b)
}

/// The plane class whose defining triangle lies wholly at `z` — a z-cap.
fn class_at_z(planes: &[WorkingPlane], z: f64) -> usize {
    planes
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[2] - z).abs() < 1e-12)
        })
        .expect("a z-cap class")
}

/// A cylinder cap alone on its class: the disk's circular outer traces as a **seated
/// circle** — and the lateral, whose rim lies on this very plane, adds its **graze** beside it.
/// The cell bricks turn the pair into a disk `+1` / contour `−1` pair whose labels say "body on
/// the cap's inside", with no segments anywhere. The emitted face's outer is the circle itself
/// ([`crate::draft::Bound::Circle`]), the vocabulary the assembly consumes.
#[test]
fn a_cylinder_cap_is_a_seated_circle_and_a_disk_face() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // caps at z=-1 and z=3, clear of the box
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 3.0);

    let mut tr = Trace::default();
    trace_one_of(
        &m,
        b,
        SolidSide::B,
        wc,
        &jd,
        &cyls,
        &faces_tab,
        &surf_ix,
        &inc_b,
        &plane_ix,
        Default::default(),
        &mut tr,
    );
    assert!(
        tr.segs.is_empty(),
        "a circle owes the segment machinery nothing"
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    // The class's normal is the cap's own outward +z (the cap is its only member), so the
    // body — below z=3 — is `body_above: false`. ★ **Two contributions, not one**: the cap is
    // seated here and the lateral's rim grazes here, and on this convex cap they agree.
    assert!(planes[wc].plane.normal().as_array()[2] > 0.0);
    let mut kinds: Vec<SegKind> = tr
        .circles
        .iter()
        .inspect(|c| assert!(c.cyl == 0 && c.solid == SolidSide::B, "{:?}", tr.circles))
        .map(|c| c.kind)
        .collect();
    kinds.sort_by_key(|k| matches!(k, SegKind::Graze { .. }));
    assert!(
        matches!(
            kinds[..],
            [
                SegKind::Seated { body_above: false },
                SegKind::Graze { body_above: false },
            ]
        ),
        "{:?}",
        tr.circles
    );

    let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
    let edges = ClassEdges::of(&jd, &cyls, wc, &[], &circles, &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
    // Pseudo-half-edges 0 and 1 (no segments): the disk (+1) and its contour (−1).
    assert_eq!(cells.len(), 2, "{cells:?}");
    assert_eq!(
        (cells[0].half_edges.as_slice(), cells[0].winding),
        (&[0][..], 1)
    );
    assert_eq!(
        (cells[1].half_edges.as_slice(), cells[1].winding),
        (&[1][..], -1)
    );
    let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.root_groups,
        vec![1],
        "the contour bounds the unbounded region"
    );
    assert!(nesting.holes.is_empty());
    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    assert_eq!(labels[1], [false; 4]);
    assert_eq!(
        labels[0],
        [false, false, false, true],
        "inside the circle, B below the plane only"
    );
    // Fuse keeps below and not above across the disk → the disk is a result face, and its
    // outer boundary is the circle — no ring, no nodes.
    let (out, _, _) = emit_faces(
        BoolKind::Fuse,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "{out:?}");
    assert!(matches!(
        out[0].outer,
        crate::draft::Bound::Circle { cyl: 0 }
    ));
    assert!(out[0].inner.is_empty());
}

/// The through-drill's cap arrangement, brick by brick: the box cap is a seated 4-ring, the
/// lateral crosses the cap plane in a **transversal circle**, and the bricks nest the circle
/// as the cap cell's hole, label the disk with the cylinder straddling `W`, and emit — for
/// `Cut` — exactly one face: the cap with a circular hole. The per-kind `edge_mask`
/// difference is what the two labels measure (seated flips one side, transversal flips
/// both).
#[test]
fn a_drill_circle_is_a_hole_of_the_cap_ring() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0); // z∈[-1,3]: through both caps of the box
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        cyls,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 0.0);

    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::B,
                kind: SegKind::Transversal { .. },
                // No hole in this band, so the mark is the whole circle.
                arc: None,
            }]
        ),
        "{:?}",
        tr.circles
    );

    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    let circles = merge_circles(&tr.circles, &cyls, &Aliases::default()).unwrap();
    let split = split_at_crossings(&jd, &cyls, wc, &merged, &mut Aliases::default()).unwrap();
    assert_eq!(
        split.len(),
        4,
        "the cap ring alone — the circle owes the splitter nothing"
    );

    let edges = ClassEdges::of(&jd, &cyls, wc, &split, &circles, &[], &Aliases::default()).unwrap();
    let (cells, face_of) = walk_cells(&jd, &cyls, wc, &edges).unwrap();
    assert_eq!(cells.len(), 4, "cap ±1 and circle ±1: {cells:?}");
    let at = |he: usize| cells.iter().position(|c| c.half_edges == [he]).unwrap();
    let (disk, contour) = (at(2 * split.len()), at(2 * split.len() + 1));
    let cap = cells
        .iter()
        .position(|c| c.winding == 1 && c.half_edges.len() == 4)
        .unwrap();

    let nesting = nest_cells(&jd, &cyls, wc, &cells, &edges).unwrap();
    assert_eq!(
        nesting.holes.get(&cap).map(Vec::as_slice),
        Some(&[contour][..]),
        "the circle contour is the cap cell's hole"
    );

    let labels = label_cells(&cells, &face_of, &edges, &nesting, [false; 4]).unwrap();
    // The class is the box's **bottom** cap, so which of `W`'s two sides carries the box is
    // read off the class normal rather than assumed: the box occupies z>0.
    let up = planes[wc].plane.normal().as_array()[2] > 0.0;
    let box_side = usize::from(!up); // 0 = above `n_W`, 1 = below
    let mut cap_label = [false; 4];
    cap_label[box_side] = true;
    assert_eq!(
        labels[cap], cap_label,
        "outside the circle: box material on the z>0 side only"
    );
    assert_eq!(
        labels[contour], labels[cap],
        "a hole bounds its host's region"
    );
    let mut disk_label = cap_label;
    disk_label[2] = true;
    disk_label[3] = true;
    assert_eq!(
        labels[disk], disk_label,
        "inside the circle the cylinder straddles W, the box side is unchanged"
    );

    let (out, _, _) = emit_faces(
        BoolKind::Cut,
        &labels,
        &cells,
        &edges,
        &jd,
        wc,
        &nesting.holes,
    );
    assert_eq!(out.len(), 1, "one face: the drilled cap — {out:?}");
    assert_eq!(out[0].outer.expect_ring().len(), 4);
    assert!(
        matches!(out[0].inner[..], [crate::draft::Bound::Circle { cyl: 0 }]),
        "{:?}",
        out[0].inner
    );
}

/// The hole arm of the seated tracer: a face whose inner loop is a circle (a two-hole
/// plate's input shape — a bore `Cut` builds one in production now, so the hand-doctored
/// loops this arm was written against have company) emits its polygon segments **and** a
/// seated circle that inherits the face's own
/// body side.
#[test]
fn a_circular_hole_ring_traces_as_a_seated_circle() {
    let mut m = Model::new();
    let (a, b) = drilled(&mut m, -1.0, 4.0);
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = class_at_z(&planes, 0.0);

    let input = combinatorics::trace_input(
        &m,
        [(a, &inc_a), (b, &inc_b)],
        &surf_ix,
        faces_tab.len(),
        &jd,
        &plane_ix,
        &[],
        Default::default(),
    );
    let (fp, fl) = input.faces[0]
        .iter()
        .find(|(fp, _)| matches!(plane_ix[*fp], ClassIx::Plane(c) if c == wc))
        .expect("the box's bottom cap sits on wc");
    let doctored = combinatorics::FaceLoops {
        outer: fl.outer.clone(),
        holes: Some(vec![combinatorics::LoopRing::Circle { cyl: 0 }]),
        cycles: None,
    };
    let mut tr = Trace::default();
    trace_one(
        &[(*fp, doctored)],
        SolidSide::A,
        wc,
        &jd,
        &[],
        &faces_tab,
        &plane_ix,
        &Default::default(),
        &Aliases::default(),
        &mut tr,
    );
    assert!(tr.declined.is_empty(), "{:?}", tr.declined);
    assert_eq!(tr.segs.len(), 4, "the polygon outer still emits its edges");
    let SegKind::Seated { body_above } = tr.segs[0].kind else {
        panic!("a face on wc traces seated: {:?}", tr.segs[0]);
    };
    assert!(
        matches!(
            tr.circles[..],
            [CircleTrace {
                cyl: 0,
                solid: SolidSide::A,
                kind: SegKind::Seated { body_above: ba },
                arc: None
            }]
            if ba == body_above
        ),
        "the hole circle inherits the face's body side: {:?}",
        tr.circles
    );
}

/// The existence condition, negatively: a cap plane **outside** the lateral's rim span gets
/// no circle (and no decline — a miss, like a parallel plane), while the rim-interior cap
/// still does, and a rim-**coincident** plane carries the cap's seated circle rather than a
/// transversal one. A ghost circle here would corrupt every label on the class.
#[test]
fn no_ghost_circle_outside_the_rim_span() {
    let mut m = Model::new();
    // z∈[-1,1.5]: through the box's bottom cap, short of its top cap at z=2.
    let (a, b) = drilled(&mut m, -1.0, 2.5);
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);

    // The box's top cap at z=2 is past the rim span [−1, 1.5]: no circle, no decline.
    let top = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 2.0),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert!(top.circles.is_empty(), "{:?}", top.circles);
    assert!(top.declined.is_empty(), "{:?}", top.declined);

    // The bottom cap at z=0 is strictly inside the span: the transversal circle is there.
    let bottom = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 0.0),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert_eq!(
        bottom
            .circles
            .iter()
            .filter(|c| matches!(c.kind, SegKind::Transversal { .. }))
            .count(),
        1,
        "{:?}",
        bottom.circles
    );

    // The cylinder's own top cap at z=1.5 (inside the box): the rim-coincident plane carries
    // **two** contributions — the cap's `Seated` circle and the lateral's `Graze`, because a
    // rim is a touch, not a miss. It adds no *transversal* twin (t = span end, not strictly
    // inside). ★ On a **convex** cap the two agree about which side the body is on; the whole
    // point of carrying both is the reflex corner (a blind bore's ceiling) where they do not,
    // and `edge_mask`'s `Graze > Seated` then picks the right one.
    let cap = trace_on_class_of(
        &m,
        a,
        b,
        class_at_z(&planes, 1.5),
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    assert!(cap.declined.is_empty(), "{:?}", cap.declined);
    let sides: Vec<(bool, bool)> = cap
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some((false, body_above)),
            SegKind::Graze { body_above } => Some((true, body_above)),
            SegKind::Transversal { .. } | SegKind::Tangent { .. } => None,
        })
        .collect();
    assert_eq!(sides.len(), 2, "seated and graze: {:?}", cap.circles);
    assert!(
        sides.iter().any(|(g, _)| *g) && sides.iter().any(|(g, _)| !*g),
        "one of each kind: {:?}",
        cap.circles
    );
    assert_eq!(
        sides[0].1, sides[1].1,
        "a convex cap's two contributions agree on the body's side: {:?}",
        cap.circles
    );
}

/// ★★ **The corner the graze exists for.** On a convex cap the lateral's rim-graze and the
/// cap's seated circle say the same thing, so carrying both changes nothing — the test above
/// locks that. Here they **disagree**: at a blind bore's ceiling the plate's material is above
/// the cap while the bore's wall hangs below it, a reflex dihedral in the plane. `edge_mask`'s
/// `Graze > Seated` precedence then picks the wall's side, which is the whole mechanism of the
/// repair; before the graze existed, the seated rule flipped the wrong label bit and the next
/// boolean on that solid came back `CylinderGateUndecided`.
#[test]
fn at_a_blind_bores_ceiling_the_graze_and_the_seated_circle_disagree() {
    let mut m = Model::new();
    let plate = m.add_cuboid(Point3::from_array([0.0; 3]), Point3::from_array([10.0; 3]));
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([5.0, 5.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.0,
        3.0,
    )
    .solid;
    m.rebuild_adjacency();
    let bored = crate::boolean(&mut m, BoolKind::Cut, plate, hole).expect("bore")[0];
    // A second operand so the arrangement runs on the bored solid as an operand.
    let boss = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 10.0]),
        Point3::from_array([2.0, 2.0, 12.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, bored, boss).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let tr = trace_on_class_of(
        &m,
        bored,
        boss,
        class_at_z(&planes, 3.0), // the bore's ceiling
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    let seated: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Seated { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    let grazes: Vec<bool> = tr
        .circles
        .iter()
        .filter_map(|c| match c.kind {
            SegKind::Graze { body_above } => Some(body_above),
            _ => None,
        })
        .collect();
    assert_eq!(
        seated.len(),
        1,
        "the ceiling is seated here: {:?}",
        tr.circles
    );
    assert_eq!(grazes.len(), 1, "the wall grazes here: {:?}", tr.circles);
    assert_ne!(
        seated[0], grazes[0],
        "a reflex dihedral: the cap and the wall put the body on opposite sides — this is the \
             only place the precedence matters, and the only reason the graze is emitted at all"
    );
    // The wall hangs below the ceiling, and that is the side `edge_mask` believes.
    assert!(!grazes[0], "the bore's wall is below its ceiling");
}

/// Partial overlap is NOT merged — different endpoints mean different edges. a and b share
/// the y=1 plane; a's y=1 chord is x∈[0,2], b's is x∈[1,3] — overlapping on x∈[1,2] but not
/// coincident. They must stay separate.
#[test]
fn partial_overlap_is_not_merged() {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0, 0.0, 0.0]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([1.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let PlaneSetup {
        planes: faces_tab,
        geom: planes,
        surf_ix,
        inc_a,
        inc_b,
        plane_ix,
        standard,
        notes,
        ..
    } = plane_index_setup(&m, a, b).unwrap();
    let jd = Judge::new(&planes, standard, &notes);
    let wc = shared_cap_class(&m, a, b, &surf_ix, &faces_tab, &plane_ix);
    let tr = trace_on_class_of(
        &m,
        a,
        b,
        wc,
        &jd,
        &[],
        &faces_tab,
        &surf_ix,
        &inc_a,
        &inc_b,
        &plane_ix,
        Default::default(),
    );
    // The y=1 wall class hosts a's chord x∈[0,2] and b's chord x∈[1,3]: same wall, different
    // endpoints. After merge they remain two distinct MergedSegs (each still merging its own
    // seated≡transversal coincidence).
    let merged = merge_coincident(&jd, &tr.segs, wc, &Aliases::default());
    // The shared y=1 wall class (a face at y=1).
    let y1 = planes
        .iter()
        .position(|p| {
            p.witness_coords()
                .iter()
                .all(|q| (q.as_array()[1] - 1.0).abs() < 1e-12)
        })
        .expect("a y=1 face");
    // a's chord x∈[0,2] and b's chord x∈[1,3] ride y=1 but have different endpoints, so they
    // stay as two distinct MergedSegs. A merge that ignored extent would collapse them to one.
    let on_y1 = merged.iter().filter(|e| e.wall == y1).count();
    assert!(
        on_y1 >= 2,
        "a's and b's y=1 chords stay distinct (partial overlap not merged): {on_y1}"
    );
}
