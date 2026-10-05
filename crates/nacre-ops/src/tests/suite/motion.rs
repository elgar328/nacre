//! Transform, rotate, re-rotate, copy and mirror: the motion forest and what a moved solid carries.

use super::*;

/// Lower corner of a solid's outer-shell vertex bounding box (for translation
/// tests: a rigid move shifts it by exactly the offset).
fn bbox_lo(m: &Model, s: Handle<Solid>) -> [f64; 3] {
    let mut lo = [f64::INFINITY; 3];
    let sh = m.solid(s).outer;
    for &fh in &m.shell(sh).faces {
        for he in &m.face(fh).outer.half_edges {
            {
                for vh in m.edge(he.edge).vertices {
                    let p = m.vertex_point(vh).as_array();
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                    }
                }
            }
        }
    }
    lo
}

fn test_iso() -> (nacre_exact::Isometry, [f64; 3]) {
    use nacre_exact::Rat;
    (
        nacre_exact::Isometry::translation([
            Rat::new(7, 2).unwrap(),
            Rat::from_int(-4),
            Rat::from_int(11),
        ]),
        [3.5, -4.0, 11.0],
    )
}

/// A cylinder tool translated onto a plate, then cutting it — **the hole-pattern idiom**. The
/// offset is non-dyadic and still exact on the tool's statements, so the move is carried into them
/// (asserted): the tool is a world cylinder at its new place, as if it had been drawn there.
#[test]
fn a_translated_tool_cuts() {
    use nacre_exact::Rat;
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([7.3, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.1,
        30.0,
    )
    .solid;
    m.rebuild_adjacency();
    let tool = transform(
        &mut m,
        tool,
        &nacre_exact::Isometry::translation([
            Rat::try_from_f64(10.7).unwrap(),
            Rat::from_int(0),
            Rat::from_int(0),
        ]),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(!carries_motion(&m, tool), "the move is carried");
    let out = boolean(&mut m, BoolKind::Cut, plate, tool).expect("the moved tool cuts");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    let want = 40.0 * 40.0 * 10.0 - std::f64::consts::PI * 2.1 * 2.1 * 10.0;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// The same total move **split in two** — each step carried into the statements, so the second
/// moves what the first wrote. Same solid, so the same volume: the two roads must agree to the
/// tolerance the oracle is stated at.
#[test]
fn a_chained_translation_folds() {
    use nacre_exact::Rat;
    let shift = |m: &mut Model, s, x: f64| {
        let s = transform(
            m,
            s,
            &nacre_exact::Isometry::translation([
                Rat::try_from_f64(x).unwrap(),
                Rat::from_int(0),
                Rat::from_int(0),
            ]),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    };
    let mut m = Model::new();
    let plate = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 40.0, 10.0]),
    );
    let tool = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([7.3, 7.3, -5.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        2.1,
        30.0,
    )
    .solid;
    m.rebuild_adjacency();
    let tool = shift(&mut m, tool, 5.2);
    let tool = shift(&mut m, tool, 5.5);
    assert!(!carries_motion(&m, tool), "both moves are carried");
    let out = boolean(&mut m, BoolKind::Cut, plate, tool).expect("the chained tool cuts");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    let want = 40.0 * 40.0 * 10.0 - std::f64::consts::PI * 2.1 * 2.1 * 10.0;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// A **bored body** translated and fused onto its twin: the plane side of the same fact — the
/// move is carried into every statement, walls and bore alike. The offset moves all three axes,
/// so no wall of either body is shared (a shared wall is the contact family, another cell).
#[test]
fn a_translated_bored_body_fuses() {
    use nacre_exact::Rat;
    let bored = |m: &mut Model| {
        let plate = crate::fixtures::cuboid(
            m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([40.0, 40.0, 10.0]),
        );
        let bore = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([7.3, 7.3, -5.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            2.1,
            30.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    };
    let mut m = Model::new();
    let a = bored(&mut m);
    let b = bored(&mut m);
    let b = transform(
        &mut m,
        b,
        &nacre_exact::Isometry::translation(
            [25.7, 3.3, 2.0].map(|c| Rat::try_from_f64(c).unwrap()),
        ),
    )
    .unwrap();
    m.rebuild_adjacency();
    assert!(!carries_motion(&m, b), "the move is carried");
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the moved body fuses");
    assert_eq!(out.len(), 1);
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    // Two bored plates, their overlap counted once, and the part of b's bore that a fills back.
    let pi = std::f64::consts::PI;
    let bore = pi * 2.1 * 2.1;
    let want = 2.0 * (16000.0 - bore * 10.0) - (14.3 * 36.7 * 8.0) + bore * 8.0;
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
}

/// ★★ **Two bored plates meeting face to face fuse.** The shared wall is interior, so the result
/// keeps no face on it — and the corners that sat on it have to *dissolve*, which they can only do
/// once each cell's bottom (and top) become **one** face. Skipping that merge whenever a member
/// carries a circle hole leaves the corners corners, naming the dropped wall, and the assembly
/// refuses the whole boolean (`VertexNamesAbsentSurface` at the corner — measured). The merge
/// carries the circles through, and this pins the result: right volume, clean `validate`,
/// watertight mesh.
///
/// ★ **The merged face holds two circle holes** (one bore per cell), so the owner-assignment loop
/// actually runs rather than falling out on a single candidate.
#[test]
fn two_bored_plates_fuse_face_to_face() {
    let cell = |m: &mut Model, x0: f64| {
        let plate = crate::fixtures::cuboid(
            m,
            Point3::from_array([x0, 0.0, 0.0]),
            Point3::from_array([x0 + 20.0, 20.0, 10.0]),
        );
        m.rebuild_adjacency();
        let bore = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([x0 + 14.0, 14.0, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            2.0,
            12.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
        m.rebuild_adjacency();
        out
    };
    let mut m = Model::new();
    let a = cell(&mut m, 0.0);
    let b = cell(&mut m, 20.0);
    let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the plates fuse");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let want = 2.0 * (20.0 * 20.0 * 10.0 - std::f64::consts::PI * 4.0 * 10.0);
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} vs {want}");
    // The corner on the dropped wall is gone: no result vertex sits at (20, 20, 0).
    let reach = m.reachable();
    assert!(
        !reach.vertices.iter().any(|&vh| {
            let p = m.vertex_point(vh).as_array();
            (p[0] - 20.0).abs() < 1e-9 && (p[1] - 20.0).abs() < 1e-9 && p[2].abs() < 1e-9
        }),
        "the straight-angle corner on the shared wall dissolved"
    );
    // And the two bores' rims both survived as holes of the merged caps.
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    assert_eq!(
        uses.values().filter(|&&n| n != 2).count(),
        0,
        "the mesh is watertight"
    );
}

/// ★★★ **A 2×2 grid: an array fused in x, then moved in y and fused again.** Every move is
/// exact on the cells' statements, so each is carried: the second generation is a world-stated
/// row, like a row drawn in place. Where a row's surfaces come from two provenances with
/// different chains (a recorded motion — a `Through` face, an overflow or a history), moving it
/// puts carriers with different chains on one corner, and the corner road answers those in the
/// world; this grid no longer builds that population (every chain here is carried away).
///
/// The oracle is **proportionality**: n cells of one part must weigh exactly n times one cell.
/// That is what says nothing was lost or double-counted at the joins.
///
/// ★ **The cell is the user's part** (one pocket, two through-bores). Measured when these moves
/// were recorded: it was the smallest set that refused `CylinderGateUndecided` with the world
/// corner road off and built with it on — a plain plate, or one with a bore or two, never landed a
/// second-generation corner where the gate had to judge it, and the relation is not monotone
/// (the user's cell cut to two pockets and four bores builds, one of those pockets with two other
/// bores refuses). Carried, it locks the proportionality of the grid.
#[test]
fn a_two_by_two_grid_fuses() {
    use nacre_exact::Rat;
    let cell = |m: &mut Model, at: [f64; 3]| {
        let plate = crate::fixtures::cuboid(
            m,
            Point3::from_array(at),
            Point3::from_array([at[0] + 86.0, at[1] + 86.0, at[2] + 71.5]),
        );
        m.rebuild_adjacency();
        let pocket = crate::fixtures::cuboid(
            m,
            Point3::from_array([at[0] + 1.0, at[1] + 38.6, at[2]]),
            Point3::from_array([at[0] + 33.2, at[1] + 58.6, at[2] + 68.5]),
        );
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
        m.rebuild_adjacency();
        let bore = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([at[0] + 17.1, at[1] + 8.9, at[2]]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            2.22,
            143.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore).expect("bore 1")[0];
        m.rebuild_adjacency();
        let bore2 = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([at[0] + 46.75, at[1] + 37.35, at[2]]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            2.34,
            143.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(m, BoolKind::Cut, out, bore2).expect("bore 2")[0];
        m.rebuild_adjacency();
        out
    };
    let shift = |m: &mut Model, s, o: [f64; 3]| {
        let s = transform(
            m,
            s,
            &nacre_exact::Isometry::translation(o.map(|c| Rat::try_from_f64(c).unwrap())),
        )
        .unwrap();
        m.rebuild_adjacency();
        s
    };
    // One cell, for the oracle.
    let one = {
        let mut m = Model::new();
        let c = cell(&mut m, [0.0; 3]);
        nacre_props::mass_props(&m, c).unwrap().volume
    };
    // The grid: (uc ∪ uc→x) ∪ (that ∪ →y).
    let mut m = Model::new();
    let a = cell(&mut m, [0.0; 3]);
    let b = cell(&mut m, [0.0; 3]);
    let b = shift(&mut m, b, [86.0, 0.0, 0.0]);
    let row = boolean(&mut m, BoolKind::Fuse, a, b).expect("the x pair fuses")[0];
    m.rebuild_adjacency();
    let c = cell(&mut m, [0.0; 3]);
    let d = cell(&mut m, [0.0; 3]);
    let d = shift(&mut m, d, [86.0, 0.0, 0.0]);
    let row2 = boolean(&mut m, BoolKind::Fuse, c, d).expect("the second row fuses")[0];
    m.rebuild_adjacency();
    // ★ The second generation: a row that already mixes two provenances, moved again.
    let row2 = shift(&mut m, row2, [0.0, 86.0, 0.0]);
    assert!(!carries_motion(&m, row2), "the move is carried");
    let out = boolean(&mut m, BoolKind::Fuse, row, row2).expect("the grid fuses");
    assert_eq!(out.len(), 1, "one solid");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    assert!(
        (v - 4.0 * one).abs() <= 1e-9 * one,
        "{v} vs {} (4 x {one})",
        4.0 * one
    );
    let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for (_, tri) in mesh.triangles.iter() {
        for k in 0..3 {
            let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
            *uses.entry((x.min(y), x.max(y))).or_default() += 1;
        }
    }
    assert_eq!(
        uses.values().filter(|&&n| n != 2).count(),
        0,
        "the mesh is watertight"
    );

    // ★★ **The grid's corners are a datum's carriers too, so the ledger says so here.** Every
    // move was carried, so every corner answers in the world (its carriers each state it) and
    // none in a frame. A datum through three of them gets a **name**, which is the road
    // `ThroughStatement` takes when the three meets agree on a frame.
    let mut world: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
    let mut framed = 0usize;
    {
        let sol = m.solid(out[0]).clone();
        let mut seen: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
        for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                for he in &m.face(fh).outer.half_edges {
                    let vh = m.he_start(*he);
                    if seen.contains(&vh) {
                        continue;
                    }
                    seen.push(vh);
                    match m.vertex_meet(vh) {
                        Some((_, None)) => world.push(vh),
                        Some((_, Some(_))) => framed += 1,
                        None => {}
                    }
                }
            }
        }
    }
    world.sort_by_key(|v| v.index());
    assert!(
        world.len() >= 3 && framed == 0,
        "a carried grid answers every corner in the world (world {}, framed {framed})",
        world.len()
    );
    let Ok(OpOutput::DatumPlane { plane, .. }) = apply(
        &mut m,
        &Operation::DatumPlane {
            def: DatumDef::ThroughVertices([world[0], world[1], world[2]]),
        },
    ) else {
        panic!("a datum through three world-answered corners")
    };
    assert!(
        m.surface_name.contains_key(&plane),
        "the datum is named, not judged"
    );
}

/// ★★★★ **A part does not care which side of itself a boss sits on.**
///
/// With the axis exactly on a wall, the `+x` and `+y` walls built and the `-x` and `-y` walls
/// answered `RingOrientation` — a **SuspectedDefect** name, on an input a script writes by
/// accident. The cause was one cross product reading two spellings of the same fact: the chord's
/// sense is made against **stored** frames, and the class's *canonical* name opposes it on half
/// the classes, so the product was right only where the two agreed. On the other half the chord's
/// sense flipped, one cell of six (the overhang below the plate) came back wound the other way,
/// and the walk found two outer contours where a one-component class has one.
///
/// ★ **The lock the code asked for.** That site's own comment said the `wc` half "is still
/// lock-invisible … so that half stands on the convention argument, recorded rather than
/// assumed". This is the fixture that sees it: the four walls have to weigh the same, and the
/// value is hand-derived rather than copied from a run.
#[test]
fn a_boss_on_any_wall_weighs_the_same() {
    // plate 4×4×2 = 32; the boss (r = 0.5, h = 4) runs clear through it, half of it inside the
    // material, so the shared part is half a cylinder of the plate's height.
    let plate = 32.0;
    let boss = std::f64::consts::PI * 0.25 * 4.0;
    let shared = 0.5 * std::f64::consts::PI * 0.25 * 2.0;
    let run = |base: [f64; 3], kind: BoolKind| -> (f64, usize) {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(&mut m, kind, a, b).expect("a boss on a wall builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty());
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight"
        );
        let sol = m.solid(out[0]).clone();
        let faces = std::iter::once(&sol.outer)
            .chain(sol.cavities.iter())
            .map(|&sh| m.shell(sh).faces.len())
            .sum();
        (nacre_props::mass_props(&m, out[0]).unwrap().volume, faces)
    };
    let mut fuse_faces: Vec<usize> = Vec::new();
    // ★★ **All four walls, all three kinds** — each a different relation to the seam's `θ = 0`,
    // and every row exact against the same closed forms.
    for base in [
        [0.0, 2.0, -1.0], // -x wall
        [4.0, 2.0, -1.0], // +x
        [2.0, 0.0, -1.0], // -y
        [2.0, 4.0, -1.0], // +y
    ] {
        let (v, f) = run(base, BoolKind::Fuse);
        assert!(
            (v - (plate + boss - shared)).abs() < 1e-12,
            "fuse {base:?}: {v}"
        );
        fuse_faces.push(f);
        let (v, _) = run(base, BoolKind::Cut);
        assert!((v - (plate - shared)).abs() < 1e-12, "cut {base:?}: {v}");
        let (v, _) = run(base, BoolKind::Common);
        assert!((v - shared).abs() < 1e-12, "common {base:?}: {v}");
    }
    // ★★ **The mirror invariant a volume cannot state.** Four congruent solids must also be built
    // out of the same number of faces: a volume can come out right while a face is split or merged
    // away, and a face wound the wrong way is exactly the shape of what went wrong here.
    assert!(
        fuse_faces.windows(2).all(|w| w[0] == w[1]),
        "the four walls give congruent solids, so their face counts must agree: {fuse_faces:?}"
    );
    // ★ The control that says the defect needed the cylinder: the same straddle with a box tool
    // built on every wall before this fix and must go on doing so.
    for (lo, hi) in [
        ([-0.5, 1.5, -1.0], [0.5, 2.5, 3.0]),
        ([3.5, 1.5, -1.0], [4.5, 2.5, 3.0]),
        ([1.5, -0.5, -1.0], [2.5, 0.5, 3.0]),
        ([1.5, 3.5, -1.0], [2.5, 4.5, 3.0]),
    ] {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = crate::fixtures::cuboid(&mut m, Point3::from_array(lo), Point3::from_array(hi));
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Fuse, a, b).expect("the box straddle builds");
        m.rebuild_adjacency();
        let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
        assert!((v - 35.0).abs() < 1e-12, "box {lo:?}: {v}");
    }
}

/// **The budget is derived, not chosen.** A circle sampled into `n` chords is approximated by the
/// inscribed regular `n`-gon, whose area is `sinc(2π/n) = 1 − (2π/n)²/6` of the true one; a
/// cylinder's lateral loses the arc-versus-chord ratio `sinc(π/n) = 1 − (π/n)²/6`, four times
/// less. Both have `n ≥ 360°/Δθ = 180` (`circle_segments`' angular term holds whatever the
/// radius), so the worst is the disk's **2.0e-4**, and a **relative 1e-3** leaves five times that
/// while sitting two hundred times below the defect this exists to catch. The error takes
/// **either sign**: a plate whose circular bite is inscribed comes out slightly *larger*.
#[test]
fn the_mesh_covers_the_faces_it_approximates() {
    let check = mesh_covers_faces;
    let plate = |m: &mut Model| {
        crate::fixtures::cuboid(
            m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        )
    };
    // The plate alone is the instrument's own control: nothing is curved, so both readings are
    // exact to rounding and a failure here would be the instrument's, not the mesh's.
    {
        let mut m = Model::new();
        let a = plate(&mut m);
        m.rebuild_adjacency();
        check("plate", &m, &[a]);
    }
    // ★ The pinched caps: a stud tangent to the cube's wall. Blind (base anchored,
    // one cap pinched) and through (`center` anchored — the user's own script — both caps
    // pinched). These are the first solids whose caps are triangulated across a bridge, and
    // this oracle is the only thing that would see a bridged cap come out the wrong shape.
    for (name, base_z) in [("blind stud", 0.0), ("through stud", -1.0)] {
        let mut m = Model::new();
        let cube = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-0.5; 3]),
            Point3::from_array([0.5; 3]),
        );
        let stud = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array([0.3, 0.0, base_z]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.2,
            2.0,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(&mut m, BoolKind::Fuse, cube, stud).expect("the tangent stud fuses");
        m.rebuild_adjacency();
        check(name, &m, &out);
    }
    for (name, base, h) in [
        ("wall -x", [0.0, 2.0, -1.0], 4.0),
        ("wall +x", [4.0, 2.0, -1.0], 4.0),
        ("wall -y", [2.0, 0.0, -1.0], 4.0),
        ("wall +y", [2.0, 4.0, -1.0], 4.0),
        ("corner", [4.0, 4.0, -1.0], 4.0),
        ("through", [2.0, 2.0, -1.0], 4.0),
        ("on top", [2.0, 2.0, 2.0], 1.0),
        ("flush", [2.0, 2.0, 0.0], 3.0),
        // The half-height wall boss: its upper cap sits inside the plate, and
        // its mirror with the lower cap inside — one lateral face each, whose chain
        // rim's arc runs **across** θ = 0 here (the reference direction is −y, outside the
        // plate), where the boss on the x = 40 wall in `bands` has a ruling *on* θ = 0.
        ("half wall", [2.0, 0.0, -1.0], 2.0),
        ("half wall, cap below", [2.0, 0.0, 1.0], 2.0),
    ] {
        for kind in [BoolKind::Fuse, BoolKind::Cut, BoolKind::Common] {
            let mut m = Model::new();
            let a = plate(&mut m);
            let b = crate::fixtures::cylinder_with_seam(
                &mut m,
                Point3::from_array(base),
                Vector3::from_array([0.0, 0.0, 1.0]),
                Vector3::from_array([0.0, -1.0, 0.0]),
                0.5,
                h,
            )
            .solid;
            m.rebuild_adjacency();
            // A population this ladder does not build yet is not this test's business; what it
            // *does* build has to be drawn correctly.
            let Ok(out) = boolean(&mut m, kind, a, b) else {
                continue;
            };
            if out.is_empty() {
                continue;
            }
            m.rebuild_adjacency();
            check(&format!("{name} {kind:?}"), &m, &out);
        }
    }
}

/// ★★★ **A boss standing on a wall has ONE lateral face, not three.**
///
/// Emitted in as many pieces as the arrangement cuts it into — the band under the plate, the
/// half-band beside it, the band above — a lateral surface draws two circles between those pieces
/// that bound **nothing** (the surface runs smooth across them): a line ringing a boss that has
/// none. The emitter builds the lateral as one **region** of its chart, so what comes out is a
/// band with one notch punched out of its side.
///
/// The negative controls are the point of the test: a boss standing on the middle of the plate
/// really *is* two lateral faces (the plate interrupts it), and a bore really is one. Both must
/// come out with their face counts untouched, or the pass is erasing boundaries that exist.
#[test]
fn a_boss_on_a_wall_has_one_lateral_face() {
    let run = |base: [f64; 3], h: f64, kind: BoolKind| -> (f64, usize, usize, Vec<usize>) {
        let mut m = Model::new();
        let a = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let b = crate::fixtures::cylinder_with_seam(
            &mut m,
            Point3::from_array(base),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            h,
        )
        .solid;
        m.rebuild_adjacency();
        let out = boolean(&mut m, kind, a, b).expect("the boss builds");
        assert_eq!(out.len(), 1, "one solid");
        m.rebuild_adjacency();
        assert!(nacre_validate::validate(&m).is_empty(), "validate {base:?}");
        let mesh = nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).expect("tess");
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for (_, tri) in mesh.triangles.iter() {
            for k in 0..3 {
                let (x, y) = (tri.vertices[k].index(), tri.vertices[(k + 1) % 3].index());
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        assert_eq!(
            uses.values().filter(|&&n| n != 2).count(),
            0,
            "the mesh is watertight {base:?}"
        );
        let sol = m.solid(out[0]).clone();
        let (mut total, mut lateral, mut holes) = (0usize, 0usize, Vec::new());
        for sh in std::iter::once(sol.outer).chain(sol.cavities.iter().copied()) {
            for fh in m.shell(sh).faces.clone() {
                total += 1;
                let f = m.face(fh);
                if matches!(m.surface_cache(f.surface), nacre_geom::Surface::Cylinder(_)) {
                    lateral += 1;
                    // Holes only: an inner loop that is one closed edge is the upper rim.
                    let rim = |l: &nacre_topo::Loop| {
                        matches!(l.half_edges[..], [he] if {
                            let [a, b] = m.edge(he.edge).vertices;
                            a == b
                        })
                    };
                    holes.push(f.inner.iter().filter(|l| !rim(l)).count());
                }
            }
        }
        (
            nacre_props::mass_props(&m, out[0]).unwrap().volume,
            total,
            lateral,
            holes,
        )
    };
    // The plate is 32; the boss (r = 0.5, h = 4) runs clear through it with half its section
    // buried, so a wall boss weighs `32 + π/4·4 − ½·π/4·2` and a corner one buries a quarter.
    let quarter = std::f64::consts::PI * 0.25;
    // ★★ **The four congruent walls come out alike.** `ref_dir` is a function of the axis alone,
    // so the four walls stand in four different relations to θ = 0; with no seam edge that
    // relation reaches nothing — the notch is a hole in the one lateral face, and there are
    // exactly two faces fewer, on every one of them.
    for (name, base, faces, notch_holes) in [
        ("-x", [0.0, 2.0, -1.0], 10, 1),
        ("+x", [4.0, 2.0, -1.0], 10, 1),
        ("-y", [2.0, 0.0, -1.0], 10, 1),
        ("+y", [2.0, 4.0, -1.0], 10, 1),
    ] {
        let (v, total, lateral, holes) = run(base, 4.0, BoolKind::Fuse);
        assert!(
            (v - (32.0 + quarter * 4.0 - 0.5 * quarter * 2.0)).abs() < 1e-12,
            "{name}: {v}"
        );
        assert_eq!(lateral, 1, "{name}: one lateral face");
        assert_eq!(holes, vec![notch_holes], "{name}: the notch is a hole");
        assert_eq!(
            total, faces,
            "{name}: exactly two faces fewer than the three pieces"
        );
    }
    // The corner: two walls cut the cylinder, so the notch's four corners name **two different**
    // wall classes — and the merge is the same one, the notch again a hole.
    let (v, total, lateral, holes) = run([4.0, 4.0, -1.0], 4.0, BoolKind::Fuse);
    assert!(
        (v - (32.0 + quarter * 4.0 - 0.25 * quarter * 2.0)).abs() < 1e-12,
        "corner: {v}"
    );
    assert_eq!((lateral, holes, total), (1, vec![1], 9), "corner");
    // ★★ **The negative controls — seams that are real must survive untouched.**
    for (name, base, h, kind, lateral, total) in [
        // The plate genuinely interrupts this one: two lateral faces, and neither has a hole.
        ("through boss", [2.0, 2.0, -1.0], 4.0, BoolKind::Fuse, 2, 10),
        ("bore", [2.0, 2.0, -1.0], 4.0, BoolKind::Cut, 1, 7),
        ("boss on top", [2.0, 2.0, 2.0], 1.0, BoolKind::Fuse, 1, 8),
        ("flush boss", [2.0, 2.0, 0.0], 3.0, BoolKind::Fuse, 1, 8),
    ] {
        let (_, t, l, holes) = run(base, h, kind);
        assert_eq!((l, t), (lateral, total), "{name} must not move");
        assert!(holes.iter().all(|&n| n == 0), "{name}: no notch");
    }
}

/// A rational translation supersedes a cuboid: rigid, so volume/area are
/// invariant and the bounding box shifts by exactly the offset; validate/tess/
/// STEP all accept the moved solid, and the input drops from `live_solids`.
#[test]
fn transform_translate_cuboid() {
    let (iso, off) = test_iso();
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let lo0 = bbox_lo(&m, c);

    let c2 = transform(&mut m, c, &iso).unwrap();
    m.rebuild_adjacency();

    assert_eq!(m.live_solids().to_vec(), vec![c2], "input superseded");
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-12,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-12, "area invariant");
    let lo1 = bbox_lo(&m, c2);
    for k in 0..3 {
        assert!(
            (lo1[k] - (lo0[k] + off[k])).abs() < 1e-12,
            "bbox shifted by offset"
        );
    }
    assert!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok(),
        "moved solid tessellates"
    );
    assert!(
        nacre_step::to_step(&m, fixtures::STAMP)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP"),
        "moved solid exports to STEP"
    );
}

/// Transforming a boolean *result* keeps its seam vertices' three-plane definitions — the
/// failure this guards is a vertex losing its definition in the move.
///
/// A vertex carries no motion of its own: its definition names three planes, and the planes
/// carry the motion (a recorded node, or planes remapped in place when the motion records
/// nothing). Either way the truth survives; only its spelling depends on whether the motion was
/// worth recording.
#[test]
fn transform_translate_preserves_discovered_definition() {
    let (iso, _) = test_iso();
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    let before = nacre_props::mass_props(&m, r).unwrap().volume;
    let disc = count_discovered(&m, r);
    assert!(disc > 0, "the Cut result must have realized seam vertices");

    let r2 = transform(&mut m, r, &iso).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    // Every seam vertex still names its three-plane definition — directly, or through the
    // base of the motion that moved it.
    let named = boundary_verts(&m, r2)
        .into_iter()
        // A measured tolerance survives an exact move and is dropped by a recorded one
        // so what "still named" means here is the definition: every seam
        // vertex names three planes, whichever road it travelled.
        .filter(|&vh| matches!(*m.vertex(vh), Vertex::ThreePlane(_)))
        .count();
    assert_eq!(named, disc, "seam definitions preserved");
    let after = nacre_props::mass_props(&m, r2).unwrap().volume;
    assert!((after - before).abs() < 1e-12, "volume invariant");
}

/// Replay determinism: the same construction and motions reproduce the same geometry **and** the
/// same handle, down to the index — for a rigid motion with a translation in it, a lone tilted
/// rotation, and a rotation of a rotated solid.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is legal
/// only because `Handle`'s equality is its index — the very premise `replay` relies on.
/// `tests/invariants/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn a_chain_of_motions_is_deterministic() {
    use nacre_exact::Axis;
    let chains: [(&str, Vec<nacre_exact::Isometry>); 3] = [
        ("rigid motion", vec![test_iso().0]),
        ("30° about Z", vec![rot30()]),
        (
            "30° about Z, then 45° about X",
            vec![rot_iso(Axis::Z, 30), rot_iso(Axis::X, 45)],
        ),
    ];
    for (name, chain) in &chains {
        let build = || {
            let mut m = Model::new();
            let mut s = crate::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 3.0, 4.0]),
            );
            for iso in chain {
                s = transform(&mut m, s, iso).unwrap();
            }
            (bbox_lo(&m, s), s)
        };
        assert_eq!(
            build(),
            build(),
            "{name}: same ops → same geometry and handle"
        );
    }
}

/// The bit-identity guard, on the one producer that can break it.
///
/// A definition and its cached coordinate must agree **exactly**, and they do only because the
/// replay performs the very same float operations the producer did. A reflection is the step
/// where that is easiest to lose — `WitnessPoint::mirror` must walk the same `2c − x` that
/// `AxisMirror::point` just walked — and a chain that ends in one is what this checks.
#[test]
fn a_mirrored_rotated_vertex_reconstructs_from_its_definition() {
    use nacre_exact::{Axis, Rat};
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let r = transform(&mut m, c, &rot30()).unwrap();
    m.rebuild_adjacency();
    let mirrored = crate::transform::mirror(&mut m, r, Axis::X, Rat::from_int(0)).unwrap();
    m.rebuild_adjacency();

    let shell = m.solid(mirrored).outer;
    let mut checked = 0;
    for &fh in &m.shell(shell).faces.clone() {
        for he in &m.face(fh).outer.half_edges.clone() {
            for vh in m.edge(he.edge).vertices.iter() {
                // The image keeps its definition, and the definition reproduces the
                // coordinate: solve the corner's three planes in the frame their names are
                // stated in, replay the chain the *faces* record (the vertex has no
                // motion of its own), and the answer is the stored coordinate bit for bit.
                let Vertex::ThreePlane(tri) = *m.vertex(*vh) else {
                    unreachable!("a cuboid corner is a three-plane point")
                };
                let motion_of = |h| match m.surface(h) {
                    nacre_topo::Surface::Plane { motion, .. } => *motion,
                    nacre_topo::Surface::Cylinder { motion, .. } => *motion,
                };
                // The caps were fixed by the Z turn and by the X mirror (normal ⊥ both),
                // so they are world-stated; the moved carriers share the one leaf whose
                // chain the fixed caps provably survive — the reconstruction below runs
                // through the walls' recorded [rotate, mirror] history unchanged.
                let leaves: Vec<_> = tri.iter().filter_map(|&h| motion_of(h)).collect();
                let rotation = *leaves
                    .first()
                    .expect("a mirrored image's corner names a moved carrier");
                assert!(
                    leaves.iter().all(|&l| l == rotation),
                    "the corner's moved planes share one motion leaf"
                );
                let mut coeffs = [[Rat::from_int(0); 4]; 3];
                for (o, h) in coeffs.iter_mut().zip(tri) {
                    *o = *m.surface_name.get(&h).unwrap().narrow().unwrap();
                }
                let base = nacre_exact::three_planes_rat(coeffs)
                    .expect("three distinct planes meet in a point");
                let replayed = crate::rotated_vertex::replay_chain_coord(
                    &m,
                    [base[0].to_f64(), base[1].to_f64(), base[2].to_f64()],
                    rotation,
                )
                .expect("the root coordinate lifts to an exact rational");
                // ★ The stored coordinate is the realization of the definition; the
                // f64 replay agrees with it to its own rounding, not bit for bit.
                let stored = m.vertex_point(*vh).as_array();
                assert!(
                    matches!(m.vertex_cache(*vh), nacre_topo::PointCache::Bounded { .. }),
                    "a mirrored, rotated corner is realized from its definition"
                );
                let (realized, _) = crate::realize_vertex(&m, *vh, crate::Precision::NearestF64)
                    .expect("the def road realizes what it solves")
                    .to_f64()
                    .expect("decided");
                assert_eq!(
                    realized, stored,
                    "the cache is the realization, bit for bit"
                );
                for k in 0..3 {
                    assert!(
                        (replayed[k] - stored[k]).abs() <= 1e-9 * (1.0 + stored[k].abs()),
                        "f64 replay within the construction epsilon: {replayed:?} vs {stored:?}"
                    );
                }
                checked += 1;
            }
        }
    }
    assert!(checked > 0, "walked no vertices");
}

/// Replay determinism for the one additive operation: a copy reproduces the same
/// geometry *and* the same handle index, and leaves the same live set behind it.
///
/// The `assert_eq!` below compares handles minted by two *different* `Model`s, which is
/// legal only because `Handle`'s equality is its index — the very premise `replay` now
/// relies on. `tests/invariants/replay.rs` measures that premise directly instead of assuming it.
#[test]
fn copy_is_deterministic() {
    let build = || {
        let mut m = Model::new();
        let c = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 3.0, 4.0]),
        );
        let twin = crate::transform::copy(&mut m, c).unwrap();
        (bbox_lo(&m, twin), twin, m.live_solids().to_vec())
    };
    assert_eq!(build(), build(), "same ops → same geometry and handles");
}

/// A genuinely tilted rigid rotation: 30° about Z through the rational axis
/// point (1,1,0). Non-90° and non-axis-aligned, so it records a motion node (unlike the 90°
/// family, which stays exact).
fn rot30() -> nacre_exact::Isometry {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(1), Rat::from_int(1), Rat::from_int(0)],
        angle: Angle::from_deg(Rat::from_int(30)).unwrap(),
    })
}

/// A non-90° rotation genuinely tilts the solid: rigid (volume/area invariant),
/// validate/tess/STEP clean, a known corner lands at its exact rotated image, the
/// faces record their motion (`solid_is_rotated`), and a boolean against it
/// runs (a *mixed*-rotation cut: rotated `c2` minus an axis-aligned `d` it contains, so
/// `d` becomes a cavity).
#[test]
fn transform_rotate_cuboid_tilts_and_cuts() {
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let c2 = transform(&mut m, c, &rot30()).unwrap();
    m.rebuild_adjacency();

    assert_eq!(m.live_solids().to_vec(), vec![c2], "input superseded");
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-9,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-9, "area invariant");
    assert!(
        nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok(),
        "rotated solid tessellates"
    );
    assert!(
        nacre_step::to_step(&m, fixtures::STAMP)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP"),
        "rotated solid exports to STEP"
    );

    assert!(solid_is_rotated(&m, c2), "the faces record their rotation");

    // Corner (0,0,0) rotates about pivot (1,1) by 30°: dx=dy=-1, so
    // x' = 1 - cos30 + sin30, y' = 1 - sin30 - cos30, z' = 0.
    let (c30, s30) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let want = [1.0 - c30 + s30, 1.0 - s30 - c30, 0.0];
    let pts = outer_points(&m, c2);
    assert!(
        pts.iter()
            .any(|p| p.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9)),
        "corner (0,0,0) rotated to its exact image {want:?}; got {pts:?}"
    );

    // A mixed-rotation cut now runs: the axis-aligned `d` sits inside the rotated `c2`, so
    // `Cut(c2, d)` leaves `c2` with `d` carved out as a cavity (volume 24 − 1 = 23).
    let d = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5; 3]),
        Point3::from_array([1.5; 3]),
    );
    let r = boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "mixed cut is valid"
    );
    assert_eq!(
        m.solid(r).cavities.len(),
        1,
        "the contained box is a cavity"
    );
    let vol = nacre_props::mass_props(&m, r).unwrap().volume;
    assert!((vol - 23.0).abs() < 1e-9, "volume {vol}");
}

/// A 90° rotation about Z is axis-aligned and exact: the solid records no motion
/// (`solid_is_rotated` false), volume is exact, and a following cut succeeds and
/// validates.
#[test]
fn transform_rotate_90_is_exact_and_allows_boolean() {
    use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
    let rot90 = Isometry::rotation(Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    });
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c2 = transform(&mut m, c, &rot90).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    assert!(!solid_is_rotated(&m, c2), "90° stays exact (no motion)");
    assert_eq!(
        nacre_props::mass_props(&m, c2).unwrap().volume,
        24.0,
        "exact volume"
    );

    // c rotated 90° about origin occupies x∈[-3,0], y∈[0,2], z∈[0,4].
    // Cut with d = [-1,0.5,1]-[0.5,1.5,2]: overlap volume 1 → 24 − 1 = 23.
    let d = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, 0.5, 1.0]),
        Point3::from_array([0.5, 1.5, 2.0]),
    );
    let r = boolean(&mut m, BoolKind::Cut, c2, d).expect("exact rotation → boolean allowed");
    m.rebuild_adjacency();
    assert!(nacre_validate::validate(&m).is_empty());
    let vol = nacre_props::mass_props(&m, r[0]).unwrap().volume;
    assert!((vol - 23.0).abs() < 1e-9, "cut volume {vol}");
}

/// Rotating a boolean *result*: validate stays clean, volume is invariant, and a seam vertex
/// still names three planes — the moved ones (the seam definition is preserved through the
/// rotation).
#[test]
fn transform_rotate_boolean_result_keeps_discovered_base() {
    let (mut m, a, b) = two_boxes();
    let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
    assert!(count_discovered(&m, r) > 0, "Cut result has realized seams");
    let before = nacre_props::mass_props(&m, r).unwrap().volume;

    let r2 = transform(&mut m, r, &rot30()).unwrap();
    m.rebuild_adjacency();

    assert!(nacre_validate::validate(&m).is_empty());
    assert!(
        solid_is_rotated(&m, r2),
        "rotated result records its motion"
    );
    let after = nacre_props::mass_props(&m, r2).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");

    // A seam vertex of the *moved* result still names three planes, and those planes are
    // the moved ones (the vertex follows its faces' motion — there is no base vertex
    // left to chase).
    let sh = m.solid(r2).outer;
    let own_surfaces: std::collections::HashSet<_> = m
        .shell(sh)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .collect();
    let mut found_moved_carrier = false;
    for &fh in &m.shell(sh).faces {
        for he in &m.face(fh).outer.half_edges {
            for vh in m.edge(he.edge).vertices {
                let Vertex::ThreePlane(tri) = *m.vertex(vh) else {
                    continue;
                };
                // The carriers are the result's own planes — never a superseded
                // pre-move twin. (A restated cap's handle *equals* its pre-move
                // handle by intern-back, which is why this is a subset check on the
                // live faces' surfaces rather than "all carriers moved".)
                for h in tri {
                    assert!(
                        own_surfaces.contains(&h),
                        "a vertex carrier is not one of the solid's own planes"
                    );
                }
                if tri.iter().any(|&h| {
                    !matches!(
                        m.surface(h),
                        nacre_topo::Surface::Plane { motion: None, .. }
                    )
                }) {
                    found_moved_carrier = true;
                }
            }
        }
    }
    assert!(
        found_moved_carrier,
        "a moved result's vertices name its moved planes"
    );
}

/// A rotation `Transform` flows through `apply` and the result records its motion.
#[test]
fn transform_rotate_op_applies() {
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let out = apply(
        &mut m,
        &Operation::Transform {
            solid: c,
            isometry: rot30(),
        },
    )
    .unwrap();
    let OpOutput::Transform { solid } = out else {
        panic!("expected Transform output, got {out:?}");
    };
    assert_eq!(m.live_solids().to_vec(), vec![solid]);
    assert!(solid_is_rotated(&m, solid));
}

fn boundary_verts(m: &Model, s: Handle<Solid>) -> Vec<Handle<Vertex>> {
    let mut seen = std::collections::HashSet::new();
    let mut vs = Vec::new();
    let sh = m.solid(s).outer;
    for &fh in &m.shell(sh).faces {
        for he in &m.face(fh).outer.half_edges {
            {
                for vh in m.edge(he.edge).vertices {
                    if seen.insert(vh) {
                        vs.push(vh);
                    }
                }
            }
        }
    }
    vs
}

/// Whether any wall of `s` is a **rotated image** — asked of the surfaces that record it, not
/// of the vertices.
fn solid_is_rotated(m: &Model, s: Handle<Solid>) -> bool {
    m.shell(m.solid(s).outer).faces.iter().any(|&fh| {
        matches!(
            m.surface(m.face(fh).surface),
            nacre_topo::Surface::Plane {
                motion: Some(_),
                ..
            } | nacre_topo::Surface::Cylinder {
                motion: Some(_),
                ..
            }
        )
    })
}

/// Chain: `(node_count, axes-root-to-leaf)` of the **fullest** face history of `s` — the
/// walls', which go through every motion. Not `faces[0]`: that is a cap, and a cap a
/// rotation fixes is restated world-side (no motion, or a shorter chain
/// begun by a later non-fixing motion) — the forest's story lives on the faces that
/// genuinely moved. `None` when no face carries a rotation.
///
/// ★ There is no third field `base_is_rotated` — "does this vertex's base vertex itself
/// carry a motion?", the one-hop invariant that kept a replay from applying the same motion
/// twice. There is no base vertex any more (a moved corner is the intersection of its moved
/// planes), so double application is **unrepresentable** rather than merely untrue: the type
/// absorbed the invariant, and a probe field that could only ever read `false` would be
/// theatre.
fn forest_probe(m: &Model, s: Handle<Solid>) -> Option<(usize, Vec<nacre_exact::Axis>)> {
    let mut best: Option<Vec<nacre_exact::Axis>> = None;
    for &fh in &m.shell(m.solid(s).outer).faces {
        let &nacre_topo::Surface::Plane {
            motion: Some(rotation),
            ..
        } = m.surface(m.face(fh).surface)
        else {
            continue;
        };
        let mut axes = Vec::new();
        let mut cur = Some(rotation);
        while let Some(h) = cur {
            let n = m.motion(h);
            if let nacre_topo::Motion::Rotate { axis, .. } = n.motion {
                axes.push(axis);
            }
            cur = n.parent;
        }
        axes.reverse();
        if best.as_ref().is_none_or(|b| axes.len() > b.len()) {
            best = Some(axes);
        }
    }
    best.map(|axes| (axes.len(), axes))
}

/// **A chained boolean must not lose exactness.**
///
/// A rotated solid's face coordinates are rounded, so its planes are truthful only through an
/// exact *definition* (`FaceInfo::tri_pt3` built from a rotation history). A boolean's *result*
/// is just as rotated as its operands — but the result carries no rotation provenance, so
/// describing its faces by a tol-0 witness of the rounded triangle treats a rounded copy as the
/// truth, and one wall becomes two plane classes on the next operation.
///
/// The invariant: **every face of a boolean between rotated operands is described by a
/// rotation definition, not by its rounded coordinates.**
#[test]
fn a_chained_boolean_keeps_its_faces_exact() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([-1.0, -1.0, 0.0]),
        Point3::from_array([1.0, 1.0, 3.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5, -0.2, 1.0]),
        Point3::from_array([4.0, 0.2, 3.0]),
    );
    m.rebuild_adjacency();
    let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    // The operands hold up: a rotated solid's faces do carry definitions. The caps'
    // *truth* is world-stated since the invariant-plane restatement, but the judging
    // table re-chains a chain-fixed plane (the mirror takes the strongest description),
    // so the whole table reads rotated — exactly as it did when the twin carried the
    // motion itself.
    for (name, s) in [("operand a", a), ("operand b", b)] {
        let planes = crate::planes::collect_planes(&m, s).expect("planes");
        assert!(
            planes.iter().all(|f| f.plane().rotated),
            "{name}: a rotated operand's faces must be described by their rotation"
        );
    }
    let r = boolean(&mut m, BoolKind::Fuse, a, b).expect("fuse")[0];
    m.rebuild_adjacency();
    let planes = crate::planes::collect_planes(&m, r).expect("planes");
    let described = planes.iter().filter(|f| f.plane().rotated).count();
    assert_eq!(
        described,
        planes.len(),
        "the result of a rotated boolean is rotated too: {described}/{} faces carry a \
             definition, the rest are rounded coordinates declared exact",
        planes.len()
    );
}

/// **The same motion, applied twice, is the same node.**
///
/// The motion handle is the canonical name of "which motion", and judgments use it to decide
/// whether a whole judgement can be answered exactly in the pre-motion frame. Without
/// interning, turning two solids by the same 30° makes two nodes, their shared motion stops
/// cancelling, and rotating a model turns its exact questions into assumed ones — which is
/// what `a_shared_rotation_still_assumes_nothing` measures. (The identity it replaced was a
/// 64-bit hash of the chain's contents, where a collision would have answered a *different*
/// question with full confidence.)
#[test]
fn the_same_motion_applied_twice_is_one_node() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([2.0, 0.0, 0.0]),
        Point3::from_array([3.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let a = transform(&mut m, a, &rot_iso(Axis::X, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::X, 30)).unwrap();
    m.rebuild_adjacency();
    fn leaf(m: &Model, s: Handle<Solid>) -> Handle<nacre_topo::MotionNode> {
        let sh = m.solid(s).outer;
        let fh = m.shell(sh).faces[0];
        match m.surface(m.face(fh).surface) {
            nacre_topo::Surface::Plane {
                motion: Some(motion),
                ..
            } => *motion,
            other => panic!("a rotated solid's faces record their motion, got {other:?}"),
        }
    }
    assert_eq!(leaf(&m, a), leaf(&m, b), "one motion, one node");
    // …and a *different* motion is a different node, or the identity would be worthless.
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([5.0, 0.0, 0.0]),
        Point3::from_array([6.0, 1.0, 1.0]),
    );
    m.rebuild_adjacency();
    let c = transform(&mut m, c, &rot_iso(Axis::X, 31)).unwrap();
    m.rebuild_adjacency();
    assert_ne!(
        leaf(&m, a),
        leaf(&m, c),
        "different motions must not share a node"
    );
}

/// **A boolean result rotated again continues its history — per wall.**
///
/// One solid does not have one rotation history. A result's vertices carry no motion, so
/// asking them "what rotation is this solid at" answers `None` and the next rotation would
/// start a fresh root — replaying a pre-first-rotation witness through only the *second*
/// rotation, which is a plane that does not exist. And its walls can come from operands
/// rotated by different angles, so there is no single answer to give. Each surface therefore
/// chains from its own leaf.
#[test]
fn a_rerotated_boolean_result_continues_each_walls_history() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let a = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([1.0, 1.0, 0.0]),
        Point3::from_array([3.0, 3.0, 2.0]),
    );
    m.rebuild_adjacency();
    // Different angles, so the two operands' walls carry genuinely different histories.
    let a = transform(&mut m, a, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let b = transform(&mut m, b, &rot_iso(Axis::Z, 50)).unwrap();
    m.rebuild_adjacency();
    let r = boolean(&mut m, BoolKind::Fuse, a, b).unwrap()[0];
    m.rebuild_adjacency();
    let r = transform(&mut m, r, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    let mut leaves = std::collections::HashSet::new();
    let mut caps = 0usize;
    let sh = m.solid(r).outer;
    for &fh in &m.shell(sh).faces {
        let s = m.face(fh).surface;
        match m.surface(s) {
            &nacre_topo::Surface::Plane {
                motion: Some(rotation),
                ..
            } => {
                assert_eq!(
                    crate::rotated_vertex::motion_chain(&m, rotation)
                        .expect("an axis-aligned history holds no frame")
                        .len(),
                    2,
                    "both rotations, once each"
                );
                leaves.insert(rotation);
            }
            // Every rotation in this chain is about Z, so the z-caps are fixed by all
            // of them and stay restated — the only faces allowed to carry nothing.
            nacre_topo::Surface::Plane { motion: None, .. } => {
                let n = m.surface_name.get(&s).unwrap().narrow().unwrap();
                assert!(
                    n[0] == nacre_exact::Rat::from_int(0) && n[1] == nacre_exact::Rat::from_int(0),
                    "only a Z-fixed cap may go without a history here"
                );
                caps += 1;
            }
            other => panic!("unexpected truth {other:?}"),
        }
    }
    assert_eq!(leaves.len(), 2, "the two operands' histories stay apart");
    assert!(caps > 0, "the restated caps are present in the result");
}

/// **A copy of a rotated solid is still exactly defined.**
///
/// `copy` is `transform` under a *zero* translation, so a surface rule that reads "any
/// translation makes a rotated plane inexpressible" swallows it — and then a copy of a rotated
/// solid cannot take part in a boolean at all, though its geometry is bit-identical to the
/// original's. The identity is not a translation.
#[test]
fn a_copy_of_a_rotated_solid_answers_like_the_original() {
    use nacre_exact::Axis;
    let build = |use_copy: bool| {
        let mut m = Model::new();
        let hub = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([-1.0, -1.0, 0.0]),
            Point3::from_array([1.0, 1.0, 3.0]),
        );
        let fin = crate::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.5, -0.2, 1.0]),
            Point3::from_array([4.0, 0.2, 3.0]),
        );
        m.rebuild_adjacency();
        let mut fin = transform(&mut m, fin, &rot_iso(Axis::Z, 30)).unwrap();
        m.rebuild_adjacency();
        if use_copy {
            fin = crate::transform::copy(&mut m, fin).unwrap();
            m.rebuild_adjacency();
        }
        let r = boolean(&mut m, BoolKind::Fuse, hub, fin).expect("fuse");
        nacre_props::mass_props(&m, r[0]).expect("props").volume
    };
    // Bit-identical, not merely close: the copy walks the same definitions.
    assert_eq!(build(true), build(false));
}

/// A boolean produces seam vertices at exact axis-aligned intersections. A quadrantal
/// rotation is realized exactly (an inexact 90° leaves an ~8e-17 f64 residual, which
/// `VertexOffSurface` catches), so the residual stays 0 and validate is clean — both for a
/// pure 90° rotation and for a rigid 90°+translation (the offset cancels in
/// vertex−plane, so it does not reintroduce a residual).
#[test]
fn boolean_result_rotated_90_validates() {
    use nacre_exact::Axis;
    for iso in [rot_iso(Axis::Z, 90), rigid_iso(Axis::Z, 90, [5, -3, 2])] {
        let (mut m, a, b) = two_boxes();
        let r = boolean_one(&mut m, BoolKind::Cut, a, b).unwrap();
        let before = nacre_props::mass_props(&m, r).unwrap().volume;
        let r2 = transform(&mut m, r, &iso).unwrap();
        m.rebuild_adjacency();
        let vs = nacre_validate::validate(&m);
        assert!(
            vs.is_empty(),
            "exact 90° realization → validate clean: {vs:?}"
        );
        let after = nacre_props::mass_props(&m, r2).unwrap().volume;
        assert!((after - before).abs() < 1e-12, "volume invariant");
    }
}

/// A 90°-family rotation lands a cuboid's vertices exactly on the axis-aligned grid
/// (no ~6e-17 spurious offset): the corner (2,3,4) rotated 90° about Z maps to
/// exactly (-3,2,4).
#[test]
fn rotate_90_lands_vertices_exactly() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c2 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
    m.rebuild_adjacency();
    let pts = outer_points(&m, c2);
    // (2,3,4) about Z by 90°: (x,y)→(-y,x) → (-3,2,4). Bit-exact.
    assert!(
        pts.iter().any(|p| *p == [-3.0, 2.0, 4.0]),
        "corner lands exactly on the grid; got {pts:?}"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// Re-rotating about the same axis chains a second forest node onto the first
/// (same-axis nodes are not bundled): the leaf's parent is the earlier rotation, and the
/// solid remains a rigid (volume/area-invariant) rotated solid that validate/tess/STEP accept
/// and a boolean runs against.
#[test]
fn rerotate_same_axis_chains() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap();
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::Z])),
        "two Z nodes chained"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap();
    assert!(
        (after.volume - before.volume).abs() < 1e-9,
        "volume invariant"
    );
    assert!((after.area - before.area).abs() < 1e-9, "area invariant");
    assert!(nacre_tess::tessellate(&m, &nacre_tess::TessConfig::default()).is_ok());
    assert!(
        nacre_step::to_step(&m, fixtures::STAMP)
            .unwrap()
            .contains("MANIFOLD_SOLID_BREP")
    );
    assert!(solid_is_rotated(&m, c2), "re-rotated solid stays rotated");
    // A boolean against the chain-rotated solid runs: the axis-aligned
    // `d` inside the re-rotated `c2` is carved out, and the result is a valid solid.
    let d = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.5; 3]),
        Point3::from_array([1.5; 3]),
    );
    boolean_one(&mut m, BoolKind::Cut, c2, d).unwrap();
    m.rebuild_adjacency();
    assert!(
        nacre_validate::validate(&m).is_empty(),
        "chain-rotated cut is valid"
    );
}

/// Re-rotating about a different axis chains a node whose parent is the first
/// rotation (root → Z → X), `base` still the root.
#[test]
fn rerotate_different_axis_chains() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap().volume;
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::X])),
        "chain root→Z→X"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c2).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");
}

/// An **exact** (90°-family) rotation applied to an already-rotated solid is still
/// recorded as a chain node — the composite is inexact (an ancestor is), so the
/// forest must stay complete. The solid stays rotated.
#[test]
fn rerotate_exact_after_inexact_records_node() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 90)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c2),
        Some((2, vec![Axis::Z, Axis::X])),
        "the exact 90°X is recorded as a chain node, not dropped"
    );
    assert!(
        solid_is_rotated(&m, c2),
        "composite is inexact → still rotated"
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A translation between two same-axis rotations forces a chain (this cell never
/// bundles anyway): the forest records both rotations, `base` stays the root, and
/// the result is rigid and valid — sound with no adjacency guard (each rotation is
/// its own node).
#[test]
fn rerotate_across_translation_chains() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let before = nacre_props::mass_props(&m, c).unwrap().volume;
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &translate_iso([5, -3, 2])).unwrap();
    m.rebuild_adjacency();
    let c3 = transform(&mut m, c2, &rot_iso(Axis::Z, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c3),
        Some((2, vec![Axis::Z, Axis::Z])),
        "both rotations recorded; the intervening translation is not in the forest"
    );
    assert!(nacre_validate::validate(&m).is_empty());
    let after = nacre_props::mass_props(&m, c3).unwrap().volume;
    assert!((after - before).abs() < 1e-9, "volume invariant");
}

/// A three-axis chain records three nodes root→Z→X→Y.
#[test]
fn rerotate_deep_chain() {
    use nacre_exact::Axis;
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 3.0, 4.0]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    let c2 = transform(&mut m, c1, &rot_iso(Axis::X, 45)).unwrap();
    m.rebuild_adjacency();
    let c3 = transform(&mut m, c2, &rot_iso(Axis::Y, 20)).unwrap();
    m.rebuild_adjacency();

    assert_eq!(
        forest_probe(&m, c3),
        Some((3, vec![Axis::Z, Axis::X, Axis::Y])),
    );
    assert!(nacre_validate::validate(&m).is_empty());
}

/// A fresh rotation of an unmoved solid: an inexact angle records a single root node; a
/// 90°-family angle records none.
#[test]
fn fresh_rotation_of_constructed_unchanged() {
    use nacre_exact::Axis;
    // inexact → one root node over the cuboid.
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 30)).unwrap();
    m.rebuild_adjacency();
    assert_eq!(forest_probe(&m, c1), Some((1, vec![Axis::Z])));

    // exact 90° → no node.
    let mut m = Model::new();
    let c = crate::fixtures::cuboid(
        &mut m,
        Point3::from_array([0.0; 3]),
        Point3::from_array([1.0; 3]),
    );
    let c1 = transform(&mut m, c, &rot_iso(Axis::Z, 90)).unwrap();
    m.rebuild_adjacency();
    assert!(!solid_is_rotated(&m, c1), "fresh 90° records no motion");
    assert!(nacre_validate::validate(&m).is_empty());
    assert_eq!(forest_probe(&m, c1), None, "no rotation node");
}
