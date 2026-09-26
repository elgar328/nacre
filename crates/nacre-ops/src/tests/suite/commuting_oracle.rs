//! The commuting oracle: a boolean commutes with every rigid motion of the group.

use super::*;

/// One motion per sign class, for the always-on run; the ignored run takes the whole group.
const MOTION_SUBSET: [&str; 6] = [
    "t(-4,-4,-2)",
    "rz90",
    "rx90",
    "ry90",
    "rz90+t",
    "t(7/11,3/10,1/4)",
];

/// The families the oracle runs: the boss corpus, an enclosed boss, a planar pair and the
/// second-operation families (the first boolean's result against a tool).
fn oracle_families() -> Vec<(&'static str, Build)> {
    let mut v: Vec<(&'static str, Build)> = Vec::new();
    for &(fam, base, h) in BOSS_FAMILIES.iter() {
        v.push((fam, Box::new(move || boss_family(base, h))));
    }
    v.push(("enclosed", Box::new(|| boss_family([2.0, 2.0, 0.5], 1.0))));
    v.push(("planar", Box::new(planar_pair)));
    v.push(("bore offset", Box::new(|| bored_plate_and_slab(12.0))));
    v.push(("bore axis", Box::new(|| bored_plate_and_slab(10.0))));
    v.push(("through mid", Box::new(fused_through_and_mid_slab)));
    v.push(("recorded block", Box::new(recorded_block_and_box)));
    v
}

/// A block whose every motion is recorded — its `Through` face, `fixtures::through_block` — and a
/// box across one of its edges. The family where the group's quarter turns, translations and
/// rigid motions all become chain nodes on one operand while the other carries them into its
/// points.
fn recorded_block_and_box() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = super::fixtures::through_block(&mut m, 1.0);
    let b = m.add_cuboid(
        Point3::from_array([1.0, 2.5, 0.5]),
        Point3::from_array([2.0, 3.5, 1.5]),
    );
    m.rebuild_adjacency();
    (m, a, b)
}

/// The census's `rot` pair: a unit-ish box and a post through it.
fn planar_pair() -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let a = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([2.0, 2.0, 2.0]),
    );
    let b = m.add_cuboid(
        Point3::from_array([0.4, 0.4, -1.0]),
        Point3::from_array([1.6, 1.6, 3.0]),
    );
    m.rebuild_adjacency();
    (m, a, b)
}

/// `bands`' bored plate (40×20×5, bore r 3 at (8, 10)) and a slab `y ≥ y0` whose wall face
/// crosses the bore — offset from the axis at `y0 = 12`, through it at `y0 = 10`.
fn bored_plate_and_slab(y0: f64) -> (Model, Handle<Solid>, Handle<Solid>) {
    let mut m = Model::new();
    let plate = m.add_cuboid(
        Point3::from_array([0.0; 3]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    let hole = crate::fixtures::cylinder_with_seam(
        &mut m,
        Point3::from_array([8.0, 10.0, -1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, -1.0, 0.0]),
        3.0,
        7.0,
    )
    .solid;
    m.rebuild_adjacency();
    let holed = boolean(&mut m, BoolKind::Cut, plate, hole).expect("the bore cuts")[0];
    m.rebuild_adjacency();
    let slab = m.add_cuboid(
        Point3::from_array([0.0, y0, 0.0]),
        Point3::from_array([40.0, 20.0, 5.0]),
    );
    m.rebuild_adjacency();
    (m, holed, slab)
}

/// The through boss fused onto its plate, and the crossing census's `mid` slab — the tool that
/// parts the result in two.
fn fused_through_and_mid_slab() -> (Model, Handle<Solid>, Handle<Solid>) {
    let (mut m, plate, boss) = boss_family([2.0, 2.0, -1.0], 4.0);
    let fused = boolean(&mut m, BoolKind::Fuse, plate, boss).expect("the through boss fuses")[0];
    m.rebuild_adjacency();
    let slab = m.add_cuboid(
        Point3::from_array([-1.0, -1.0, 0.5]),
        Point3::from_array([5.0, 5.0, 1.5]),
    );
    m.rebuild_adjacency();
    (m, fused, slab)
}

/// The unmoved answer's volume for the families the crossing census does not already lock —
/// the closed forms, so an error that is *equivariant* under motion cannot pass as commuting.
fn expected_unmoved_volume(fam: &str, kind: BoolKind) -> Option<f64> {
    use std::f64::consts::PI;
    let seg = |d: f64, r: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let bored = 4000.0 - 45.0 * PI;
    Some(match (fam, kind) {
        ("enclosed", BoolKind::Fuse) => 32.0,
        ("enclosed", BoolKind::Cut) => 32.0 - PI * 0.25,
        ("enclosed", BoolKind::Common) => PI * 0.25,
        ("planar", BoolKind::Fuse) => 8.0 + 5.76 - 2.88,
        ("planar", BoolKind::Cut) => 8.0 - 2.88,
        ("planar", BoolKind::Common) => 2.88,
        // The slab's box less the bore's `y ≥ 12` segment (d = 2, r = 3), taken from the bored plate.
        // ★ The two fuses build (the slab's wall halves the bore's rim; the pieces
        // the split keeps and the witnesses the nesting has now read it): the bored plate plus
        // the slab's share of the bore.
        ("bore offset", BoolKind::Fuse) => bored + 5.0 * seg(2.0, 3.0),
        ("bore offset", BoolKind::Cut) => bored - (1600.0 - 5.0 * seg(2.0, 3.0)),
        ("bore offset", BoolKind::Common) => 1600.0 - 5.0 * seg(2.0, 3.0),
        ("bore axis", BoolKind::Fuse) => bored + 22.5 * PI,
        ("bore axis", BoolKind::Cut) => bored - (2000.0 - 22.5 * PI),
        ("bore axis", BoolKind::Common) => 2000.0 - 22.5 * PI,
        // The fused body is `32 + π/2`; the slab holds `16` of it and `20` beside it.
        ("through mid", BoolKind::Fuse) => 32.0 + PI / 2.0 + 20.0,
        ("through mid", BoolKind::Cut) => 16.0 + PI / 2.0,
        ("through mid", BoolKind::Common) => 16.0,
        _ => return None,
    })
}

/// The unmoved answer's **name** where it is not a build — a reject the oracle keeps stable by
/// name, so the day it builds is noticed and the row gets its closed form.
fn expected_unmoved_name(fam: &str, kind: BoolKind) -> Option<&'static str> {
    // ★ Empty: the bore fuses build and are locked by volume above. The door stays for the next
    // family a reject is the honest answer for.
    let _ = (fam, kind);
    None
}

fn answers_agree(a: &Answer, b: &Answer) -> Result<(), String> {
    if a.name != b.name {
        return Err(format!("{} -> {}", a.name, b.name));
    }
    if a.volumes.len() != b.volumes.len()
        || a.volumes
            .iter()
            .zip(&b.volumes)
            .any(|(x, y)| (x - y).abs() > 1e-9 * (1.0 + x.abs()))
    {
        return Err(format!("volume {:?} -> {:?}", a.volumes, b.volumes));
    }
    if a.counts != b.counts {
        return Err(format!("counts {:?} -> {:?}", a.counts, b.counts));
    }
    if !b.valid {
        return Err("moved result does not validate".into());
    }
    Ok(())
}

/// The result's vertex coordinates as bits, sorted, each tagged by its definition's kind
/// (0 three-plane, 1 seam, 2 pierce) — a `Vec` rather than a hash so the first difference can be
/// named. Traversal order does not matter; the sort removes it.
fn sorted_vertex_bits(m: &Model, solids: &[Handle<Solid>]) -> Vec<(u8, [u64; 3])> {
    let mut out = Vec::new();
    for &s in solids {
        let src = m.solid(s).clone();
        for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
            for &fh in &m.shell(sh).faces {
                let face = m.face(fh);
                for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                    for he in &lp.half_edges {
                        for &vh in m.edge(he.edge).vertices.iter() {
                            let kind = match m.vertex(vh) {
                                Vertex::ThreePlane(_) => 0u8,
                                Vertex::OnSeam(_) => 1,
                                Vertex::Pierce { .. } => 2,
                            };
                            let p = m.vertex_point(vh).as_array();
                            out.push((kind, [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]));
                        }
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// What the commuting diagram's digest said for one cell — a statistic the ledger prints; which
/// of its arms is a lock is decided by measurement.
#[derive(Debug, Clone, PartialEq)]
enum Digest {
    NotAsked,
    /// Moving the unmoved result was refused (`OriginNotOnSolid`): nothing to compare bits with.
    Refused,
    Identical,
    /// Equal once pierce vertices are left out.
    PiercelessIdentical,
    Differs(String),
}

fn panic_text(p: Box<dyn std::any::Any + Send>) -> String {
    p.downcast_ref::<String>()
        .cloned()
        .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_else(|| "<non-string panic>".into())
}

/// **The oracle.** `Ok((outcome, digest))`; the unmoved answer is handed back too so the caller
/// can lock it against a closed form.
fn boolean_commutes(
    build: &Build,
    kind: BoolKind,
    iso: &nacre_exact::Isometry,
    class: MotionClass,
) -> (Answer, Outcome, Digest) {
    let (mut m0, a0, b0) = build();
    let r0 = boolean(&mut m0, kind, a0, b0);
    m0.rebuild_adjacency();
    let ans0 = answer(&m0, &r0);
    let (mut m, a, b) = build();
    let a = match transform(&mut m, a, iso) {
        Ok(x) => x,
        Err(e) => return (ans0, Outcome::InputUntransportable(e), Digest::NotAsked),
    };
    let b = match transform(&mut m, b, iso) {
        Ok(x) => x,
        Err(e) => return (ans0, Outcome::InputUntransportable(e), Digest::NotAsked),
    };
    m.rebuild_adjacency();
    let r = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        boolean(&mut m, kind, a, b)
    })) {
        Ok(r) => r,
        Err(p) => return (ans0, Outcome::Panicked(panic_text(p)), Digest::NotAsked),
    };
    m.rebuild_adjacency();
    let ans = answer(&m, &r);
    let outcome = match answers_agree(&ans0, &ans) {
        Ok(()) => Outcome::Commutes,
        Err(why) => Outcome::Diverged(why),
    };
    // The commuting diagram's other road: the unmoved result, moved.
    let digest = match (&r0, &r, class) {
        (Ok(v0), Ok(v), MotionClass::Quadrantal | MotionClass::ExactRigid) if !v0.is_empty() => {
            let mut moved0 = Vec::new();
            let mut refused = false;
            for &s in v0 {
                let twin = crate::transform::copy(&mut m0, s).expect("a copy of a live solid");
                match transform(&mut m0, twin, iso) {
                    Ok(t) => moved0.push(t),
                    Err(_) => refused = true,
                }
            }
            if refused {
                Digest::Refused
            } else {
                m0.rebuild_adjacency();
                let x = sorted_vertex_bits(&m0, &moved0);
                let y = sorted_vertex_bits(&m, v);
                if x == y {
                    Digest::Identical
                } else {
                    let pierceless = |v: &[(u8, [u64; 3])]| -> Vec<(u8, [u64; 3])> {
                        v.iter().filter(|e| e.0 != 2).copied().collect()
                    };
                    if pierceless(&x) == pierceless(&y) {
                        Digest::PiercelessIdentical
                    } else {
                        let first = x
                            .iter()
                            .zip(&y)
                            .find(|(p, q)| p != q)
                            .map(|(p, q)| format!("{p:?} vs {q:?}"))
                            .unwrap_or_else(|| format!("lengths {} vs {}", x.len(), y.len()));
                        Digest::Differs(first)
                    }
                }
            }
        }
        _ => Digest::NotAsked,
    };
    (ans0, outcome, digest)
}

/// A «site» that is not an assertion: the cell diverges **silently** — no check catches it and
/// the moved answer is simply a different answer. The most important rows of the ledger.
const DIVERGES: &str = "<diverges silently>";

/// **The ledger of cells the kernel does not commute on** — `(family, motions, sites)`,
/// for all three kinds (measured: every divergence here is the same for Fuse, Cut and Common).
/// A known cell must fail at one of `sites` (a set, because a boolean's classes are built in
/// parallel and two of them may die at different sites; which payload surfaces is not
/// determined). Commuting instead is red — the row must then be removed — and so is failing
/// somewhere else.
///
/// ★ The ledger is empty: every cell of the whole group commutes, the recorded quarter turns
/// included (the cylinder families refused them by name until a cylinder's chain folded to a
/// world statement). It stays here as the shape the next divergence is written in — a reject
/// under the reason's name, a silent divergence under [`DIVERGES`].
const KNOWN: &[(&str, &[&str], &[&str])] = &[];

/// The count lock: how many cells `KNOWN` names (three kinds per motion) — zero, since every
/// cell commutes.
const KNOWN_CELLS: usize = 0;

/// **Cells whose digest differs today, by name** — `(family, motions)`, for Fuse and Cut (the
/// Common of that family keeps no such vertex). They commute: volumes, counts and validity agree;
/// only the bits of one vertex do not. Under a recorded chain that folds, a vertex whose true
/// coordinate is exactly `0` is realized by replaying the chain, and the replay cannot decide the
/// zero (`RealizeError::Undecided` — the paid door as well), so the cache keeps the construction's
/// figure and the vertex's other coordinates sit an ulp off the nearest `f64` the other road
/// realizes (measured: `1.4399999999999997` against `1.44`). The chain folds to exact rationals,
/// and a realization read from the fold would decide the zero outright — the day it does, these
/// cells come out identical, go red here, and leave.
const DIGEST_DIFFERS: &[(&str, &[&str])] = &[(
    "recorded block",
    &["rx90", "rx180", "rx270", "ry90", "ry180", "ry270", "ry90+t"],
)];

/// The count lock on [`DIGEST_DIFFERS`]: two kinds per motion.
const DIGEST_DIFFERS_CELLS: usize = 14;

fn digest_differs(fam: &str, motion: &str, kind: BoolKind) -> bool {
    kind != BoolKind::Common
        && DIGEST_DIFFERS
            .iter()
            .any(|d| d.0 == fam && d.1.contains(&motion))
}

fn known_sites(fam: &str, motion: &str) -> Option<&'static [&'static str]> {
    KNOWN
        .iter()
        .find(|k| k.0 == fam && k.1.contains(&motion))
        .map(|k| k.2)
}

fn run_commuting_oracle(labels: &[&str]) {
    assert_eq!(
        KNOWN.iter().map(|k| 3 * k.1.len()).sum::<usize>(),
        KNOWN_CELLS,
        "the ledger of known cells changed size"
    );
    assert_eq!(
        DIGEST_DIFFERS.iter().map(|d| 2 * d.1.len()).sum::<usize>(),
        DIGEST_DIFFERS_CELLS,
        "the digest scoreboard changed size"
    );
    let group: Vec<_> = motion_group()
        .into_iter()
        .filter(|g| labels.contains(&g.0.as_str()))
        .collect();
    assert_eq!(
        group.len(),
        labels.len(),
        "every requested motion label exists"
    );
    let fams = oracle_families();
    for k in KNOWN {
        assert!(
            fams.iter().any(|f| f.0 == k.0),
            "KNOWN names a family: {}",
            k.0
        );
        for mn in k.1 {
            assert!(
                motion_group().iter().any(|g| g.0 == *mn),
                "KNOWN names a motion: {mn}"
            );
        }
    }
    for d in DIGEST_DIFFERS {
        assert!(
            fams.iter().any(|f| f.0 == d.0),
            "the scoreboard names a family: {}",
            d.0
        );
        for mn in d.1 {
            assert!(
                motion_group().iter().any(|g| g.0 == *mn),
                "the scoreboard names a motion: {mn}"
            );
        }
    }
    let kinds = [
        ("Fuse", BoolKind::Fuse),
        ("Cut", BoolKind::Cut),
        ("Common", BoolKind::Common),
    ];
    let mut failures: Vec<String> = Vec::new();
    // Per class: identical, pierceless-identical, differs (first), refused, not asked.
    let mut stats: HashMap<MotionClass, (usize, usize, usize, usize, Option<String>)> =
        HashMap::new();
    let mut known_hit = 0usize;
    for (fam, build) in &fams {
        for (kn, kind) in kinds {
            for (mn, iso, class) in &group {
                let (ans0, outcome, digest) = boolean_commutes(build, kind, iso, *class);
                if let Some(want) = expected_unmoved_volume(fam, kind) {
                    let got: f64 = ans0.volumes.iter().sum();
                    if (got - want).abs() > 1e-9 * (1.0 + want.abs()) {
                        failures.push(format!(
                            "{fam} {kn}: unmoved {} volume {got} vs {want}",
                            ans0.name
                        ));
                    }
                }
                if let Some(want) = expected_unmoved_name(fam, kind)
                    && ans0.name != want
                {
                    failures.push(format!("{fam} {kn}: unmoved {} vs {want}", ans0.name));
                }
                let known = known_sites(fam, mn);
                match (known, &outcome) {
                    (None, Outcome::Commutes) => {}
                    (Some(sites), Outcome::Panicked(msg))
                        if sites.iter().any(|s| msg.contains(s)) =>
                    {
                        known_hit += 1;
                    }
                    (Some(sites), Outcome::Diverged(_)) if sites.contains(&DIVERGES) => {
                        known_hit += 1;
                    }
                    // A **named** reject on the moved road: the site is the reason's name, and
                    // only that reason matches — a different reason, or a different answer, is
                    // red. `DIVERGES` would take any divergence, a wrong volume included.
                    (Some(sites), Outcome::Diverged(why))
                        if sites
                            .iter()
                            .any(|s| why.ends_with(&format!("-> Rejected({s})"))) =>
                    {
                        known_hit += 1;
                    }
                    _ => failures.push(format!(
                        "{fam} {kn} {mn}: {} (known sites: {known:?})",
                        outcome.text()
                    )),
                }
                // ★ The commuting diagram's digest is a lock where measurement said it holds
                // (the whole group: every quadrantal cell bit-identical, 600/612; every
                // exact-rigid cell identical or identical without its pierce vertices, 338/340 —
                // the pierce vertices of the offset bosses round once more under a translation,
                // as their `a + b√c` predicts; the rest are the `DIGEST_DIFFERS` scoreboard).
                if digest_differs(fam, mn, kind) {
                    if !matches!(digest, Digest::Differs(_)) {
                        failures.push(format!(
                            "{fam} {kn} {mn}: digest {digest:?}, which the scoreboard says \
                             differs — the row must go"
                        ));
                    }
                } else if known.is_none() {
                    let ok = match (class, &digest) {
                        (MotionClass::Quadrantal, Digest::Identical) => true,
                        (
                            MotionClass::ExactRigid,
                            Digest::Identical | Digest::PiercelessIdentical,
                        ) => true,
                        (MotionClass::Recorded, Digest::NotAsked) => true,
                        // Nothing to compare: no solid came out (a reject, or an empty Common).
                        (_, Digest::NotAsked) => ans0.volumes.is_empty(),
                        // Moving the unmoved result was refused — counted, and the volume arm
                        // above is what remains.
                        (_, Digest::Refused) => true,
                        _ => false,
                    };
                    if !ok {
                        failures.push(format!("{fam} {kn} {mn}: digest {digest:?}"));
                    }
                }
                if !matches!(digest, Digest::Identical | Digest::NotAsked) {
                    eprintln!("EXPD {fam} {kn} {mn}: {digest:?}");
                }
                let e = stats.entry(*class).or_insert((0, 0, 0, 0, None));
                match digest {
                    Digest::Identical => e.0 += 1,
                    Digest::PiercelessIdentical => e.1 += 1,
                    Digest::Differs(first) => {
                        e.2 += 1;
                        e.4.get_or_insert(format!("{fam} {kn} {mn}: {first}"));
                    }
                    Digest::Refused => e.3 += 1,
                    Digest::NotAsked => {}
                }
            }
        }
    }
    for (class, (i, b, d, r, first)) in &stats {
        eprintln!(
            "EXPM class {class:?} identical {i} / pierceless-identical {b} / differs {d} (first: {}) / refused {r}",
            first.as_deref().unwrap_or("-")
        );
    }
    eprintln!(
        "EXPM cells {} known-hit {known_hit} failures {}",
        fams.len() * kinds.len() * group.len(),
        failures.len()
    );
    assert!(
        failures.is_empty(),
        "the boolean does not commute with rigid motion on {} cells:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// One motion per sign class — the always-on lock.
#[test]
fn the_boolean_commutes_with_rigid_motion() {
    run_commuting_oracle(&MOTION_SUBSET);
}

/// The whole group — the ignored sweep's lock.
#[test]
#[ignore = "the whole motion group (~15 × 66 booleans); the subset runs always"]
fn the_boolean_commutes_with_every_motion_of_the_group() {
    let group = motion_group();
    let labels: Vec<&str> = group.iter().map(|g| g.0.as_str()).collect();
    run_commuting_oracle(&labels);
}

/// **The transport law, as a lock**: `transform(rigid(R, t)) ≡ transform(T) ∘ transform(R)` —
/// on the offset boss below and on a block whose every motion is recorded.
///
/// The fixture is the counterexample to an exactness probe that tests the translation on the
/// pre-turn coordinates — a boss whose seam vertex `(4.8, 0.5)`
/// is exact before the turn (`0.8`, `8.5` under `t = (−4, 8, 0)` alone) and rounds after it
/// (`4.8 + 8 = 12.8` is not an f64): one rigid operation must leave the model the two operations
/// leave — the turn absorbed into the statements, the translation recorded, the chain folding to
/// the same world cylinder — and the boolean with a plate moved the same way must not be able to
/// tell the roads apart.
#[test]
fn a_rigid_motion_behaves_as_its_two_operations() {
    use nacre_exact::{Angle, Isometry, Rat, Rotation};
    let build = |m: &mut Model| -> (Handle<Solid>, Handle<Solid>) {
        let plate = m.add_cuboid(
            Point3::from_array([0.0; 3]),
            Point3::from_array([4.0, 4.0, 2.0]),
        );
        let boss = crate::fixtures::cylinder_with_seam(
            m,
            Point3::from_array([4.3, 0.5, -1.0]),
            Vector3::from_array([0.0, 0.0, 1.0]),
            Vector3::from_array([0.0, -1.0, 0.0]),
            0.5,
            4.0,
        )
        .solid;
        m.rebuild_adjacency();
        (plate, boss)
    };
    let turn = Rotation {
        axis: Axis::Z,
        pivot: [Rat::from_int(0); 3],
        angle: Angle::from_deg(Rat::from_int(90)).unwrap(),
    };
    let shift = [Rat::from_int(-4), Rat::from_int(8), Rat::from_int(0)];
    // Road A: one rigid operation.
    let mut a = Model::new();
    let (pa, ba) = build(&mut a);
    let rigid = Isometry::rigid(turn, shift);
    let pa = transform(&mut a, pa, &rigid).unwrap();
    let ba = transform(&mut a, ba, &rigid).unwrap();
    a.rebuild_adjacency();
    // Road B: the turn, then the translation.
    let mut b = Model::new();
    let (pb, bb) = build(&mut b);
    let (rot, tr) = (Isometry::rotation(turn), Isometry::translation(shift));
    let pb = transform(&mut b, pb, &rot).unwrap();
    let pb = transform(&mut b, pb, &tr).unwrap();
    let bb = transform(&mut b, bb, &rot).unwrap();
    let bb = transform(&mut b, bb, &tr).unwrap();
    b.rebuild_adjacency();
    // Both record the translation (the datum rounds after the turn), and the chain folds to
    // one world cylinder on both roads.
    assert!(carries_motion(&a, ba) && carries_motion(&b, bb));
    let lateral = |m: &Model, s: Handle<Solid>| -> nacre_topo::CylinderDef {
        let sol = m.solid(s);
        let surf = m
            .shell(sol.outer)
            .faces
            .iter()
            .map(|&fh| m.face(fh).surface)
            .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
            .expect("the boss has a lateral");
        m.world_cylinder_def(surf)
            .expect("the chain folds to a world cylinder")
    };
    let (da, db) = (lateral(&a, ba), lateral(&b, bb));
    assert_eq!(da.origin(), db.origin());
    assert_eq!(da.dir(), db.dir());
    assert_eq!(
        da.origin().map(|x| x.to_f64()),
        [-4.5, 12.3, -1.0],
        "the folded cylinder stands where the motion put it"
    );
    assert_eq!(sorted_vertex_bits(&a, &[ba]), sorted_vertex_bits(&b, &[bb]));
    // The boolean cannot tell the roads apart.
    let ra = boolean(&mut a, BoolKind::Cut, pa, ba).expect("road A cuts");
    a.rebuild_adjacency();
    let rb = boolean(&mut b, BoolKind::Cut, pb, bb).expect("road B cuts");
    b.rebuild_adjacency();
    assert_eq!(sorted_vertex_bits(&a, &ra), sorted_vertex_bits(&b, &rb));
    assert!(nacre_validate::validate(&a).is_empty() && nacre_validate::validate(&b).is_empty());

    // The recorded road: a block whose every motion is a node (its `Through` face) takes the
    // same law — one rigid operation and its two halves leave every face the same chain, the
    // same world plane and the same vertex bits.
    let recorded = |split: bool| -> (Model, Handle<Solid>) {
        let mut m = Model::new();
        let s = super::fixtures::through_block(&mut m, 1.0);
        let s = if split {
            let s = transform(&mut m, s, &rot).unwrap();
            transform(&mut m, s, &tr).unwrap()
        } else {
            transform(&mut m, s, &rigid).unwrap()
        };
        m.rebuild_adjacency();
        (m, s)
    };
    let faces_said = |m: &Model, s: Handle<Solid>| -> Vec<String> {
        let mut out: Vec<String> = m
            .shell(m.solid(s).outer)
            .faces
            .iter()
            .map(|&f| {
                let h = m.face(f).surface;
                let mut chain = Vec::new();
                let mut at = m.plane_motion(h);
                while let Some(node) = at {
                    chain.push(match m.motion(node).motion {
                        nacre_topo::Motion::Rotate { .. } => "R",
                        nacre_topo::Motion::Translate { .. } => "T",
                        nacre_topo::Motion::Mirror { .. } => "M",
                        nacre_topo::Motion::Frame { .. } => "F",
                    });
                    at = m.motion(node).parent;
                }
                format!("{:?} {chain:?}", m.world_plane_name(h))
            })
            .collect();
        out.sort();
        out
    };
    let ((ma, sa), (mb, sb)) = (recorded(false), recorded(true));
    let said = faces_said(&ma, sa);
    assert!(
        said.iter().any(|f| f.ends_with(r#"["T", "R"]"#)),
        "the recorded road records the turn and the translation: {said:?}"
    );
    assert_eq!(said, faces_said(&mb, sb), "one chain shape, one world");
    assert_eq!(
        sorted_vertex_bits(&ma, &[sa]),
        sorted_vertex_bits(&mb, &[sb])
    );
}

/// **The half-recorded chain's own fixture**: the offset boss under `rz90 + t(5, −3, 2)`. Its
/// seam vertex `4.8 + 5` rounds *before* the turn and is exact *after* it, so the whole motion
/// carries and nothing is recorded — and `world_cylinder_def`'s postcondition (truth == cache)
/// holds.
#[test]
fn the_offset_boss_under_a_rigid_motion_keeps_one_cylinder() {
    let (mut m, plate, boss) = boss_family([4.3, 2.0, -1.0], 4.0);
    let iso = rigid_iso(Axis::Z, 90, [5, -3, 2]);
    let plate = transform(&mut m, plate, &iso).unwrap();
    let boss = transform(&mut m, boss, &iso).unwrap();
    m.rebuild_adjacency();
    assert!(
        !carries_motion(&m, boss),
        "the whole motion is exact on this data"
    );
    let sol = m.solid(boss);
    let surf = m
        .shell(sol.outer)
        .faces
        .iter()
        .map(|&fh| m.face(fh).surface)
        .find(|&h| matches!(m.surface_cache(h), nacre_geom::Surface::Cylinder(_)))
        .expect("the boss has a lateral");
    let def = m.world_cylinder_def(surf).expect("a world cylinder");
    assert_eq!(def.origin().map(|x| x.to_f64()), [3.0, 1.3, 1.0]);
    let out = boolean(&mut m, BoolKind::Cut, plate, boss).expect("the moved boss cuts");
    m.rebuild_adjacency();
    assert_eq!(out.len(), 1);
    let v = nacre_props::mass_props(&m, out[0]).unwrap().volume;
    let want = 32.0 - 0.223648;
    assert!((v - want).abs() < 1e-5, "{v} vs {want}");
}
