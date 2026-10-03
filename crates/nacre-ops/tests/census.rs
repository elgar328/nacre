//! **The bit census — a baseline to diff across a change that moves coordinates.**
//!
//! Prints one line per case: the exact bits of volume/area/centroid **and a hash of the result's
//! sorted vertex coordinates**. The derived quantities alone are not enough: a 1-ULP coordinate
//! change can cancel out of a volume integral, so a census of volumes can pass while the geometry
//! moved. The coordinate hash is what makes "nothing changed" checkable.
//!
//! A dump to diff, and two locks asserted as each row is recorded — plane senses and meshing (see
//! `record`). The gates run it by name, in three builds (`overview.md` 「관문」). Run it on two commits
//! and `diff`:
//!
//! ★ **World-plane seeding**: the seeds intern with every origin-touching producer, so the
//! survivor's f64 plane cache — and with it the `in:` plane digest — is the seed's on the
//! origin-touching lines. `stat seeded_hits` below is the falsifiability bridge for that
//! population.
//!
//! ```text
//! cargo test -p nacre-ops --release --test census -- --ignored --nocapture | grep '^c '
//! ```
//!
//! ★★★ **Diff it across *profiles* too, not only across commits.** Drop `--release` and the same
//! lines must come out — they do, measured (398 of them at one count; the claim is *the same
//! lines*, not a number, because the corpus grows). That is not a formality: `Angle`'s f64 route is
//! `(deg.to_f64() * PI / 180.0).cos()`, and LLVM evaluates that at compile time wherever it can see
//! the angle, one ulp away from what libm returns at run time. Two builds disagreeing here would
//! mean the coordinates a model stores depend on how it was compiled.
//!
//! ★★ **What keeps them equal is that an angle crosses the model store** before anything realizes
//! it — written into an `Operation`, pushed, read back — and no optimiser propagates a constant
//! through a heap structure. So this is a real check with a real way to fail, and no in-process
//! test can express it: a test runs in one profile. It belongs here, next to the diff it extends.
//!
//! ★★★★★ **What a census can and cannot see.** Every case is built from coordinates written here,
//! so it sees a change that moves a coordinate — and it is **blind by construction** to a change in
//! a population it does not contain. The first 130 lines are all *short* decimals, which is why a
//! regression that closed the exact path for every arbitrarily-tilted plane crossed this file
//! bit-identical, twice. The `fw` (full-width coordinates) and `tp` (tilted sketch plane) families
//! exist for that reason, and so does `mot` (the boss corpus under rigid motion):
//! without it, no row holds a **rotated** cylinder boolean at all. **Read "bit-identical"
//! as evidence only about the population present.**

use nacre_math::Point3;
use nacre_ops::{BoolKind, OpOutput, Operation, apply, boolean};
#[path = "support/stated.rs"]
mod stated;
use nacre_exact::{Angle, Axis, Isometry, Rat, Rotation};
use nacre_ops::{DatumDef, SketchFrame, SketchPlane};
use nacre_store::Handle;
use nacre_topo::{Model, PointCache, Solid};
use stated::*;

/// State `plane` as a datum and hand back the frame it implies — the two steps a caller takes
/// when the plane is not one the model already holds (a seed, or a face's).
fn datum_frame(m: &mut Model, plane: SketchPlane) -> SketchFrame {
    match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating a plane: {other:?}"),
    }
}

/// A cylinder the way an application states one — a datum plane, a whole circle about its origin,
/// an extrude ([`nacre_ops::fixtures::cylinder_with_seam`]) — with the seam on
/// `axis.any_perpendicular()`.
///
/// ★ **The seam is stated, not left to the frame**, because rows here were built around where it
/// falls (`arc straddle` puts it on the pierce point) and this corpus is frozen by holding its own
/// fixture. The axes stated here are integer vectors or a 3-4-5 triple, so the seam and the frame
/// are exact decimals.
fn cylinder(
    m: &mut Model,
    base: Point3,
    axis: nacre_math::Vector3,
    radius: f64,
    height: f64,
) -> Handle<Solid> {
    let seam = axis.any_perpendicular().expect("a nonzero axis");
    nacre_ops::fixtures::cylinder_with_seam(m, base, axis, seam, radius, height).solid
}

fn xf(m: &mut Model, s: Handle<Solid>, iso: Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: iso,
        },
    )
    .expect("transform") else {
        unreachable!("transform yields Transform output")
    };
    m.rebuild_adjacency();
    solid
}

fn mirror(m: &mut Model, s: Handle<Solid>, axis: Axis, offset: Rat) -> Handle<Solid> {
    let OpOutput::Mirror { solid } = apply(
        m,
        &Operation::Mirror {
            solid: s,
            axis,
            offset,
        },
    )
    .expect("mirror") else {
        unreachable!("mirror yields Mirror output")
    };
    m.rebuild_adjacency();
    solid
}

/// A hash of the solid's **plane caches** — origin and unit normal — sorted by bits, the other
/// half of what a boolean actually reads.
///
/// **Measured, not assumed:** a boolean's result vertices are all recomputed from plane triples
/// (a measured `Vertex::ThreePlane`), so the result carries *no* trace of the operands' vertex
/// coordinates. A change that moves operand vertices but leaves the planes alone is therefore
/// invisible in the result — which is exactly what happened the first time this census was used.
/// Recording the operands is what makes the census see the change it exists to see.
fn plane_digest(m: &Model, s: Handle<Solid>) -> (usize, u64) {
    use std::hash::{Hash, Hasher};
    let mut bits: Vec<[u64; 6]> = Vec::new();
    let src = m.solid(s).clone();
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &m.shell(sh).faces {
            match m.surface_cache(m.face(fh).surface) {
                nacre_geom::Surface::Plane(pl) => {
                    let (o, n) = (pl.origin().as_array(), pl.normal().as_array());
                    bits.push([o[0], o[1], o[2], n[0], n[1], n[2]].map(f64::to_bits));
                }
                nacre_geom::Surface::Cylinder(_) => {}
            }
        }
    }
    bits.sort_unstable();
    bits.dedup();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bits.hash(&mut h);
    (bits.len(), h.finish())
}

/// The operands' own state, before the boolean consumes them: vertices *and* planes.
fn operands(m: &Model, a: Handle<Solid>, b: Handle<Solid>) -> String {
    let f = |s| {
        let (vn, vh) = coord_digest(m, s);
        let (pn, ph) = plane_digest(m, s);
        format!("v{vn}h{vh:016x}p{pn}h{ph:016x}")
    };
    format!("{}+{}", f(a), f(b))
}

/// A hash of the solid's vertex coordinates, **sorted by bits** so it does not depend on traversal
/// or handle order — only on the set of coordinates the operation produced.
fn coord_digest(m: &Model, s: Handle<Solid>) -> (usize, u64) {
    use std::hash::{Hash, Hasher};
    let mut bits: Vec<[u64; 3]> = Vec::new();
    let src = m.solid(s).clone();
    for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
        for &fh in &m.shell(sh).faces {
            let face = m.face(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for &vh in m.edge(he.edge).vertices.iter() {
                        let p = m.vertex_point(vh).as_array();
                        bits.push([p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]);
                    }
                }
            }
        }
    }
    bits.sort_unstable();
    bits.dedup();
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bits.hash(&mut h);
    (bits.len(), h.finish())
}

/// Running totals of [`nacre_ops::audit_plane_senses`] over every row — `agree`, `unmeasured`,
/// `mirrored`, `turned`. This binary's dump is its only test, so a process-global sum is the
/// corpus's.
/// Result-vertex coordinates the realization read as `+0.0` by the coincidence rule
/// (`Realized::to_f64`) — the ones whose cache is `0.0` with a nonzero bound. An exact zero carries
/// a zero bound, and no other decided coordinate is `0.0`, so this counts the rule's answers
/// without a counter in the product.
static ZERO_BY_COINCIDENCE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

static SENSE: [std::sync::atomic::AtomicUsize; 4] =
    [const { std::sync::atomic::AtomicUsize::new(0) }; 4];

/// Circle edges the rim lock reads: `[stated, unstated]` — a stated rim's centre is checked
/// against the truth; an unstated one (a cylinder without a world statement, or a cap without a
/// narrow world name) still takes the caches' `f64` meet, and is counted so the population shows.
static CIRCLES: [std::sync::atomic::AtomicUsize; 2] =
    [const { std::sync::atomic::AtomicUsize::new(0) }; 2];

/// The centre a circle edge's truth names, as the nearest `f64` — `None` where the cylinder has no
/// world statement or its cap no narrow world name.
///
/// ★ **Spelled here, not borrowed.** The axis `o + t·m` meets the cap `n·x + d = 0` at
/// `t = −(n·o + d)/(n·m)`, written out in this test's own `Rat` arithmetic from the truth's public
/// statements, so the kernel's meet is checked against something it did not compute.
fn circle_centre_truth(
    m: &Model,
    cyl: Handle<nacre_topo::Surface>,
    cap: Handle<nacre_topo::Surface>,
) -> Option<[f64; 3]> {
    let def = m.world_cylinder_def(cyl)?;
    let k = *m.world_plane_name(cap)?.narrow()?;
    let (o, d) = (def.origin(), def.dir());
    let dot = |v: &[Rat; 3]| -> Option<Rat> {
        k[0].checked_mul(v[0])?
            .checked_add(k[1].checked_mul(v[1])?)?
            .checked_add(k[2].checked_mul(v[2])?)
    };
    let nm = dot(&d)?;
    let t = Rat::from_int(0)
        .checked_sub(dot(&o)?.checked_add(k[3])?)?
        .checked_mul(Rat::new(nm.denom(), nm.numer())?)?;
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = o[i].checked_add(t.checked_mul(d[i])?)?.to_f64();
    }
    Some(c)
}

/// Every live circle edge of `m` against its truth: the tags of the stated circles whose cached
/// centre is not the truth's nearest `f64`, and of the circles whose frame (normal, `ref_dir`,
/// radius) is not the cylinder cache's to the bit.
fn circle_audit(m: &Model) -> (Vec<String>, Vec<String>) {
    let (mut centre, mut frame) = (Vec::new(), Vec::new());
    let mut edges: Vec<_> = m.reachable().edges.into_iter().collect();
    edges.sort_by_key(|e| e.index());
    for eh in edges {
        let nacre_geom::Curve::Circle(c) = m.edge_curve(eh) else {
            continue;
        };
        let [s0, s1] = m.edge(eh).surfaces;
        let (cyl, cap) = match m.surface(s0) {
            nacre_topo::Surface::Cylinder { .. } => (s0, s1),
            nacre_topo::Surface::Plane { .. } => (s1, s0),
        };
        let nacre_geom::Surface::Cylinder(cc) = m.surface_cache(cyl) else {
            unreachable!("a circle edge has a cylinder carrier")
        };
        let bits = |v: nacre_math::Vector3| v.as_array().map(f64::to_bits);
        if bits(c.normal()) != bits(cc.axis().direction())
            || bits(c.ref_dir()) != bits(cc.ref_dir())
            || c.radius().to_bits() != cc.radius().to_bits()
        {
            frame.push(format!("e{}", eh.index()));
        }
        match circle_centre_truth(m, cyl, cap) {
            Some(t) => {
                CIRCLES[0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if t.map(f64::to_bits) != c.center().as_array().map(f64::to_bits) {
                    centre.push(format!(
                        "e{} {:?} vs {:?}",
                        eh.index(),
                        c.center().as_array(),
                        t
                    ));
                }
            }
            None => {
                CIRCLES[1].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }
    (centre, frame)
}

fn record(
    tag: &str,
    m: &Model,
    inputs: &str,
    out: &Result<Vec<Handle<Solid>>, nacre_ops::BoolError>,
) {
    // ★ The rim lock: every circle edge's frame is its cylinder cache's, and every stated rim's
    // centre is the truth's nearest `f64` — the edge cache's twin of the vertex lock below.
    let (centre, frame) = circle_audit(m);
    assert!(
        centre.is_empty(),
        "{tag}: rim centres that are not the truth's nearest f64: {centre:?}"
    );
    assert!(
        frame.is_empty(),
        "{tag}: rim frames that are not the cylinder cache's: {frame:?}"
    );
    // ★ The plane-sense lock: every live plane's cache faces the way its truth's sense says.
    let audit = nacre_ops::audit_plane_senses(m);
    assert!(
        audit.disagree.is_empty(),
        "{tag}: planes whose cache opposes their stated sense: {:?}",
        audit.disagree
    );
    for (slot, n) in SENSE
        .iter()
        .zip([audit.agree, audit.unmeasured, audit.mirrored, audit.turned])
    {
        slot.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
    }
    // ★ The mesh lock: every solid a boolean returns here can be drawn. The lib suite asserts the
    // same at the boolean's exit (`tess_census`, `cfg(test)`), which this corpus does not pass through.
    if let Ok(v) = out
        && !v.is_empty()
    {
        let mesh = nacre_tess::tessellate(m, &nacre_tess::TessConfig::default());
        assert!(
            mesh.is_ok(),
            "{tag}: a boolean built a solid the mesher refuses: {:?}",
            mesh.err()
        );
    }
    print!("c {tag} in:{inputs} ");
    match out {
        Err(e) => println!("ERR {e:?}"),
        Ok(v) if v.is_empty() => println!("EMPTY"),
        Ok(v) => {
            // Every body, in the order the engine returned them — the count is part of the answer.
            let mut parts = Vec::new();
            for &s in v {
                let p = nacre_props::mass_props(m, s).expect("props");
                // ★ **A curved result has no centroid to print, and that is a record, not a
                // crash.** `centroid` refuses a face it cannot integrate a first moment over
                // (`CentroidOfCurvedFace`), which the census met the moment cylinder booleans
                // started returning solids. Volume and area still measure the whole body, so the
                // row keeps its comparable numbers and marks the missing field rather than
                // taking the process down.
                let (n, d) = coord_digest(m, s);
                let centroid = match nacre_props::centroid(m, s) {
                    Ok(c) => format!(
                        "{:016x},{:016x},{:016x}",
                        c[0].to_bits(),
                        c[1].to_bits(),
                        c[2].to_bits()
                    ),
                    Err(nacre_props::PropsError::CentroidOfCurvedFace) => "curved".to_string(),
                    Err(e) => panic!("centroid: {e:?}"),
                };
                parts.push(format!(
                    "{:016x}/{:016x}/{centroid}/v{n}h{d:016x}",
                    p.volume.to_bits(),
                    p.area.to_bits(),
                ));
            }
            println!("{} {}", v.len(), parts.join(" "));
            // ★ **The cache is the realization**: every result vertex the push funnel
            // could realize on the ladder's first rung is `Bounded` with that value bit for bit.
            // Asked again here with the funnel's own question, so the row cannot pass by measuring
            // nothing.
            //
            // ★★ And the refusals are counted **by name**, not as one lump. "The cache road stopped
            // and a paid road can answer" (`ceiling`) and "there is no road" (`unrealized`) are
            // different reports — one is a cost, the other is a gap in the kernel — and a single
            // `kept=` would hide the second behind the first as the population moves.
            let mut seen = std::collections::HashSet::new();
            let (mut realized, mut kept, mut ceiling) = (0usize, 0usize, 0usize);
            for &s in v {
                let src = m.solid(s).clone();
                for &sh in std::iter::once(&src.outer).chain(src.cavities.iter()) {
                    for &fh in &m.shell(sh).faces {
                        let face = m.face(fh);
                        for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                            for he in &lp.half_edges {
                                for &vh in m.edge(he.edge).vertices.iter() {
                                    if !seen.insert(vh) {
                                        continue;
                                    }
                                    // ★★ The contract is now **answer ⇔ variant**, three ways: the
                                    // road answers and the cache is `Bounded` with that very
                                    // coordinate; the road stops for bits or cost and the cache is
                                    // `Ceiling`; the road has no way there at all and the cache is
                                    // `Unrealized`. Written out rather than as "not Bounded", so a
                                    // vertex cannot drift between the two refusals unnoticed.
                                    match nacre_ops::realize_cache(m, m.vertex(vh)) {
                                        Ok((c, b)) => {
                                            realized += 1;
                                            let zeros = (0..3)
                                                .filter(|&k| c[k] == 0.0 && !b[k].is_zero())
                                                .count();
                                            ZERO_BY_COINCIDENCE.fetch_add(
                                                zeros,
                                                std::sync::atomic::Ordering::Relaxed,
                                            );
                                            assert!(
                                                matches!(m.vertex_cache(vh), PointCache::Bounded { coord, .. } if coord.as_array() == c),
                                                "{tag}: vertex {} is not the realization: {:?} vs {c:?}",
                                                vh.index(),
                                                m.vertex_cache(vh)
                                            );
                                        }
                                        Err(d) => {
                                            kept += 1;
                                            let want_ceiling = matches!(
                                                d,
                                                nacre_ops::CacheDecline::CostCap
                                                    | nacre_ops::CacheDecline::Cannot(
                                                        nacre_ops::RealizeError::Undecided
                                                    )
                                            );
                                            let is_ceiling = matches!(
                                                m.vertex_cache(vh),
                                                PointCache::Ceiling { .. }
                                            );
                                            ceiling += usize::from(is_ceiling);
                                            assert_eq!(
                                                want_ceiling,
                                                is_ceiling,
                                                "{tag}: vertex {} declined with {d:?} but the cache says {:?}",
                                                vh.index(),
                                                m.vertex_cache(vh)
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            println!(
                "r {tag} realized={realized} ceiling={ceiling} unrealized={}",
                kept - ceiling
            );
        }
    }
}

const KINDS: [(&str, BoolKind); 3] = [
    ("fuse", BoolKind::Fuse),
    ("cut", BoolKind::Cut),
    ("common", BoolKind::Common),
];

#[test]
#[ignore = "census dump with plane-sense and mesh locks; the gates run it (--ignored --nocapture)"]
fn measure_census() {
    // ── Axis-aligned, no motion at all: the untouched baseline.
    let boxes: [([f64; 3], [f64; 3]); 6] = [
        ([0.5, 0.5, 0.5], [2.5, 2.5, 2.5]),
        ([0.25, 0.25, 0.25], [0.75, 0.75, 0.75]),
        ([1.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
        ([0.3, -1.0, 0.3], [0.7, 2.0, 0.7]),
        ([-1.0, 0.2, 0.2], [0.5, 0.8, 0.8]),
        ([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]),
    ];
    for (kn, k) in KINDS {
        for (i, (lo, hi)) in boxes.iter().enumerate() {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([1.0; 3]),
            );
            let b = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array(*lo),
                Point3::from_array(*hi),
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("aa {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── Rotated: one operand turned about each axis by several angles.
    for ax in [Axis::X, Axis::Y, Axis::Z] {
        for deg in [7i128, 30, 45, 90, 123] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([2.0, 2.0, 2.0]),
                );
                let b = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.4, 0.4, -1.0]),
                    Point3::from_array([1.6, 1.6, 3.0]),
                );
                m.rebuild_adjacency();
                let b = xf(
                    &mut m,
                    b,
                    Isometry::rotation(Rotation {
                        axis: ax,
                        pivot: [Rat::from_int(1); 3],
                        angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
                    }),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("rot {ax:?} {deg} {kn}"), &m, &inputs, &out);
            }
        }
    }

    // ── **Translated.**
    // Dyadic offsets and non-dyadic ones (whose `f64` images round), on both operands so the
    // *relative* placement varies too. Every one is carried into the statements.
    let offsets: [(i128, i128); 10] = [
        (1, 1),
        (3, 1),
        (1, 2),
        (-5, 4),
        (7, 8),
        (1, 3),
        (2, 5),
        (3, 10),
        (7, 11),
        (13, 23),
    ];
    for (n, d) in offsets {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0; 3]),
            );
            let b = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([1.0, 1.0, 1.0]),
                Point3::from_array([3.0, 3.0, 3.0]),
            );
            m.rebuild_adjacency();
            let t = Rat::new(n, d).expect("offset");
            let a = xf(
                &mut m,
                a,
                Isometry::translation([t, Rat::from_int(0), Rat::from_int(0)]),
            );
            let b = xf(&mut m, b, Isometry::translation([t, t, Rat::from_int(0)]));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("tr {n}/{d} {kn}"), &m, &inputs, &out);
        }
        // Chained: the same total move split in two, which f64 accumulation and exact folding
        // disagree about.
        let mut m = Model::new();
        let a = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0; 3]),
        );
        let b = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([1.0, 1.0, 1.0]),
            Point3::from_array([3.0, 3.0, 3.0]),
        );
        m.rebuild_adjacency();
        let t = Rat::new(n, d).expect("offset");
        let mut a = a;
        for _ in 0..2 {
            a = xf(
                &mut m,
                a,
                Isometry::translation([t, Rat::from_int(0), Rat::from_int(0)]),
            );
        }
        let b = xf(
            &mut m,
            b,
            Isometry::translation([
                t.checked_add(t).expect("double"),
                Rat::from_int(0),
                Rat::from_int(0),
            ]),
        );
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, BoolKind::Fuse, a, b);
        m.rebuild_adjacency();
        record(&format!("tr2 {n}/{d} fuse"), &m, &inputs, &out);
    }

    // ── Mirrored, about dyadic and non-dyadic planes — every one carried into the statements.
    for (n, d) in [(0i128, 1i128), (1, 2), (1, 3), (7, 22), (5, 7)] {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([1.0, 0.0, 0.0]),
                Point3::from_array([2.0, 1.0, 1.0]),
            );
            let b = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([-1.0, 0.25, 0.25]),
                Point3::from_array([1.5, 0.75, 0.75]),
            );
            m.rebuild_adjacency();
            let a = mirror(&mut m, a, Axis::X, Rat::new(n, d).expect("plane"));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("mir {n}/{d} {kn}"), &m, &inputs, &out);
        }
    }

    // ── Rotate then place.
    for deg in [7i128, 30, 45] {
        let mut m = Model::new();
        let base = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([-2.0, -2.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        let tool = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([-0.5, -0.5, -1.0]),
            Point3::from_array([0.5, 0.5, 2.0]),
        );
        m.rebuild_adjacency();
        let tool = xf(
            &mut m,
            tool,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
            }),
        );
        let place = Isometry::translation([Rat::from_int(5), Rat::from_int(-3), Rat::from_int(2)]);
        let tool = xf(&mut m, tool, place);
        let base = xf(&mut m, base, place);
        let inputs = operands(&m, base, tool);
        let out = boolean(&mut m, BoolKind::Cut, base, tool);
        m.rebuild_adjacency();
        record(&format!("rt {deg}"), &m, &inputs, &out);
    }

    // ── **Constructed by `extrude`, with the dimension split across steps.**
    //
    // Every family above builds its boxes with `fixtures::cuboid` — one floor and one height each,
    // so none of them can see a change to how construction *accumulates* over steps.
    //
    // Here one operand is raised in a single step and the other in two that should add to the
    // same height. In `f64` they do not (`1.1 + 6.6 != 7.7`), and the difference lands in the
    // operand digests below — which is the point of the family.
    for (kn, k) in KINDS {
        for (i, (whole, first, second)) in [(7.7, 1.1, 6.6), (3.0, 0.1, 2.9), (1.0, 0.3, 0.7)]
            .iter()
            .enumerate()
        {
            let mut m = Model::new();
            let a = ex(&mut m, 0.0, [0.0, 0.0], [2.0, 2.0], *whole);
            let b1 = ex(&mut m, 0.0, [1.0, 1.0], [3.0, 3.0], *first);
            let b2 = ex(&mut m, *first, [1.0, 1.0], [3.0, 3.0], *second);
            let b = boolean(&mut m, BoolKind::Fuse, b1, b2).expect("stack fuses")[0];
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("ex {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── ★★★★★ **Full-width coordinates — the population every family above is blind to.**
    //
    // Every corner and dimension above is a *short* decimal (`0.5`, `2.5`, `7.7`), and every
    // rational derived from one stays far inside `i128`. So a change to what the kernel can and
    // cannot **name** leaves all 130 lines bit-identical, and this census reports nothing.
    //
    // ★★★★★ That is not hypothetical: a regression that closed the exact path for every
    // arbitrarily-tilted plane shipped, and the census was bit-identical across it — **twice**,
    // because the same blindness had already been recorded once.
    //
    // ★★ **Measured, and it picks the right target.** A cuboid on seventeen-digit corners
    // overflows `plane_through_points`' *intermediates* on two of its six faces — the
    // narrow-route population the wide derivation (`plane_name_big`) names, since their
    // canonical answers are small. (The genuinely wide-vessel population lives in computed ring
    // coordinates and is counted by `WIDE_PLANES`.) The proptests are where this population lives
    // (`stacked_boxes_merge_volumes` and friends generate arbitrary `f64`), and a proptest
    // cannot be a census line: it has no fixed coordinates to diff. These constants are those
    // coordinates, pinned.
    let fw: [([f64; 3], [f64; 3]); 3] = [
        (
            [-2.8374652839472, 1.0937465283947, -0.5837465283947],
            [1.4738264859372, 2.9384756293847, 0.8473625849372],
        ),
        // A stacked pair sharing one interface plane — the coplanar-contact route, on coordinates
        // whose interface plane the narrow route cannot name.
        (
            [-2.8374652839472, 1.0937465283947, 0.8473625849372],
            [1.4738264859372, 2.9384756293847, 2.1937465283947],
        ),
        // Overlapping, so the result's planes come from both operands.
        (
            [-1.1937465283947, 1.9384756293847, -0.1837465283947],
            [2.8473625849372, 3.4738264859372, 1.4937465283947],
        ),
    ];
    for (kn, k) in KINDS {
        for (i, (lo, hi)) in fw.iter().enumerate() {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array(fw[0].0),
                Point3::from_array(fw[0].1),
            );
            let b = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array(*lo),
                Point3::from_array(*hi),
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("fw {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── **A prism on a tilted sketch plane** — the `Motion::Frame` path, which no family above
    // reaches either. Its ring lives in the plane's own frame, so its walls *are* nameable; what
    // it exercises is the frame machinery, not the width limit.
    for (i, n) in [
        [0.3141592653589793, -0.2718281828459045, 1.0],
        [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
        [0.1, 0.2, 0.30000000000000004],
    ]
    .iter()
    .enumerate()
    {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = tilted_prism(&mut m, *n, 0.0, 1.0, 0.75);
            // A second prism on the *same* tilted plane, offset in the sketch and raised further:
            // its base cap and `a`'s are one plane, and its walls can be coplanar with `a`'s.
            let b = tilted_prism(&mut m, *n, 0.5, 1.6, 1.25);
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("tp {kn} {i}"), &m, &inputs, &out);
        }
    }

    // ── **A framed sketch on an `n·n`-overflow wall** (`wf`) — a population which
    // no family above contains: `fw` is all-narrow (measured) and `tp`'s walls are
    // in-frame narrow. A prism on a fully tilted, exactly-orthonormal *decimal* frame has walls
    // whose names run ~110 bits — narrow, with squared lengths past `i128` — and a pad or
    // pocket on such a wall must not fall to the f64 path. Without this family a regression in
    // the wide-frame road would cross this file bit-identical.
    {
        use nacre_math::{Point2, Vector3};
        let p2 = |x: f64, y: f64| Point2::from_array([x, y]);
        for (kn, pad) in [("pocket", false), ("pad", true)] {
            let mut m = Model::new();
            let plane = nacre_ops::SketchPlane::from_axes(
                Point3::from_array([0.1234567890123456, 0.2345678901234567, 0.3456789012345678]),
                Vector3::from_array([0.6, 0.8, 0.0]),
                Vector3::from_array([-0.48, 0.36, 0.8]),
            );
            let __f105 = datum_frame(&mut m, plane);
            let OpOutput::Extrude { solid, .. } = apply(
                &mut m,
                &Operation::Extrude {
                    frame: __f105,
                    profile: nacre_ops::Profile2d::polygon(vec![
                        p2(0.1111111111111111, 0.1234567890123456),
                        p2(4.123456789012345, 0.2345678901234567),
                        p2(3.9876543210987654, 3.1234567890123459),
                        p2(0.2222222222222222, 2.765432109876543),
                    ])
                    .unwrap(),
                    dist: 2.5,
                },
            )
            .expect("the wf base prism") else {
                unreachable!()
            };
            m.rebuild_adjacency();
            // The family qualifies itself: a wall whose name is narrow but whose squared
            // lengths overflow — the exact population the narrow frame derivation declines.
            let wall = *m
                .shell(m.solid(solid).outer)
                .faces
                .iter()
                .find(|&&f| {
                    let s = m.face(f).surface;
                    m.surface_name
                        .get(&s)
                        .and_then(|n| n.narrow())
                        .is_some_and(|c| nacre_exact::plane_frame_default(*c).is_none())
                })
                .expect("the wf population vanished — retune the constants");
            // A small square centred on the wall, in its own sketch frame.
            let sp = nacre_ops::face_plane(&m, wall).expect("planar");
            let d = nacre_props::face_props(&m, wall).unwrap().centroid - sp.origin();
            let (cu, cv) = (d.dot(sp.x_axis()), d.dot(sp.y_axis()));
            let profile = nacre_ops::Profile2d::polygon(vec![
                p2(cu - 0.3, cv - 0.3),
                p2(cu + 0.3, cv - 0.3),
                p2(cu + 0.3, cv + 0.3),
                p2(cu - 0.3, cv + 0.3),
            ])
            .unwrap();
            let inputs = operands(&m, solid, solid);
            // The feature as an application builds it — the wall's frame, the tool swept off it
            // (into the wall for a pocket), the boolean with the wall's solid first. Written out
            // here rather than through a shared fixture, so the corpus does not move with one.
            let frame = nacre_ops::face_sketch_frame(&m, wall).expect("the wall's frame");
            let Ok(OpOutput::Extrude { solid: tool, .. }) = apply(
                &mut m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist: if pad { 0.4 } else { -0.4 },
                },
            ) else {
                panic!("the wf tool must build")
            };
            let kind = if pad { BoolKind::Fuse } else { BoolKind::Cut };
            let Ok(OpOutput::Boolean { solids }) = apply(
                &mut m,
                &Operation::Boolean {
                    kind,
                    a: solid,
                    b: tool,
                },
            ) else {
                panic!("the wf feature must build")
            };
            m.rebuild_adjacency();
            record(&format!("wf {kn}"), &m, &inputs, &Ok(solids));
        }
    }

    // ── **Contact — two bodies that meet along a line or at a point, and one pair that cannot part.**
    //
    // ★ The population this file was blind to. Before the contact work every one of these was a
    // reject, and census carried **no reject line at all** (measured: 152 lines, 7 `EMPTY`, zero
    // `ERR`), so the whole family was invisible to the gate that reads every commit. It records
    // *both* halves on purpose: the separable contacts coming back as the bodies they are, and the
    // pinch that stays a reject because the material loops around it.
    for (kn, k) in KINDS {
        // Two cubes sharing exactly the vertical line x=1, y=1.
        let mut m = Model::new();
        let a = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0; 3]),
        );
        let b = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([1.0, 1.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        m.rebuild_adjacency();
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, k, a, b);
        m.rebuild_adjacency();
        record(&format!("ct edge {kn}"), &m, &inputs, &out);

        // Two cubes sharing exactly the point (1,1,1).
        let mut m = Model::new();
        let a = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([1.0; 3]),
        );
        let b = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([1.0; 3]),
            Point3::from_array([2.0; 3]),
        );
        m.rebuild_adjacency();
        let inputs = operands(&m, a, b);
        let out = boolean(&mut m, k, a, b);
        m.rebuild_adjacency();
        record(&format!("ct corner {kn}"), &m, &inputs, &out);
    }
    // ★★ **The pinch that must stay a reject.** A and B meet only along the line x=2, y=2, but a
    // bridge overlaps both, so the material loops around the contact: cut it there and one piece
    // remains, not two. A line that ever stops saying `ERR` here is the separation going too far.
    {
        let mut m = Model::new();
        let a = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0, 0.0, 0.0]),
            Point3::from_array([2.0, 2.0, 1.0]),
        );
        let bridge = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([1.0, 0.3, 0.2]),
            Point3::from_array([3.0, 3.0, 0.8]),
        );
        let b = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([2.0, 2.0, 0.0]),
            Point3::from_array([4.0, 4.0, 1.0]),
        );
        m.rebuild_adjacency();
        let ab = boolean(&mut m, BoolKind::Fuse, a, bridge).expect("a and its bridge overlap");
        m.rebuild_adjacency();
        let inputs = operands(&m, ab[0], b);
        let out = boolean(&mut m, BoolKind::Fuse, ab[0], b);
        m.rebuild_adjacency();
        record("ct ring-pinch fuse", &m, &inputs, &out);
    }
    // ── Cylinders: the truth rides beside the cache. Solo and moved bodies are digest
    // lines (a boolean never runs); the boolean rows record the answer, refusals included — a
    // reject string is part of the corpus, so a change in what is admitted moves these lines
    // between ERR and results *in the diff*, not silently.
    {
        use nacre_math::Vector3;
        let solo = |m: &Model, s: Handle<Solid>| {
            let (vn, vh) = coord_digest(m, s);
            let (pn, ph) = plane_digest(m, s);
            format!("v{vn}h{vh:016x}p{pn}h{ph:016x}")
        };
        let mut m = Model::new();
        let c = cylinder(
            &mut m,
            Point3::from_array([0.5, -1.25, 2.0]),
            nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
            1.5,
            2.5,
        );
        m.rebuild_adjacency();
        println!("c cyl solo {}", solo(&m, c));
        // An exact 90° turn (truth transported, nothing recorded) and an inexact 31° turn
        // (node recorded, def carried verbatim) — both roads in the corpus.
        let turned = xf(
            &mut m,
            c,
            Isometry::rotation(Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(90)).expect("angle"),
            }),
        );
        println!("c cyl turn90 {}", solo(&m, turned));
        let leaned = xf(
            &mut m,
            turned,
            Isometry::rotation(Rotation {
                axis: Axis::X,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(31)).expect("angle"),
            }),
        );
        println!("c cyl turn31 {}", solo(&m, leaned));

        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([2.0, 2.0, 2.0]),
            );
            let b = cylinder(
                &mut m,
                Point3::from_array([1.0, 1.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("cyl {kn}"), &m, &inputs, &out);
        }
        // ── **Cut rims**: a boss whose circle a boundary segment cuts — the arc
        // population. Three placements: straddling the plate's top edge (seam ≡ pierce), turned
        // over the corner (the seam splits the wrap arc), and hung under the bottom edge (the
        // cut circle is the band's hi end).
        for (pn, origin, axis) in [
            ("straddle", [4.0, 2.0, 2.0], [0.0, 0.0, 1.0]),
            ("turned", [4.0, 0.25, 2.0], [1.0, 0.0, 0.0]),
            ("hung", [4.0, 2.0, -1.0], [0.0, 0.0, 1.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([4.0, 4.0, 2.0]),
                );
                let b = cylinder(
                    &mut m,
                    Point3::from_array(origin),
                    Vector3::from_array(axis),
                    0.5,
                    1.0,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arc {pn} {kn}"), &m, &inputs, &out);
            }
        }
        // ── **A bored plate, then a straddling boss** (chaining wall 3): the second boolean's
        // nesting reads a bitten top ring — pierce corners and an arc step — so the containment
        // parity runs the mixed road instead of the rational chart. The first boolean is fixed
        // (`cut` the through-bore); the second varies by kind.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let bore = cylinder(
                &mut m,
                Point3::from_array([1.0, 1.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            let bored = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            let boss = cylinder(
                &mut m,
                Point3::from_array([4.0, 2.0, 2.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                1.0,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, bored, boss);
            let out = boolean(&mut m, k, bored, boss);
            m.rebuild_adjacency();
            record(&format!("arc bored straddle {kn}"), &m, &inputs, &out);
        }
        // ── **The rulings road** (the gate's record-and-pass arm): a boss whose axis
        // lies exactly on the plate's wall plane builds (through, the corner's two walls, an
        // asymmetric station); offset (`0 < d < r`)
        // and tangent (`d = r`) both pass the gate too, and the tangent row's
        // three kinds are where the *verdict* answers (`Fuse` builds, `Cut` pinches).
        for (pn, base, h) in [
            ("through", [40.0, 20.0, -10.0], 50.0),
            ("corner", [40.0, 40.0, -10.0], 50.0),
            ("offmid", [40.0, 10.0, -10.0], 50.0),
            ("offset", [38.0, 20.0, -10.0], 50.0),
            ("tangent", [35.0, 20.0, -10.0], 50.0),
            ("half", [40.0, 20.0, -10.0], 20.0),
            ("flush", [40.0, 20.0, 0.0], 20.0),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([40.0, 40.0, 20.0]),
                );
                let b = cylinder(
                    &mut m,
                    Point3::from_array(base),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    5.0,
                    h,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("rul {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **Face-to-face contact** (`ct2`): two cells meeting on a full wall, which is what a
    // pattern of parts is made of. The shared wall is interior, so the corners on it must
    // dissolve — and they only can once each cell's caps merge, which is what the coplanar
    // merge's circle carry opened. Bored, pocketed+bored, and a three-cell chain (the second
    // fuse meets a result that already carries merged caps).
    {
        let cell = |m: &mut Model, x0: f64, pocket: bool| {
            let plate = nacre_ops::fixtures::cuboid(
                m,
                Point3::from_array([x0, 0.0, 0.0]),
                Point3::from_array([x0 + 20.0, 20.0, 10.0]),
            );
            m.rebuild_adjacency();
            let mut out = plate;
            if pocket {
                let p = nacre_ops::fixtures::cuboid(
                    m,
                    Point3::from_array([x0 + 2.0, 2.0, 4.0]),
                    Point3::from_array([x0 + 8.0, 8.0, 10.0]),
                );
                m.rebuild_adjacency();
                out = boolean(m, BoolKind::Cut, out, p).expect("the pocket cuts")[0];
                m.rebuild_adjacency();
            }
            let bore = cylinder(
                m,
                Point3::from_array([x0 + 14.0, 14.0, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.0,
                12.0,
            );
            m.rebuild_adjacency();
            out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            out
        };
        for (pn, pocket) in [("bored", false), ("pocketed", true)] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = cell(&mut m, 0.0, pocket);
                let b = cell(&mut m, 20.0, pocket);
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("ct2 {pn} {kn}"), &m, &inputs, &out);
            }
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = cell(&mut m, 0.0, false);
            let b = cell(&mut m, 20.0, false);
            let ab = boolean(&mut m, BoolKind::Fuse, a, b).expect("the first pair fuses")[0];
            m.rebuild_adjacency();
            let c = cell(&mut m, 40.0, false);
            let inputs = operands(&m, ab, c);
            let out = boolean(&mut m, k, ab, c);
            m.rebuild_adjacency();
            record(&format!("ct2 chain {kn}"), &m, &inputs, &out);
        }
    }
    // ── **Translated cylinders** (`trc`): cylinders moved by a non-dyadic offset, carried into
    // their statements — world cylinders where they land. Four placements: a
    // tool cutting, a boss fusing, a bored body fused onto its twin (the plane side), and a
    // bore translated **onto** another bore — whose axes then coincide, which the cylinder-pair
    // rule refuses by name (`CylinderPairContact`); before the road opened it never got that
    // far and said `CylinderGateUndecided` instead.
    {
        let off = |x: f64, y: f64, z: f64| {
            Isometry::translation([x, y, z].map(|c| Rat::try_from_f64(c).expect("rational")))
        };
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let tool = cylinder(
                &mut m,
                Point3::from_array([7.3, 7.3, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let tool = xf(&mut m, tool, off(10.7, 0.0, 0.0));
            let inputs = operands(&m, plate, tool);
            let out = boolean(&mut m, k, plate, tool);
            m.rebuild_adjacency();
            record(&format!("trc tool {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let boss = cylinder(
                &mut m,
                Point3::from_array([10.0, 20.0, 5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                3.0,
                8.0,
            );
            m.rebuild_adjacency();
            let boss = xf(&mut m, boss, off(15.3, 0.1, 0.0));
            let inputs = operands(&m, plate, boss);
            let out = boolean(&mut m, k, plate, boss);
            m.rebuild_adjacency();
            record(&format!("trc boss {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let bored = |m: &mut Model| {
                let plate = nacre_ops::fixtures::cuboid(
                    m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([40.0, 40.0, 10.0]),
                );
                let bore = cylinder(
                    m,
                    Point3::from_array([7.3, 7.3, -5.0]),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    2.1,
                    30.0,
                );
                m.rebuild_adjacency();
                let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
                m.rebuild_adjacency();
                out
            };
            let a = bored(&mut m);
            let b = bored(&mut m);
            let b = xf(&mut m, b, off(25.7, 3.3, 2.0));
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("trc body {kn}"), &m, &inputs, &out);
        }
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let plate = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([40.0, 40.0, 10.0]),
            );
            let bore = cylinder(
                &mut m,
                Point3::from_array([18.0, 20.0, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let holed = boolean(&mut m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            let twin = cylinder(
                &mut m,
                Point3::from_array([7.3, 20.0, -5.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                30.0,
            );
            m.rebuild_adjacency();
            let twin = xf(&mut m, twin, off(10.7, 0.0, 0.0)); // lands on the bore's own axis
            let inputs = operands(&m, holed, twin);
            let out = boolean(&mut m, k, holed, twin);
            m.rebuild_adjacency();
            record(&format!("trc onaxis {kn}"), &m, &inputs, &out);
        }
    }
    // ── **The same boss, on each of the four walls** (`wal`). A part does not care which side of
    // itself a boss sits on, and neither should the kernel — but the rulings road was built and
    // measured on the `+x` wall alone, and the corpus inherited that: every `rul` row above puts
    // the axis on a max-side wall. These rows put the same solid on all four, so the table can
    // see a rule that holds on one side and not the other.
    //
    // The box twins are the control: the same straddle with a cuboid tool builds on every wall,
    // volume 35 exactly, so a row that breaks here needs the **cylinder's** edges. ★ Which *kind*
    // of cylinder edge, this table cannot say — a population fact is not a code fact.
    // Fuse alone for the twins: the control only has to say that the planar straddle builds.
    {
        let plate = |m: &mut Model| {
            let a = nacre_ops::fixtures::cuboid(
                m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            m.rebuild_adjacency();
            a
        };
        for (pn, base) in [
            ("xlo", [0.0, 2.0, -1.0]),
            ("xhi", [4.0, 2.0, -1.0]),
            ("ylo", [2.0, 0.0, -1.0]),
            ("yhi", [2.0, 4.0, -1.0]),
            ("corner-lo", [0.0, 0.0, -1.0]),
            ("corner-hi", [4.0, 4.0, -1.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = plate(&mut m);
                let b = cylinder(
                    &mut m,
                    Point3::from_array(base),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    0.5,
                    4.0,
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("wal {pn} {kn}"), &m, &inputs, &out);
            }
        }
        for (pn, lo, hi) in [
            ("box-xlo", [-0.5, 1.5, -1.0], [0.5, 2.5, 3.0]),
            ("box-xhi", [3.5, 1.5, -1.0], [4.5, 2.5, 3.0]),
            ("box-ylo", [1.5, -0.5, -1.0], [2.5, 0.5, 3.0]),
            ("box-yhi", [1.5, 3.5, -1.0], [2.5, 4.5, 3.0]),
        ] {
            let mut m = Model::new();
            let a = plate(&mut m);
            let b =
                nacre_ops::fixtures::cuboid(&mut m, Point3::from_array(lo), Point3::from_array(hi));
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, BoolKind::Fuse, a, b);
            m.rebuild_adjacency();
            record(&format!("wal {pn} fuse"), &m, &inputs, &out);
        }
    }
    // ── **The boss corpus under rigid motion** (`mot`): the production-side rows of the commuting
    // oracle (`tests.rs::the_boolean_commutes_with_rigid_motion`), one per sign class the oracle
    // names — ∥ wall classes with `frame_sign = −1` (a max-side wall put on a seed plane by a
    // translation; a wall turned onto one by rz90), ⊥ classes with `(axis_up, frame)`
    // `= (true, −1)` (corner-lo with the axis turned to −y) and `(false, −1)` (the top cap put on
    // z = 0), and the offset boss under a rigid motion and a non-dyadic translation, whose `f64`
    // images round while the statements move exactly. ★ Rows were added only where **both**
    // profiles dump them: the transport row (`offset-out` under `rz90 + t(5,−3,2)`) rests on the
    // transport law — without it, `world_cylinder_def`'s postcondition (a `debug_assert`) takes the
    // dev census down there while release answers «disjoint» (Fuse 2 bodies, Cut the plate
    // untouched, Common empty): the silent wrong answer the law closes.
    {
        let boss = |m: &mut Model, base: [f64; 3]| -> (Handle<Solid>, Handle<Solid>) {
            let plate = nacre_ops::fixtures::cuboid(
                m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let b = cylinder(
                m,
                Point3::from_array(base),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                4.0,
            );
            m.rebuild_adjacency();
            (plate, b)
        };
        let rot = |ax: Axis, deg: i128| {
            Isometry::rotation(Rotation {
                axis: ax,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).expect("angle"),
            })
        };
        let t = |x: i128, y: i128, z: i128| {
            Isometry::translation([Rat::from_int(x), Rat::from_int(y), Rat::from_int(z)])
        };
        let rows: Vec<(&str, [f64; 3], Isometry)> = vec![
            ("par wall+x t", [4.0, 2.0, -1.0], t(-4, -4, -2)),
            ("par wall-y rz90", [2.0, 0.0, -1.0], rot(Axis::Z, 90)),
            ("perp corner-lo rx90", [0.0, 0.0, -1.0], rot(Axis::X, 90)),
            ("perp wall+y t", [2.0, 4.0, -1.0], t(-4, -4, -2)),
            (
                "xport offset-out rz90+t",
                [4.3, 2.0, -1.0],
                Isometry::rigid(
                    Rotation {
                        axis: Axis::Z,
                        pivot: [Rat::from_int(0); 3],
                        angle: Angle::from_deg(Rat::from_int(90)).expect("angle"),
                    },
                    [Rat::from_int(5), Rat::from_int(-3), Rat::from_int(2)],
                ),
            ),
            (
                "rec through t",
                [2.0, 2.0, -1.0],
                Isometry::translation([
                    Rat::new(7, 11).expect("rat"),
                    Rat::new(3, 10).expect("rat"),
                    Rat::new(1, 4).expect("rat"),
                ]),
            ),
        ];
        for (name, base, iso) in &rows {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (p, b) = boss(&mut m, *base);
                let p = xf(&mut m, p, *iso);
                let b = xf(&mut m, b, *iso);
                let inputs = operands(&m, p, b);
                let out = boolean(&mut m, k, p, b);
                m.rebuild_adjacency();
                record(&format!("mot {name} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A translation, then a quarter turn** (`mot h+…`): the boss corpus again, moved by
    // `(1/10, 1/10, 1/10)` first. Both motions are carried into the statements, so these rows
    // leave no node; a recorded quarter turn is the commuting oracle's `recorded block` family.
    {
        for (name, axis) in [("h+rx90", Axis::X), ("h+rz90", Axis::Z)] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let plate = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([4.0, 4.0, 2.0]),
                );
                let boss = cylinder(
                    &mut m,
                    Point3::from_array([2.0, 2.0, -1.0]),
                    nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                    0.5,
                    4.0,
                );
                m.rebuild_adjacency();
                let prefix = Isometry::translation([Rat::new(1, 10).expect("1/10"); 3]);
                let turn = Isometry::rotation(Rotation {
                    axis,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(90)).expect("angle"),
                });
                let plate = xf(&mut m, plate, prefix);
                let boss = xf(&mut m, boss, prefix);
                let plate = xf(&mut m, plate, turn);
                let boss = xf(&mut m, boss, turn);
                let inputs = operands(&m, plate, boss);
                let out = boolean(&mut m, k, plate, boss);
                m.rebuild_adjacency();
                record(&format!("mot {name} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A cap that lies in another face's plane** (`cap`): the population where a *disk* is a
    // face of the result and the face around it holds the same circle as a hole. The two are
    // adjacent across that circle and nothing else — a disk has no nodes — so before the merge
    // learned to read circles, this table could not have held the answer: every row here either
    // was refused outright or came back with a plane split in two.
    //
    // ★ **Why the census had nothing to say about that change**: it contained no such fixture at
    // all (measured — the diff over the whole table was empty), which is the same lesson three
    // cells running: **the corpus does not contain the population a cell opens until the cell
    // adds it.** These rows are that addition. The bore row is the negative control: its circle
    // is shared with a *cylinder* face, which is not in the plane's group, so it must stay a hole.
    for (pn, base, h) in [
        ("ontop", [2.0, 2.0, 2.0], 1.0), // a boss standing on the top face (contact)
        ("sunk", [2.0, 2.0, 1.0], 1.0),  // a boss whose cap is flush with the top, inside
        ("through", [2.0, 2.0, 0.0], 3.0), // a boss whose cap is flush with the bottom
        ("bore", [2.0, 2.0, -1.0], 4.0), // the negative control: a real hole
    ] {
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = nacre_ops::fixtures::cuboid(
                &mut m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([4.0, 4.0, 2.0]),
            );
            let b = cylinder(
                &mut m,
                Point3::from_array(base),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                0.5,
                h,
            );
            m.rebuild_adjacency();
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("cap {pn} {kn}"), &m, &inputs, &out);
        }
    }
    // ── **A grid built in two generations** (`xy`): an array joined along x, then joined again
    // along y. What makes this its own population is not the second axis but the second
    // *generation* — a row is a body whose surfaces come from two provenances. Every move here is
    // carried into the statements, so the rows are world-stated; where a row's provenances carry
    // different chains (a recorded motion), moving it puts them on one corner and the corner road
    // answers in the world. The `ct2` family above builds each cell at its place instead.
    //
    // The cells are deliberately feature-poor (one bore) — this table runs twice at every gate,
    // and the user-scale cell costs tens of seconds. What the rows have to pin is the road, and
    // one bore already puts a cylinder gate on a twice-moved wall.
    {
        let cell = |m: &mut Model| {
            let plate = nacre_ops::fixtures::cuboid(
                m,
                Point3::from_array([0.0; 3]),
                Point3::from_array([20.0, 20.0, 10.0]),
            );
            let bore = cylinder(
                m,
                Point3::from_array([6.3, 6.3, -1.0]),
                nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                2.1,
                12.0,
            );
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Cut, plate, bore).expect("the bore cuts")[0];
            m.rebuild_adjacency();
            out
        };
        let row = |m: &mut Model, n: usize, step: [f64; 3]| {
            let mut acc = cell(m);
            for i in 1..n {
                let c = cell(m);
                let c = xf(
                    m,
                    c,
                    Isometry::translation(
                        step.map(|v| Rat::try_from_f64(v * i as f64).expect("rational")),
                    ),
                );
                acc = boolean(m, BoolKind::Fuse, acc, c).expect("the row fuses")[0];
                m.rebuild_adjacency();
            }
            acc
        };
        // x-then-y and y-then-x: the same grid by two routes. The second row is a *moved result*,
        // which is the shape the third door exists for.
        for (an, first, second) in [
            ("grid", [20.0, 0.0, 0.0], [0.0, 20.0, 0.0]),
            ("yx", [0.0, 20.0, 0.0], [20.0, 0.0, 0.0]),
        ] {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = row(&mut m, 2, first);
                let b = row(&mut m, 2, first);
                let b = xf(
                    &mut m,
                    b,
                    Isometry::translation(second.map(|v| Rat::try_from_f64(v).expect("rational"))),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("xy {an} {kn}"), &m, &inputs, &out);
            }
        }
        // Three in x, then y: the moved row carries three chains, not two.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = row(&mut m, 3, [20.0, 0.0, 0.0]);
            let b = row(&mut m, 3, [20.0, 0.0, 0.0]);
            let b = xf(
                &mut m,
                b,
                Isometry::translation(
                    [0.0, 20.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                ),
            );
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("xy three {kn}"), &m, &inputs, &out);
        }
        // ★★ **The row that actually needs the world road.** Everything above builds whether or
        // not the third door exists — measured, by switching the door off and watching these rows
        // come back byte-identical. A grid only reaches a corner the population gate must judge
        // when the cell carries enough features, and the relation is **not monotone** (two pockets
        // with four bores build; one of those pockets with two *other* bores does not). So this row holds
        // the smallest configuration measured to refuse `CylinderGateUndecided` with the door off,
        // taken from the user's own part: one pocket and two through-bores in an 86 mm cell.
        // Treat the numbers as pinned — nudging one moves the row out of its population.
        {
            let featured = |m: &mut Model| {
                let plate = nacre_ops::fixtures::cuboid(
                    m,
                    Point3::from_array([0.0; 3]),
                    Point3::from_array([86.0, 86.0, 71.5]),
                );
                let pocket = nacre_ops::fixtures::cuboid(
                    m,
                    Point3::from_array([1.0, 38.6, 0.0]),
                    Point3::from_array([33.2, 58.6, 68.5]),
                );
                m.rebuild_adjacency();
                let out = boolean(m, BoolKind::Cut, plate, pocket).expect("the pocket cuts")[0];
                m.rebuild_adjacency();
                let mut out = out;
                for (c, r) in [([17.1, 8.9], 2.22), ([46.75, 37.35], 2.34)] {
                    let bore = cylinder(
                        m,
                        Point3::from_array([c[0], c[1], 0.0]),
                        nacre_math::Vector3::from_array([0.0, 0.0, 1.0]),
                        r,
                        143.0,
                    );
                    m.rebuild_adjacency();
                    out = boolean(m, BoolKind::Cut, out, bore).expect("the bore cuts")[0];
                    m.rebuild_adjacency();
                }
                out
            };
            let featured_row = |m: &mut Model| {
                let a = featured(m);
                let b = featured(m);
                let b = xf(
                    m,
                    b,
                    Isometry::translation(
                        [86.0, 0.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                    ),
                );
                let out = boolean(m, BoolKind::Fuse, a, b).expect("the row fuses")[0];
                m.rebuild_adjacency();
                out
            };
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = featured_row(&mut m);
                let b = featured_row(&mut m);
                let b = xf(
                    &mut m,
                    b,
                    Isometry::translation(
                        [0.0, 86.0, 0.0].map(|v: f64| Rat::try_from_f64(v).expect("rational")),
                    ),
                );
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("xy needsworld {kn}"), &m, &inputs, &out);
            }
        }
        // The negative control: a turned row has no world statement, so it must keep taking the
        // old roads — a rotation does not fold into a translation and this door never opens for it.
        for (kn, k) in KINDS {
            let mut m = Model::new();
            let a = row(&mut m, 2, [20.0, 0.0, 0.0]);
            let b = row(&mut m, 2, [20.0, 0.0, 0.0]);
            let b = xf(
                &mut m,
                b,
                Isometry::rotation(Rotation {
                    axis: Axis::Z,
                    pivot: [Rat::from_int(0); 3],
                    angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
                }),
            );
            let inputs = operands(&m, a, b);
            let out = boolean(&mut m, k, a, b);
            m.rebuild_adjacency();
            record(&format!("xy turned {kn}"), &m, &inputs, &out);
        }
    }
    // ── **A datum through a turned solid's corners, and a boolean over it** (`dt`).
    //
    // ★★★★★ **This family exists because its absence hid a regression for 169 commits.** Every
    // datum above is `DatumDef::Stated`; nothing in this corpus stated one **through vertices**,
    // so when the invariant-plane restatement made `vertex_meet` read a turned solid's corners as
    // straddling — costing every one of them the named datum road — the census crossed it
    // bit-identical and said nothing. The rule is to put the
    // population a change touches into the ledger, and this is that population.
    //
    // ★★ **It does not stop at making the datum.** The plane a `ThroughVertices` datum mints is
    // where `collect_planes` has to choose between the exact witness road and the judged one, and
    // that fork is only reached when such a plane meets a boolean. So the datum hosts a prism and
    // the prism is an operand: a row that only *built* a datum would never walk the road worth
    // watching.
    for (kn, k) in KINDS {
        let mut m = Model::new();
        let block = nacre_ops::fixtures::cuboid(
            &mut m,
            Point3::from_array([0.0; 3]),
            Point3::from_array([2.0, 2.0, 2.0]),
        );
        let turned = xf(
            &mut m,
            block,
            Isometry::rotation(nacre_exact::Rotation {
                axis: Axis::Z,
                pivot: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(37)).unwrap(),
            }),
        );
        // Three corners of the turned block, in handle order so the statement is replay-stable.
        let mut vs: Vec<Handle<nacre_topo::Vertex>> = Vec::new();
        {
            let sol = m.solid(turned);
            'pick: for &sh in std::iter::once(&sol.outer).chain(sol.cavities.iter()) {
                for &fh in &m.shell(sh).faces {
                    for &he in &m.face(fh).outer.half_edges {
                        let vh = m.he_start(he);
                        if !vs.contains(&vh) {
                            vs.push(vh);
                        }
                        if vs.len() == 3 {
                            break 'pick;
                        }
                    }
                }
            }
        }
        vs.sort_by_key(|v| v.index());
        let Ok(OpOutput::DatumPlane { frame, .. }) = apply(
            &mut m,
            &Operation::DatumPlane {
                def: DatumDef::ThroughVertices([vs[0], vs[1], vs[2]]),
            },
        ) else {
            panic!("a datum through three corners of a turned block")
        };
        let OpOutput::Extrude { solid: tool, .. } = apply(
            &mut m,
            &Operation::Extrude {
                frame,
                profile: nacre_ops::Profile2d::polygon(vec![
                    nacre_math::Point2::from_array([-0.6, -0.6]),
                    nacre_math::Point2::from_array([0.6, -0.6]),
                    nacre_math::Point2::from_array([0.6, 0.6]),
                    nacre_math::Point2::from_array([-0.6, 0.6]),
                ])
                .unwrap(),
                dist: 1.4,
            },
        )
        .expect("a prism on the datum") else {
            unreachable!("extrude yields Extrude output")
        };
        m.rebuild_adjacency();
        let inputs = operands(&m, turned, tool);
        let out = boolean(&mut m, k, turned, tool);
        m.rebuild_adjacency();
        record(&format!("dt {kn}"), &m, &inputs, &out);
    }

    // ★★★★★ **The link that turns "interning explains it" into something falsifiable.**
    // Every plane with points is named (wide ones in the arbitrary-precision vessel), so a
    // coordinate can move only when wide planes *merge* — a `c ` line that moves must come with
    // a nonzero count here. If the coordinates move and this stays zero, the cause is something
    // else and the diff is not explained.
    println!(
        "stat wide_planes {}",
        nacre_topo::WIDE_PLANES.load(std::sync::atomic::Ordering::Relaxed)
    );
    // The falsifiability bridge for seeding-shaped changes. A `c ` diff that the wide
    // counter cannot explain (name-collision populations are narrow) must come with a nonzero
    // count here instead — pushes that interned onto a seeded world plane.
    println!(
        "stat seeded_hits {}",
        nacre_topo::SEEDED_HITS.load(std::sync::atomic::Ordering::Relaxed)
    );
    // ── **Arc profiles**: sketched circles and arcs extruded, then met by a box. The
    // half disk's chord wall *crosses* its cylinder (pierce corners `Lo`/`Hi`), the annulus and
    // the bored plate have only whole circles; the slot's straight walls are *tangent* to its half
    // cylinders (its own fillet-like joints, no tangency row), and the box's wall `x = 1` touches
    // the left end's circle on the side the slot does not have — nothing touches, and all three
    // build.
    {
        let p2 = |x: f64, y: f64| nacre_math::Point2::from_array([x, y]);
        let sketch = |m: &mut Model, edges: Vec<Stated>, dist: f64| -> Handle<Solid> {
            let profile = stated(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, Axis::Z);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the arc profile extrudes") else {
                unreachable!()
            };
            solid
        };
        type Shape = Box<dyn Fn(&mut Model) -> Handle<Solid>>;
        let shapes: [(&str, Shape); 4] = [
            (
                "halfdisk",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            line(p2(2.0, 5.0), p2(2.0, -1.0)),
                            arc_turns(p2(2.0, 2.0), p2(2.0, -1.0), 2),
                        ],
                        3.0,
                    )
                }),
            ),
            (
                "annulus",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![circle(p2(2.0, 2.0), 3.0), circle(p2(2.0, 2.0), 1.5)],
                        3.0,
                    )
                }),
            ),
            (
                "bored",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            line(p2(-1.0, -1.0), p2(5.0, -1.0)),
                            line(p2(5.0, -1.0), p2(5.0, 5.0)),
                            line(p2(5.0, 5.0), p2(-1.0, 5.0)),
                            line(p2(-1.0, 5.0), p2(-1.0, -1.0)),
                            circle(p2(2.0, 2.0), 1.0),
                        ],
                        3.0,
                    )
                }),
            ),
            (
                "slot",
                Box::new(move |m| {
                    sketch(
                        m,
                        vec![
                            line(p2(0.0, 1.0), p2(4.0, 1.0)),
                            arc_turns(p2(4.0, 2.0), p2(4.0, 1.0), 2),
                            line(p2(4.0, 3.0), p2(0.0, 3.0)),
                            arc_turns(p2(0.0, 2.0), p2(0.0, 3.0), 2),
                        ],
                        3.0,
                    )
                }),
            ),
        ];
        for (sn, build) in &shapes {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let a = build(&mut m);
                let b = nacre_ops::fixtures::cuboid(
                    &mut m,
                    Point3::from_array([1.0, 0.0, 1.0]),
                    Point3::from_array([6.0, 4.0, 5.0]),
                );
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arcprofile {sn} {kn}"), &m, &inputs, &out);
            }
        }
    }

    // ── **Arc walls**: the user's first assembly and its controls. A filleted plate
    // with two holes (XY), a standing plate with a slot window (ZX, moved to y ∈ [3, 4]) and two
    // triangular gussets (YZ) met three walls at once — same-solid parallel cylinder pairs, the
    // oblique gusset plane, and the fillets' tangent rulings. Every row is a pair the assembly's
    // fold visits or a control that isolates one wall.
    {
        // Helpers as items, not closures: the shape builders below are boxed and `move`d, and a
        // closure they borrowed would not live long enough.
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64) -> Handle<Solid> {
            let profile = stated(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the arc-wall profile extrudes") else {
                unreachable!()
            };
            solid
        }
        // An exact translation by decimals — the anchor arithmetic a script hands the kit.
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn far_box(m: &mut Model) -> Handle<Solid> {
            nacre_ops::fixtures::cuboid(
                m,
                Point3::from_array([20.0, 20.0, 20.0]),
                Point3::from_array([21.0, 21.0, 21.0]),
            )
        }
        fn circle_prism(m: &mut Model, c: [f64; 2], r: f64, dist: f64) -> Handle<Solid> {
            prism(m, Axis::Z, vec![circle(p2(c[0], c[1]), r)], dist)
        }
        // p1: the plate, bottom corners filleted r 2 with the holes centred on the fillet axes.
        let p1 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-1.5, -4.0), p2(1.5, -4.0)),
                    arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                    line(p2(3.5, -2.0), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -2.0)),
                    arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1),
                    circle(p2(-1.5, -2.0), 1.0),
                    circle(p2(1.5, -2.0), 1.0),
                ],
                1.0,
            )
        };
        // p2: the standing plate on ZX (sketch u = z, v = x) with a 2×1 slot window, moved to y = 3.
        let p2p = |m: &mut Model| {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    line(p2(1.0, -3.5), p2(6.0, -3.5)),
                    line(p2(6.0, -3.5), p2(6.0, 3.5)),
                    line(p2(6.0, 3.5), p2(1.0, 3.5)),
                    line(p2(1.0, 3.5), p2(1.0, -3.5)),
                    line(p2(3.0, -0.5), p2(4.0, -0.5)),
                    arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2),
                    line(p2(4.0, 0.5), p2(3.0, 0.5)),
                    arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        };
        // p3: a triangular gusset on YZ (sketch u = y, v = z), moved to x = 1.5 as the script does.
        let p3 = |m: &mut Model| {
            let s = prism(
                m,
                Axis::X,
                vec![
                    line(p2(0.0, 1.0), p2(3.0, 1.0)),
                    line(p2(3.0, 1.0), p2(3.0, 6.0)),
                    line(p2(3.0, 6.0), p2(0.0, 1.0)),
                ],
                1.0,
            );
            shift(m, s, [1.5, 0.0, 0.0])
        };
        // Controls.
        let fillet15 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-2.0, -4.0), p2(2.0, -4.0)),
                    arc_turns(p2(2.0, -2.5), p2(2.0, -4.0), 1),
                    line(p2(3.5, -2.5), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -2.5)),
                    arc_turns(p2(-2.0, -2.5), p2(-3.5, -2.5), 1),
                ],
                1.0,
            )
        };
        let slot30 = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(0.0, -5.0), p2(30.0, -5.0)),
                    arc_turns(p2(30.0, 0.0), p2(30.0, -5.0), 2),
                    line(p2(30.0, 5.0), p2(0.0, 5.0)),
                    arc_turns(p2(0.0, 0.0), p2(0.0, 5.0), 2),
                ],
                2.0,
            )
        };
        let holes_plate = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-3.5, -4.0), p2(3.5, -4.0)),
                    line(p2(3.5, -4.0), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -4.0)),
                    circle(p2(-1.5, -2.0), 1.0),
                    circle(p2(1.5, -2.0), 1.0),
                ],
                1.0,
            )
        };
        let rrect = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(2.0, 0.0), p2(8.0, 0.0)),
                    arc_turns(p2(8.0, 2.0), p2(8.0, 0.0), 1),
                    line(p2(10.0, 2.0), p2(10.0, 6.0)),
                    arc_turns(p2(8.0, 6.0), p2(10.0, 6.0), 1),
                    line(p2(8.0, 8.0), p2(2.0, 8.0)),
                    arc_turns(p2(2.0, 6.0), p2(2.0, 8.0), 1),
                    line(p2(0.0, 6.0), p2(0.0, 2.0)),
                    arc_turns(p2(2.0, 2.0), p2(0.0, 2.0), 1),
                ],
                3.0,
            )
        };
        let tube = |m: &mut Model| {
            prism(
                m,
                Axis::Z,
                vec![circle(p2(2.0, 2.0), 3.0), circle(p2(2.0, 2.0), 1.5)],
                3.0,
            )
        };
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            ("p1p2", Box::new(move |m| (p1(m), p2p(m)))),
            ("p1p3", Box::new(move |m| (p1(m), p3(m)))),
            ("p2p3", Box::new(move |m| (p2p(m), p3(m)))),
            ("fillet15-far", Box::new(move |m| (fillet15(m), far_box(m)))),
            ("slot-far", Box::new(move |m| (slot30(m), far_box(m)))),
            ("holes-gusset", Box::new(move |m| (holes_plate(m), p3(m)))),
            (
                "rrect-box",
                Box::new(move |m| {
                    let a = rrect(m);
                    let b = nacre_ops::fixtures::cuboid(
                        m,
                        Point3::from_array([3.0, 2.0, 1.0]),
                        Point3::from_array([7.0, 6.0, 5.0]),
                    );
                    (a, b)
                }),
            ),
            (
                "bushing",
                Box::new(move |m| {
                    let a = tube(m);
                    let pin = circle_prism(m, [2.0, 2.0], 1.0, 5.0);
                    (a, shift(m, pin, [0.0, 0.0, -1.0]))
                }),
            ),
            (
                "stacked",
                Box::new(move |m| {
                    let boss = circle_prism(m, [0.0, 0.0], 2.0, 2.0);
                    let pin = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    (boss, shift(m, pin, [0.0, 0.0, 2.0]))
                }),
            ),
            (
                "stacked-same",
                Box::new(move |m| {
                    let a = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    let b = circle_prism(m, [0.0, 0.0], 1.0, 2.0);
                    (a, shift(m, b, [0.0, 0.0, 2.0]))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("arcwalls {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **Four-plane operand vertices**: a gusset whose apex lands exactly on a wall's top
    // edge makes a result vertex where **four faces** meet. Fed back as an operand, that vertex
    // named once per face would carry four names, one of them a triple whose planes share a
    // line — and the next boolean would refuse by whichever symptom its build order met first
    // (`DegenerateWitness` in the parallel build, `FourPlane` sequentially; measured). The
    // vertex names itself from its incident classes and every pair builds. The rows hold
    // both operand orders, the user's four-part fold in both fold orders, the
    // near misses (apex above and below the edge), a box that shares only the top plane's class,
    // and the mirrored fold (class numbers permuted).
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64) -> Handle<Solid> {
            let profile = stated(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the four-plane profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn mirror_x(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
            let OpOutput::Mirror { solid } = apply(
                m,
                &Operation::Mirror {
                    solid: s,
                    axis: Axis::X,
                    offset: Rat::from_int(0),
                },
            )
            .expect("the mirror applies") else {
                unreachable!()
            };
            solid
        }
        // A fuse inside a fixture: the operand a later row feeds back.
        fn fuse(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Fuse, a, b).expect("the fixture's own fuse builds");
            assert_eq!(out.len(), 1, "the fixture's own fuse is one body");
            out[0]
        }
        // The standing wall: 7 × 5 × 1 at y ∈ [3, 4], z ∈ [1, 6] (the slot plate without its slot).
        fn wall(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    line(p2(1.0, -3.5), p2(6.0, -3.5)),
                    line(p2(6.0, -3.5), p2(6.0, 3.5)),
                    line(p2(6.0, 3.5), p2(1.0, 3.5)),
                    line(p2(1.0, 3.5), p2(1.0, -3.5)),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        // The user's slot plate: the wall with its 2 × 1 slot window.
        fn slot_plate(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    line(p2(1.0, -3.5), p2(6.0, -3.5)),
                    line(p2(6.0, -3.5), p2(6.0, 3.5)),
                    line(p2(6.0, 3.5), p2(1.0, 3.5)),
                    line(p2(1.0, 3.5), p2(1.0, -3.5)),
                    line(p2(3.0, -0.5), p2(4.0, -0.5)),
                    arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2),
                    line(p2(4.0, 0.5), p2(3.0, 0.5)),
                    arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        // The user's filleted, bored plate.
        fn plate(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-1.5, -4.0), p2(1.5, -4.0)),
                    arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                    line(p2(3.5, -2.0), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -2.0)),
                    arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1),
                    circle(p2(-1.5, -2.0), 1.0),
                    circle(p2(1.5, -2.0), 1.0),
                ],
                1.0,
            )
        }
        // The gusset on YZ (u = y, v = z), apex at `(3, apex)`, moved to `x`.
        fn gusset(m: &mut Model, x: f64, apex: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::X,
                vec![
                    line(p2(0.0, 1.0), p2(3.0, 1.0)),
                    line(p2(3.0, 1.0), p2(3.0, apex)),
                    line(p2(3.0, apex), p2(0.0, 1.0)),
                ],
                1.0,
            );
            shift(m, s, [x, 0.0, 0.0])
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            (
                "wall+g",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 6.0))
                }),
            ),
            (
                "wall+g swapped",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let b = fuse(m, w, g1);
                    (gusset(m, -2.4, 6.0), b)
                }),
            ),
            (
                "user4 1234",
                Box::new(|m| {
                    let p1 = plate(m);
                    let p2s = slot_plate(m);
                    let ab = fuse(m, p1, p2s);
                    let g1 = gusset(m, 1.4, 6.0);
                    let abc = fuse(m, ab, g1);
                    (abc, gusset(m, -2.4, 6.0))
                }),
            ),
            (
                "user4 1342",
                Box::new(|m| {
                    let p1 = plate(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let ac = fuse(m, p1, g1);
                    let g2 = gusset(m, -2.4, 6.0);
                    let acd = fuse(m, ac, g2);
                    (acd, slot_plate(m))
                }),
            ),
            (
                "apex-5.9",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 5.9);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 5.9))
                }),
            ),
            (
                "apex-6.1",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.1);
                    let a = fuse(m, w, g1);
                    (a, gusset(m, -2.4, 6.1))
                }),
            ),
            (
                "far-box",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    let b = nacre_ops::fixtures::cuboid(
                        m,
                        Point3::from_array([-3.0, 0.0, 6.0]),
                        Point3::from_array([-1.0, 2.0, 7.0]),
                    );
                    (a, b)
                }),
            ),
            (
                "mirrored",
                Box::new(|m| {
                    let w = wall(m);
                    let g1 = gusset(m, 1.4, 6.0);
                    let a = fuse(m, w, g1);
                    let a = mirror_x(m, a);
                    let g2 = gusset(m, -2.4, 6.0);
                    (a, mirror_x(m, g2))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("fourplane {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A plane through a fillet's axis**: a class plane that holds a fillet's axis
    // has two rulings on the fillet, and one of them can be the fillet's own tangent ruling with
    // the plate's wall. That line then carries names from two vocabularies — the operand's tangent
    // corner (`Pierce … Double`) and the class's ruling crossing (`Pierce … Lo/Hi`) — with nothing
    // that knows they are one point, and the lateral's ruling sweep declines by name. The rows:
    // the minimal plate-and-slab in both contact shapes, an off-axis control, and the user's
    // four-part fold at the script's own gusset positions (`1.5`/`−2.5`).
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64) -> Handle<Solid> {
            let profile = stated(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the tangent-line profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        fn fuse(m: &mut Model, a: Handle<Solid>, b: Handle<Solid>) -> Handle<Solid> {
            m.rebuild_adjacency();
            let out = boolean(m, BoolKind::Fuse, a, b).expect("the fixture's own fuse builds");
            assert_eq!(out.len(), 1, "the fixture's own fuse is one body");
            out[0]
        }
        // A 7 × 8 plate with one fillet (r 2) at its bottom-right corner: axis at (1.5, −2), tangent
        // to the bottom wall `y = −4` along `x = 1.5` and to the right wall `x = 3.5` along `y = −2`.
        fn plate1(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-3.5, -4.0), p2(1.5, -4.0)),
                    arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                    line(p2(3.5, -2.0), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -4.0)),
                ],
                1.0,
            )
        }
        fn slab(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
            nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi))
        }
        // The user's parts: the filleted, bored plate, the slot plate, the gusset at `x`.
        fn plate(m: &mut Model) -> Handle<Solid> {
            prism(
                m,
                Axis::Z,
                vec![
                    line(p2(-1.5, -4.0), p2(1.5, -4.0)),
                    arc_turns(p2(1.5, -2.0), p2(1.5, -4.0), 1),
                    line(p2(3.5, -2.0), p2(3.5, 4.0)),
                    line(p2(3.5, 4.0), p2(-3.5, 4.0)),
                    line(p2(-3.5, 4.0), p2(-3.5, -2.0)),
                    arc_turns(p2(-1.5, -2.0), p2(-3.5, -2.0), 1),
                    circle(p2(-1.5, -2.0), 1.0),
                    circle(p2(1.5, -2.0), 1.0),
                ],
                1.0,
            )
        }
        fn slot_plate(m: &mut Model) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    line(p2(1.0, -3.5), p2(6.0, -3.5)),
                    line(p2(6.0, -3.5), p2(6.0, 3.5)),
                    line(p2(6.0, 3.5), p2(1.0, 3.5)),
                    line(p2(1.0, 3.5), p2(1.0, -3.5)),
                    line(p2(3.0, -0.5), p2(4.0, -0.5)),
                    arc_turns(p2(4.0, 0.0), p2(4.0, -0.5), 2),
                    line(p2(4.0, 0.5), p2(3.0, 0.5)),
                    arc_turns(p2(3.0, 0.0), p2(3.0, 0.5), 2),
                ],
                1.0,
            );
            shift(m, s, [0.0, 3.0, 0.0])
        }
        fn gusset(m: &mut Model, x: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::X,
                vec![
                    line(p2(0.0, 1.0), p2(3.0, 1.0)),
                    line(p2(3.0, 1.0), p2(3.0, 6.0)),
                    line(p2(3.0, 6.0), p2(0.0, 1.0)),
                ],
                1.0,
            );
            shift(m, s, [x, 0.0, 0.0])
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            // The slab stands on the plate; its face `x = 1.5` holds the fillet's axis. It stops
            // at `y = −3`, short of the plate's bottom wall, so that plane is the only coincidence.
            (
                "fillet-slab",
                Box::new(|m| (plate1(m), slab(m, [1.5, -3.0, 1.0], [2.5, 0.0, 2.0]))),
            ),
            // Same plane, the slab far from the plate: the class still cuts the fillet.
            (
                "fillet-slab far",
                Box::new(|m| (plate1(m), slab(m, [1.5, 10.0, 1.0], [2.5, 12.0, 2.0]))),
            ),
            // Off the axis by 0.1: an ordinary offset wall — the control.
            (
                "slab 1.6",
                Box::new(|m| (plate1(m), slab(m, [1.6, -3.0, 1.0], [2.6, 0.0, 2.0]))),
            ),
            // ★ A slab reaching the plate's bottom wall plane `y = −4`
            // — its face coplanar with the wall, not overlapping it, its bottom edge on the
            // wall's line past the tangent point. On the axis plane the slab's corner edge lies
            // on the fillet's own seam, which the rulings split does not arrange yet
            // (`RulingBoundNotYet`, frozen here by name). Off it the slab rests on the plate's
            // top across the fillet's rim, and the cut leaves the plate as it was — the rim
            // the cleaned cap holds, not the one the arrangement split (`draft::HeldRims`); the
            // wall plane plays no part (the coverage lock moves the slab off it).
            (
                "slab wall 1.5",
                Box::new(|m| (plate1(m), slab(m, [1.5, -4.0, 1.0], [2.5, 0.0, 2.0]))),
            ),
            (
                "slab wall 1.6",
                Box::new(|m| (plate1(m), slab(m, [1.6, -4.0, 1.0], [2.6, 0.0, 2.0]))),
            ),
            // The user's fold at the script's own gusset position: the third fuse of the script
            // (plate + slot plate, then the gusset at `1.5`). The fourth (`−2.5`) is the ops lock's.
            (
                "user 1.5",
                Box::new(|m| {
                    let p1 = plate(m);
                    let p2s = slot_plate(m);
                    let ab = fuse(m, p1, p2s);
                    (ab, gusset(m, 1.5))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("tangentline {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }
    // ── **A fully rounded outline**: fillet *every* corner of a plate and its cap ring
    // has no three-plane corner left — eight tangencies and nothing else. The predicate that asks
    // "is this ring inside that bore's disk" looks for a witness among three-plane names only, so
    // it finds none and refuses, and the plate cannot enter **any** boolean. The rows say exactly
    // that: three partners that share nothing but the plate (one standing on it, one overlapping,
    // one a hundred units away), and two controls that differ in one thing each — one corner left
    // sharp, and no bores. The last row carries the user's own script.
    {
        fn p2(x: f64, y: f64) -> nacre_math::Point2 {
            nacre_math::Point2::from_array([x, y])
        }
        fn prism(m: &mut Model, axis: Axis, edges: Vec<Stated>, dist: f64) -> Handle<Solid> {
            let profile = stated(edges).expect("a valid profile").remove(0);
            let frame = SketchFrame::world(m, axis);
            let OpOutput::Extrude { solid, .. } = apply(
                m,
                &Operation::Extrude {
                    frame,
                    profile,
                    dist,
                },
            )
            .expect("the round-plate profile extrudes") else {
                unreachable!()
            };
            solid
        }
        fn shift(m: &mut Model, s: Handle<Solid>, t: [f64; 3]) -> Handle<Solid> {
            let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
            xf(m, s, Isometry::translation([r(t[0]), r(t[1]), r(t[2])]))
        }
        // A 90 × 50 × 12 plate whose corners are filleted (r 5) and which carries `bores` holes
        // (d 7) at `(±38, ±18)`. `fillets` says how many of the four corners are rounded, walking
        // counter-clockwise from the bottom-right: with four the outline has **no** sharp corner.
        fn plate(m: &mut Model, fillets: usize, bores: usize) -> Handle<Solid> {
            // The corners in walk order, each with the direction in and the direction out.
            let corner = [
                ([45.0, -25.0], [1.0, 0.0], [0.0, 1.0]),
                ([45.0, 25.0], [0.0, 1.0], [-1.0, 0.0]),
                ([-45.0, 25.0], [-1.0, 0.0], [0.0, -1.0]),
                ([-45.0, -25.0], [0.0, -1.0], [1.0, 0.0]),
            ];
            let r = 5.0;
            let mut edges: Vec<Stated> = Vec::new();
            // Where the previous corner left the pen.
            let mut at = {
                let (c, _, d_out) = corner[3];
                if fillets > 3 {
                    [c[0] + r * d_out[0], c[1] + r * d_out[1]]
                } else {
                    c
                }
            };
            for (i, (c, d_in, d_out)) in corner.into_iter().enumerate() {
                if i < fillets {
                    let tin = [c[0] - r * d_in[0], c[1] - r * d_in[1]];
                    let centre = [tin[0] + r * d_out[0], tin[1] + r * d_out[1]];
                    edges.push(line(p2(at[0], at[1]), p2(tin[0], tin[1])));
                    edges.push(arc_turns(p2(centre[0], centre[1]), p2(tin[0], tin[1]), 1));
                    at = [c[0] + r * d_out[0], c[1] + r * d_out[1]];
                } else {
                    edges.push(line(p2(at[0], at[1]), p2(c[0], c[1])));
                    at = c;
                }
            }
            for &[x, y] in [[38.0, 18.0], [38.0, -18.0], [-38.0, 18.0], [-38.0, -18.0]]
                .iter()
                .take(bores)
            {
                edges.push(circle(p2(x, y), 3.5));
            }
            prism(m, Axis::Z, edges, 12.0)
        }
        // The user's rib: a stepped profile standing on the plate's top face, extruded ±20 in y.
        fn rib(m: &mut Model, x: f64) -> Handle<Solid> {
            let s = prism(
                m,
                Axis::Y,
                vec![
                    line(p2(12.0, 0.0), p2(12.0, 7.5)),
                    line(p2(12.0, 7.5), p2(27.0, 7.5)),
                    line(p2(27.0, 7.5), p2(27.0, 5.5)),
                    line(p2(27.0, 5.5), p2(62.0, 5.5)),
                    line(p2(62.0, 5.5), p2(62.0, -5.5)),
                    line(p2(62.0, -5.5), p2(27.0, -5.5)),
                    line(p2(27.0, -5.5), p2(27.0, -7.5)),
                    line(p2(27.0, -7.5), p2(12.0, -7.5)),
                    line(p2(12.0, -7.5), p2(12.0, 0.0)),
                ],
                40.0,
            );
            shift(m, s, [x, -20.0, 0.0])
        }
        fn box_at(m: &mut Model, lo: [f64; 3], hi: [f64; 3]) -> Handle<Solid> {
            nacre_ops::fixtures::cuboid(m, Point3::from_array(lo), Point3::from_array(hi))
        }
        type Pair = Box<dyn Fn(&mut Model) -> (Handle<Solid>, Handle<Solid>)>;
        let pairs: Vec<(&str, Pair)> = vec![
            // Standing on the plate's top face.
            (
                "box on top",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // Overlapping it.
            (
                "box through",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [10.0, -20.0, 6.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // ★ A hundred units away — the partner has nothing to do with it. This row is what
            // says the wall is the plate's own.
            (
                "box far",
                Box::new(|m| {
                    (
                        plate(m, 4, 4),
                        box_at(m, [200.0, -5.0, -5.0], [210.0, 5.0, 5.0]),
                    )
                }),
            ),
            // One corner left sharp — the control that differs in exactly one thing.
            (
                "one sharp corner",
                Box::new(|m| {
                    (
                        plate(m, 3, 4),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // Rounded but unbored — the other control.
            (
                "no bores",
                Box::new(|m| {
                    (
                        plate(m, 4, 0),
                        box_at(m, [10.0, -20.0, 12.0], [25.0, 20.0, 62.0]),
                    )
                }),
            ),
            // The rib the user's script builds, on one side.
            ("rib", Box::new(|m| (plate(m, 4, 4), rib(m, 17.5)))),
            // ★ The user's own script — plate, both ribs, then the cylinder that bores across
            // them. It could not be a row until the fuse built; now it is one, and the name it
            // lands on is the **next** wall, a class carrying a circle and rulings at once. That
            // is the hand-off to the cell after this one, held in the corpus so it cannot be lost.
            (
                "user script",
                Box::new(|m| {
                    let p1 = plate(m, 4, 4);
                    let r1 = rib(m, 17.5);
                    m.rebuild_adjacency();
                    let ab = boolean(m, BoolKind::Fuse, p1, r1).expect("plate + rib")[0];
                    let r2 = rib(m, -17.5);
                    m.rebuild_adjacency();
                    let abc = boolean(m, BoolKind::Fuse, ab, r2).expect("+ the second rib")[0];
                    let c = prism(m, Axis::Z, vec![circle(p2(0.0, 0.0), 10.0)], 90.0);
                    let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
                    let turned = xf(
                        m,
                        c,
                        Isometry::rotation(nacre_exact::Rotation {
                            axis: Axis::Y,
                            pivot: [r(0.0), r(0.0), r(0.0)],
                            angle: nacre_exact::Angle::from_deg(Rat::from_int(90))
                                .expect("a right angle"),
                        }),
                    );
                    // ★ The ops prism spans `0..90` along its axis, so centring the turned
                    // cylinder on the plate is a shift of **−45** — the script's own
                    // `cylinder({h: 90, center: [0,0,0]})`. With `+45` it would sit beside the
                    // plate touching only its side plane, which is a different model.
                    (abc, shift(m, turned, [-45.0, 0.0, 47.0]))
                }),
            ),
            // ★ **A disk that clears *neither* axis of the footprint.** The tool's cap
            // sits on the plate's own side plane, and there it crosses a corner fillet's tangent
            // line (5 away, radius 10) while overlapping the plate's z-span. So the wall rule's
            // answer for it is «does not clear» — the other half of the pair whose first half is
            // the user's script, where the same cap clears by both axes. Without both, the disk
            // arm could be a constant and still be green.
            (
                "disk across the tangent",
                Box::new(|m| {
                    let c = prism(m, Axis::Z, vec![circle(p2(0.0, 0.0), 10.0)], 90.0);
                    let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
                    let turned = xf(
                        m,
                        c,
                        Isometry::rotation(nacre_exact::Rotation {
                            axis: Axis::Y,
                            pivot: [r(0.0), r(0.0), r(0.0)],
                            angle: nacre_exact::Angle::from_deg(Rat::from_int(90))
                                .expect("a right angle"),
                        }),
                    );
                    (plate(m, 4, 0), shift(m, turned, [-45.0, -15.0, 6.0]))
                }),
            ),
            // ★ **Why the class-edge net's population is empty, held as a row.** Two
            // cylinders on crossing axes, one per operand, whose surfaces meet. A class could
            // otherwise carry one's circle and the other's rulings *touching*, which no split
            // cuts; the pair gate refuses the pair (`CylinderPairContact`) before any arrangement
            // runs, and this row is what would notice that stopping.
            (
                "crossing cylinders",
                Box::new(|m| {
                    let circle = |r: f64| vec![circle(p2(0.0, 0.0), r)];
                    let up = prism(m, Axis::Z, circle(5.0), 20.0);
                    let across = prism(m, Axis::Z, circle(5.0), 20.0);
                    let r = |x: f64| Rat::from_decimal(x).expect("a short decimal");
                    let turned = xf(
                        m,
                        across,
                        Isometry::rotation(nacre_exact::Rotation {
                            axis: Axis::Y,
                            pivot: [r(0.0), r(0.0), r(0.0)],
                            angle: nacre_exact::Angle::from_deg(Rat::from_int(90))
                                .expect("a right angle"),
                        }),
                    );
                    (up, shift(m, turned, [-10.0, 0.0, 10.0]))
                }),
            ),
        ];
        for (pn, build) in &pairs {
            for (kn, k) in KINDS {
                let mut m = Model::new();
                let (a, b) = build(&mut m);
                m.rebuild_adjacency();
                let inputs = operands(&m, a, b);
                let out = boolean(&mut m, k, a, b);
                m.rebuild_adjacency();
                record(&format!("roundplate {pn} {kn}"), &m, &inputs, &out);
            }
        }
    }

    // ★★★★ **The realization road for surfaces, measured.** These
    // counters cover *every* surface this process pushed (they are process-global and the
    // dump is this binary's only test): `surface_cache_differs` is exactly how many `in:` plane
    // digests may move, and a
    // `c ` row that moves without it is not explained.
    let d = nacre_topo::surface_derive_counts();
    println!("stat surface_derived {}", d.derived);
    println!("stat surface_declined {}", d.declined);
    println!("stat surface_cache_differs {}", d.differs);
    println!("stat surface_declined_unnamed {}", d.declined_unnamed);
    println!("stat surface_declined_wide {}", d.declined_wide);
    println!("stat surface_declined_motion {}", d.declined_motion);
    println!("stat surface_declined_arith {}", d.declined_arith);
    println!("stat surface_declined_cylinder {}", d.declined_cylinder);

    // The plane-sense lock's reach: how many planes it judged, how many it could not carry to the
    // world, and how many of the judged ride a reflection or a turn.
    let [agree, unmeasured, mirrored, turned] = SENSE
        .each_ref()
        .map(|a| a.load(std::sync::atomic::Ordering::Relaxed));
    println!("stat sense_agree {agree}");
    println!(
        "stat zero_by_coincidence {}",
        ZERO_BY_COINCIDENCE.load(std::sync::atomic::Ordering::Relaxed)
    );
    println!("stat sense_unmeasured {unmeasured}");
    println!("stat sense_mirrored {mirrored}");
    println!("stat sense_turned {turned}");
    let [stated, unstated] = CIRCLES
        .each_ref()
        .map(|a| a.load(std::sync::atomic::Ordering::Relaxed));
    println!("stat rims_stated {stated}");
    println!("stat rims_unstated {unstated}");
}

/// A square prism on a plane through the origin with normal `n` — the tilted twin of [`ex`].
///
/// ★ The footprint is a **U** (a notched square), so two of its walls lie on one plane. Where the
/// kernel can name that plane the two share a `Surface` handle; where it cannot they are two
/// handles for one plane, which is the state this line of work removes.
fn tilted_prism(m: &mut Model, n: [f64; 3], off: f64, size: f64, dist: f64) -> Handle<Solid> {
    use nacre_math::{Point2, Vector3};
    let p = |x: f64, y: f64| Point2::from_array([off + x * size, off + y * size]);
    let plane =
        nacre_ops::SketchPlane::from_origin_normal(Point3::origin(), Vector3::from_array(n))
            .expect("a tilted plane");
    let __g208 = datum_frame(m, plane);
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g208,
            // U-shaped: the notch's two side walls are the coplanar pair.
            profile: nacre_ops::Profile2d::polygon(vec![
                p(0.0, 0.0),
                p(3.0, 0.0),
                p(3.0, 2.0),
                p(2.0, 2.0),
                p(2.0, 1.0),
                p(1.0, 1.0),
                p(1.0, 2.0),
                p(0.0, 2.0),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("tilted extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    solid
}

/// A rectangular prism raised from `z` by `dist` — the construction path production uses.
fn ex(m: &mut Model, z: f64, lo: [f64; 2], hi: [f64; 2], dist: f64) -> Handle<Solid> {
    use nacre_math::Point2;
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let __g207 = datum_frame(
        m,
        nacre_ops::SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, z])),
    );
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame: __g207,
            profile: nacre_ops::Profile2d::polygon(vec![
                p(lo[0], lo[1]),
                p(hi[0], lo[1]),
                p(hi[0], hi[1]),
                p(lo[0], hi[1]),
            ])
            .unwrap(),
            dist,
        },
    )
    .expect("extrude") else {
        unreachable!("extrude yields Extrude output")
    };
    m.rebuild_adjacency();
    solid
}
