//! Fusing a part with rotated copies of itself — and the naming rule that made 30° fail.
//!
//! At 30° the arithmetic conspires: `sin 30° = ½` exactly, so a rotated copy's profile corner
//! lands *precisely* on one of the original's face planes, and four planes concur at that point.
//! The arrangement named it by its canonical (lexicographically first) triple, which included a
//! plane the union keeps no face on — the copy's wall, buried inside the material. The fuse then
//! succeeded while returning a solid that could not describe itself, and the failure surfaced two
//! operations later, when the next `rotateZ` refused to remap a definition pointing at a surface
//! this solid has no face on.
//!
//! So the propositions here are about **what a result says about itself**, not about volume.

use nacre_math::Point2;
use nacre_ops::{
    BoolError, BoolKind, Operation, Profile2d, RejectReason, SketchFrame, apply, boolean,
};
use nacre_scalar::{Angle, Axis, Isometry, Rat};
use nacre_store::Handle;
use nacre_topo::{Model, Solid, VertexDef};

/// The playground script's part: an L-ish plate with a notch, plus a bar across it.
fn unit(m: &mut Model) -> Handle<Solid> {
    let plate = [
        [0.0, 0.0],
        [50.0, 0.0],
        [50.0, 25.0],
        [38.0, 25.0],
        [38.0, 50.0],
        [50.0, 50.0],
        [50.0, 75.0],
        [0.0, 75.0],
    ];
    let p1 = prism(m, &plate, Axis::Z, 12.0);
    let bar = [[20.0, 12.0], [75.0, 12.0], [75.0, 37.0], [55.0, 37.0]];
    let p2 = prism(m, &bar, Axis::X, 25.0);
    let out = boolean(m, BoolKind::Fuse, p1, p2).expect("the part fuses")[0];
    m.rebuild_adjacency();
    out
}

fn prism(m: &mut Model, pts: &[[f64; 2]], axis: Axis, dist: f64) -> Handle<Solid> {
    let profile =
        Profile2d::polygon(pts.iter().map(|&p| Point2::from_array(p)).collect()).expect("profile");
    let op = Operation::Extrude {
        frame: SketchFrame::world(m, axis),
        profile,
        dist,
    };
    let nacre_ops::OpOutput::Extrude { solid, .. } = apply(m, &op).expect("extrude") else {
        panic!("extrude output")
    };
    m.rebuild_adjacency();
    solid
}

fn copy_of(m: &mut Model, s: Handle<Solid>) -> Handle<Solid> {
    let nacre_ops::OpOutput::Copy { solid } =
        apply(m, &Operation::Copy { solid: s }).expect("copy")
    else {
        panic!("copy output")
    };
    m.rebuild_adjacency();
    solid
}

fn rot_z(m: &mut Model, s: Handle<Solid>, deg: i128) -> Result<Handle<Solid>, nacre_ops::OpError> {
    let out = apply(
        m,
        &Operation::Transform {
            solid: s,
            isometry: Isometry::rotation(nacre_scalar::Rotation {
                axis: Axis::Z,
                point: [Rat::from_int(0); 3],
                angle: Angle::from_deg(Rat::from_int(deg)).unwrap(),
            }),
        },
    )?;
    m.rebuild_adjacency();
    match out {
        nacre_ops::OpOutput::Transform { solid } => Ok(solid),
        o => panic!("{o:?}"),
    }
}

/// Every vertex definition of `solid` that names a surface this solid has no face on — the
/// property `transform` needs and the one the old naming broke. Walked through the public model
/// so the test sees what any consumer would.
fn foreign_definitions(m: &Model, solid: Handle<Solid>) -> Vec<String> {
    use std::collections::HashSet;
    let s = m.solids.get(solid);
    let shells: Vec<_> = std::iter::once(s.outer)
        .chain(s.cavities.iter().copied())
        .collect();
    let mut mine: HashSet<_> = HashSet::new();
    for &sh in &shells {
        for &fh in &m.shells.get(sh).faces {
            mine.insert(m.faces.get(fh).surface);
        }
    }
    let mut bad = Vec::new();
    let mut seen: HashSet<u32> = HashSet::new();
    for &sh in &shells {
        for &fh in &m.shells.get(sh).faces {
            let face = m.faces.get(fh);
            for lp in std::iter::once(&face.outer).chain(face.inner.iter()) {
                for he in &lp.half_edges {
                    for vh in m.edges.get(he.edge).vertices.iter() {
                        if !seen.insert(vh.index()) {
                            continue;
                        }
                        let names = match &m.vertices.get(*vh).def {
                            VertexDef::ThreePlane(p) => p.to_vec(),
                            VertexDef::OnSeam(p) => p.to_vec(),
                        };
                        let missing: Vec<u32> = names
                            .iter()
                            .filter(|s| !mine.contains(s))
                            .map(|s| s.index())
                            .collect();
                        if !missing.is_empty() {
                            let at = m.vertex_point(*vh).as_array();
                            bad.push(format!(
                                "vertex {} at {at:?} names surfaces {missing:?}",
                                vh.index()
                            ));
                        }
                    }
                }
            }
        }
    }
    bad
}

/// ① The fuse that started it: it always built, and what it built was unusable. Both halves are
/// asserted, because the second is the one that was silently false.
#[test]
fn a_fused_result_names_itself_in_its_own_surfaces() {
    let mut m = Model::new();
    let u = unit(&mut m);
    let c = copy_of(&mut m, u);
    let c = rot_z(&mut m, c, 30).expect("rotating a copy");
    let fused = boolean(&mut m, BoolKind::Fuse, u, c).expect("30° fuses")[0];
    m.rebuild_adjacency();
    assert_eq!(
        foreign_definitions(&m, fused),
        Vec::<String>::new(),
        "a result must be describable in the surfaces it actually has"
    );
    // …and the consequence, which is where the failure used to appear: the next rotation of that
    // result is exactly the operation that returned `OriginNotOnSolid`.
    rot_z(&mut m, fused, 30).expect("the result can be moved");
}

/// ② The whole 30° sweep — eleven chained fuses, each on the result of the last.
#[test]
fn the_thirty_degree_sweep_runs_to_completion() {
    let mut m = Model::new();
    let u = unit(&mut m);
    // The template stays live: a boolean supersedes its operands, so `u` itself must never be
    // one of them or the next copy has nothing to copy.
    let mut part = copy_of(&mut m, u);
    for deg in (30..360).step_by(30) {
        let c = copy_of(&mut m, u);
        let c = rot_z(&mut m, c, deg as i128).expect("rotating a copy");
        let out = boolean(&mut m, BoolKind::Fuse, part, c)
            .unwrap_or_else(|e| panic!("fusing the {deg}° copy: {e:?}"));
        assert_eq!(out.len(), 1, "the star stays one body at {deg}°");
        part = out[0];
        m.rebuild_adjacency();
        assert_eq!(
            foreign_definitions(&m, part),
            Vec::<String>::new(),
            "after the {deg}° fuse"
        );
    }
}

/// ③ ★ The negative control, and the baseline for the work that comes after this fix. Three other
/// step angles reject for three *other* reasons; renaming vertices must not touch any of them.
/// If one of these moves, the fix reached further than its argument says it does.
#[test]
fn the_other_rejections_are_untouched() {
    let cases: [(usize, RejectReason); 3] = [
        // Every candidate ray for a containment test has a ring node on it.
        (40, RejectReason::NoClearRay),
        (45, RejectReason::CoplanarMerge),
        // Two bodies meeting along one line — no 2-manifold contains it.
        (120, RejectReason::NonManifoldResultEdge),
    ];
    for (step, expected) in cases {
        let mut m = Model::new();
        let u = unit(&mut m);
        let mut part = copy_of(&mut m, u);
        let mut got = None;
        for deg in (step..360).step_by(step) {
            let c = copy_of(&mut m, u);
            let Ok(c) = rot_z(&mut m, c, deg as i128) else {
                panic!("{step}°: a copy stopped being movable — that is this fix's own subject")
            };
            match boolean(&mut m, BoolKind::Fuse, part, c) {
                Ok(out) => {
                    part = out[0];
                    m.rebuild_adjacency();
                }
                Err(e) => {
                    got = Some(e);
                    break;
                }
            }
        }
        assert_eq!(
            got,
            Some(BoolError::Unsupported { reason: expected }),
            "the {step}° sweep's rejection changed"
        );
    }
}
