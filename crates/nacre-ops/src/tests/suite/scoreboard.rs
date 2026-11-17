//! **The scoreboard's reading of a boolean** — each family run as the four booleans, and the answer
//! checked as the bodies' volumes and, body by body, their face counts, or as the refusal named.

use super::*;

/// The four booleans each family is run as.
#[derive(Clone, Copy, Debug)]
pub(super) enum Op {
    Fuse,
    /// The first operand minus the second.
    CutAb,
    /// The second operand minus the first.
    CutBa,
    Common,
}

pub(super) const OPS: [Op; 4] = [Op::Fuse, Op::CutAb, Op::CutBa, Op::Common];

pub(super) fn run(m: &mut Model, op: Op, a: Handle<Solid>, b: Handle<Solid>) -> BoolResult {
    match op {
        Op::Fuse => crate::boolean(m, BoolKind::Fuse, a, b),
        Op::CutAb => crate::boolean(m, BoolKind::Cut, a, b),
        Op::CutBa => crate::boolean(m, BoolKind::Cut, b, a),
        Op::Common => crate::boolean(m, BoolKind::Common, a, b),
    }
}

pub(super) type BoolResult = Result<Vec<Handle<Solid>>, BoolError>;

/// What a family's boolean comes out as: the volume of each output solid, ascending (none for an
/// empty result), and the solid's face count — or the refusal it names.
pub(super) struct Want {
    volumes: Vec<f64>,
    /// Each body's face count, the bodies in the order of `volumes`.
    faces: Vec<usize>,
    refused: Option<RejectReason>,
}

/// The bodies' volumes and, body by body in the same order, their face counts.
pub(super) fn counted(volumes: &[f64], faces: &[usize]) -> Want {
    assert_eq!(volumes.len(), faces.len(), "one face count per body");
    Want {
        volumes: volumes.to_vec(),
        faces: faces.to_vec(),
        refused: None,
    }
}

pub(super) fn refused(reason: RejectReason) -> Want {
    Want {
        volumes: Vec::new(),
        faces: Vec::new(),
        refused: Some(reason),
    }
}

fn face_count(m: &Model, s: Handle<Solid>) -> usize {
    m.shell(m.solid(s).outer).faces.len()
}

/// **`got` is `want`** — the refusal it names, or a valid model whose output solids have the
/// volumes and, body by body, the face counts stated, and whose mesh covers their faces.
pub(super) fn check(label: &str, m: &Model, got: BoolResult, want: &Want) {
    if let Some(reason) = want.refused {
        assert!(
            matches!(&got, Err(BoolError::Rejected { reason: r, .. }) if *r == reason),
            "{label}: {got:?}"
        );
        return;
    }
    let got = got.unwrap_or_else(|e| panic!("{label}: {e:?}"));
    let bad = nacre_validate::validate(m);
    assert!(bad.is_empty(), "{label}: {bad:?}");
    let mut bodies: Vec<(f64, usize)> = got
        .iter()
        .map(|&s| {
            let v = nacre_props::mass_props(m, s).expect("props").volume;
            (v, face_count(m, s))
        })
        .collect();
    bodies.sort_by(|x, y| x.0.total_cmp(&y.0));
    let volumes: Vec<f64> = bodies.iter().map(|b| b.0).collect();
    assert_eq!(volumes.len(), want.volumes.len(), "{label}: {volumes:?}");
    for (v, w) in volumes.iter().zip(&want.volumes) {
        assert!((v - w).abs() < 1e-9, "{label}: {volumes:?}");
    }
    let faces: Vec<usize> = bodies.iter().map(|b| b.1).collect();
    assert_eq!(faces, want.faces, "{label}");
    mesh_covers_faces(label, m, &got);
}
