//! Unit tests of [`LineStar`], the one answer about the material round a line.

use super::*;
use crate::planes::SolidSide;

/// Spokes `0..n` with these faces and kept sectors.
fn star(faces: &[usize], kept: &[bool]) -> LineStar {
    LineStar::new(faces.to_vec(), kept.to_vec())
}

/// **A face runs on across the line exactly when every sector between its two halves, on the side
/// it keeps, is kept.** A plane through the line (face `0`, spokes `0` and `2`) with a second face
/// ending on it from above (face `1`, spoke `1`): above, the two sectors either side of face `1`;
/// below, one.
#[test]
fn a_face_runs_on_where_its_kept_side_is_kept_between_its_halves() {
    let faces = [0, 1, 0];
    // Kept above on both sides of face 1, not below: the plane's material passes round above.
    assert_eq!(star(&faces, &[true, true, false]).continues(0), Some(true));
    // One sector above not kept: the plane ends on the line.
    assert_eq!(
        star(&faces, &[true, false, false]).continues(0),
        Some(false)
    );
    // Kept below only: the side it keeps is one sector, kept.
    assert_eq!(star(&faces, &[false, false, true]).continues(0), Some(true));
    // Kept all round, or nowhere: the plane is not in the result.
    assert_eq!(star(&faces, &[true, true, true]).continues(0), None);
    assert_eq!(star(&faces, &[false, false, false]).continues(0), None);
}

/// **A face with one half on the line ends there**, wherever it is in the result; one with both
/// halves but only one of them in the result ends too.
#[test]
fn a_face_with_one_half_in_the_result_ends() {
    let faces = [0, 1, 0];
    assert_eq!(
        star(&faces, &[true, false, false]).continues(1),
        Some(false)
    );
    assert_eq!(star(&faces, &[true, true, false]).continues(1), None);
    // Face 0: its spoke 0 sits between `false` and `true` (in), spoke 2 between `true` and
    // `true` (out).
    let faces = [0, 1, 0, 2];
    assert_eq!(
        star(&faces, &[true, true, true, false]).continues(0),
        Some(false)
    );
}

/// **The kept sectors fall into runs round the line** — none kept, all kept, and two runs apart.
#[test]
fn the_kept_sectors_fall_into_runs() {
    let faces = [0, 1, 1, 0];
    assert_eq!(star(&faces, &[false; 4]).kept_runs(), 0);
    assert_eq!(star(&faces, &[true; 4]).kept_runs(), 1);
    assert_eq!(star(&faces, &[true, false, true, false]).kept_runs(), 2);
    assert_eq!(star(&faces, &[true, true, false, false]).kept_runs(), 1);
    assert_eq!(star(&faces, &[false, true, true, true]).kept_runs(), 1);
}

/// **A sector lies in a solid as the nearest of its spokes clockwise says** — and a solid that
/// does not bring two spokes with their material facing each other is no valid solid about a
/// line.
#[test]
fn the_sectors_are_walked_from_each_solids_two_halves() {
    let s = |face, solid, material_ccw| Spoke {
        face,
        solid,
        material_ccw,
    };
    // A's face through the line (material above: ccw from spoke 0), B's two faces ending on it
    // above with the wedge between them B's.
    let spokes = [
        s(0, SolidSide::A, true),
        s(1, SolidSide::B, true),
        s(2, SolidSide::B, false),
        s(0, SolidSide::A, false),
    ];
    assert_eq!(
        LineStar::sectors_of(&spokes),
        Some(vec![
            (true, false),
            (true, true),
            (true, false),
            (false, false)
        ])
    );
    let cut = LineStar::from_spokes(BoolKind::Cut, &spokes).expect("a valid pair");
    // `A − B` keeps A above outside B's wedge: A's face ends on the line.
    assert_eq!(cut.continues(0), Some(false));
    assert_eq!(cut.kept_runs(), 2);
    let fuse = LineStar::from_spokes(BoolKind::Fuse, &spokes).expect("a valid pair");
    assert_eq!(fuse.continues(0), Some(true));
    // Three spokes of one solid, or two with one material side, are refused.
    let mut three = spokes.to_vec();
    three.push(s(3, SolidSide::B, true));
    assert_eq!(LineStar::sectors_of(&three), None);
    let mut same = spokes;
    same[2].material_ccw = true;
    assert_eq!(LineStar::sectors_of(&same), None);
}
