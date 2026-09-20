use super::*;

fn uv(v: &[[f64; 2]]) -> Vec<P2> {
    v.to_vec()
}

/// Every piece is `v`-monotone: walking it, the sweep order rises to exactly one
/// peak and falls to exactly one trough. That is the whole proposition of this
/// module — a piece that fails it cannot be triangulated by the stack algorithm,
/// and a decomposition that produces one has not decomposed anything.
fn assert_monotone(uv: &[P2], piece: &[usize]) {
    let n = piece.len();
    assert!(n >= 3, "piece {piece:?} is not a polygon");
    let turns = (0..n)
        .filter(|&i| {
            let (p, c, q) = (
                uv[piece[(i + n - 1) % n]],
                uv[piece[i]],
                uv[piece[(i + 1) % n]],
            );
            lex_less(c, p) == lex_less(c, q)
        })
        .count();
    assert_eq!(turns, 2, "piece {piece:?} has {turns} extrema, want 2");
}

fn check(pts: &[[f64; 2]], rings: &[&[usize]]) -> Vec<Vec<usize>> {
    let uv = uv(pts);
    let pieces = decompose(&uv, rings, None).expect("decomposes");
    for p in &pieces {
        assert_monotone(&uv, p);
    }
    // The pieces partition the polygon, so their areas sum to its area.
    let area = |r: &[usize]| -> f64 {
        let mut s = 0.0;
        for k in 0..r.len() {
            let (a, b) = (uv[r[k]], uv[r[(k + 1) % r.len()]]);
            s += a[0] * b[1] - b[0] * a[1];
        }
        0.5 * s
    };
    let want: f64 = rings.iter().map(|r| area(r)).sum();
    let got: f64 = pieces.iter().map(|p| area(p)).sum();
    assert!(
        (got - want).abs() < 1e-9,
        "piece area {got} vs polygon {want} ({pieces:?})"
    );
    pieces
}

/// **The case this module exists for, and the one CAD actually produces**: every
/// edge axis-aligned, so the sweep coordinate is shared by pairs of vertices all
/// over. If the total order were not total, this is where it would show.
#[test]
fn a_rectangle_with_a_rectangular_hole() {
    let pieces = check(
        &[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 3.0],
            [0.0, 3.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [3.0, 2.0],
            [3.0, 1.0],
        ],
        &[&[0, 1, 2, 3], &[4, 5, 6, 7]],
    );
    assert_eq!(pieces.len(), 2, "one hole splits it in two: {pieces:?}");
}

/// A convex ring is already monotone, so the sweep must add nothing.
#[test]
fn a_square_is_left_alone() {
    let pieces = check(
        &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        &[&[0, 1, 2, 3]],
    );
    assert_eq!(pieces.len(), 1);
}

/// One split and one merge vertex, and no holes — the notch supplies both.
#[test]
fn a_u_shape_needs_its_diagonals() {
    check(
        &[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [2.0, 3.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 3.0],
            [0.0, 3.0],
        ],
        &[&(0..8).collect::<Vec<_>>()],
    );
}

/// **The measured failure**, transcribed: the two holes whose bridges collided on
/// one outer corner and left ear clipping with no ear at all.
#[test]
fn the_wall_with_two_windows() {
    let pieces = check(
        &[
            [0.0, 1.0],
            [0.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            [1.0, 0.25959136597258015],
            [1.0, 0.7035578716424774],
            [2.0, 0.7035578716424774],
            [2.0, 0.25959136597258015],
            [1.0, -0.703557871642477],
            [1.0, -0.25959136597258003],
            [2.0, -0.25959136597258003],
            [2.0, -0.703557871642477],
        ],
        &[&[0, 1, 2, 3], &[7, 6, 5, 4], &[11, 10, 9, 8]],
    );
    assert!(pieces.len() >= 3, "two holes need at least two cuts");
}

/// **The corpus that measured the old triangulator's 27.6% failure rate**, pointed
/// at the decomposition instead: random axis-aligned rectangles with one to four
/// non-overlapping rectangular holes on an integer grid.
///
/// Hand-picked fixtures only cover the degeneracies someone thought of. This one
/// produces shared sweep coordinates, horizontal edges and collinear vertices by
/// the thousand, which is precisely where a sweep goes wrong — and every sample
/// must come back as monotone pieces whose areas add up.
#[test]
fn random_rectilinear_faces_decompose() {
    let mut st = 0x5EED_1234_ABCD_0001u64;
    let mut lcg = move || {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        st >> 33
    };
    let mut rng = move |lo: i64, hi: i64| lo + (lcg() % ((hi - lo + 1) as u64)) as i64;
    let (w, h) = (24i64, 16i64);
    let mut samples = 0;
    for _ in 0..3000 {
        let k = rng(1, 4);
        let mut boxes: Vec<[i64; 4]> = Vec::new();
        for _ in 0..k {
            for _try in 0..20 {
                let x0 = rng(1, w - 3);
                let y0 = rng(1, h - 3);
                let x1 = x0 + rng(1, (w - 1 - x0).min(5));
                let y1 = y0 + rng(1, (h - 1 - y0).min(5));
                if boxes
                    .iter()
                    .any(|b| x0 <= b[2] && b[0] <= x1 && y0 <= b[3] && b[1] <= y1)
                {
                    continue;
                }
                boxes.push([x0, y0, x1, y1]);
                break;
            }
        }
        if boxes.is_empty() {
            continue;
        }
        let mut pts = vec![
            [0.0, 0.0],
            [w as f64, 0.0],
            [w as f64, h as f64],
            [0.0, h as f64],
        ];
        let mut holes: Vec<Vec<usize>> = Vec::new();
        for b in &boxes {
            let base = pts.len();
            // Clockwise, so the interior stays on every edge's left.
            for c in [[b[0], b[1]], [b[0], b[3]], [b[2], b[3]], [b[2], b[1]]] {
                pts.push([c[0] as f64, c[1] as f64]);
            }
            holes.push((base..base + 4).collect());
        }
        let mut rings: Vec<&[usize]> = vec![&[0, 1, 2, 3]];
        rings.extend(holes.iter().map(|h| h.as_slice()));
        check(&pts, &rings);
        samples += 1;
    }
    assert!(samples > 2500, "only {samples} samples reached the gate");
}

/// A vertex used by two rings is not a polygon with sibling holes — it is a pinch,
/// and the sweep would have no way to know which chain it was on.
#[test]
fn a_shared_vertex_is_rejected() {
    let pts = [
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 3.0],
        [0.0, 3.0],
        [1.0, 1.0],
        [1.0, 2.0],
    ];
    assert!(matches!(
        decompose(&uv(&pts), &[&[0, 1, 2, 3], &[4, 5, 0]], None),
        Err(TessError::DegenerateRing)
    ));
}
