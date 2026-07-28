//! **Constrained Delaunay by edge flipping** (Lawson 1977; Chew 1989).
//!
//! The decomposition next door produces a *correct* triangulation, not a pretty one:
//! a y-monotone piece leaves the stack algorithm almost no freedom, so its output is
//! full of slivers. Quality is a separate question with a textbook answer — among
//! **all** triangulations respecting a given set of constrained edges, the constrained
//! Delaunay one maximises the smallest angle, and any triangulation reaches it by
//! flipping the free edges that fail the incircle test, in finitely many steps.
//!
//! So correctness and quality stay orthogonal: [`super::monotone`] settles the first
//! and this settles the second, and deleting this file would leave a mesh that is still
//! right. Nothing here creates or removes a vertex or a triangle, and no constrained
//! edge is ever touched, so watertightness, the triangle count and the boundary all
//! come through untouched by construction.
//!
//! **★ Measured, because the obvious expectation was wrong.** The mean smallest angle
//! roughly doubles — on a cylinder cap at `tol` 1e-4, `0.51° → 1.05°`; at 1e-2,
//! `5.04° → 7.12°` — but the *worst* triangle on a cap does not move at all. A cap's
//! vertices are a regular polygon, and a regular polygon's vertices are **cocircular**:
//! every triangulation of them is Delaunay, `incircle` answers zero, and there is
//! nothing to flip. What flips at all there is only the residue of the arc having been
//! sampled in `f64`. Getting a fat triangle out of a 314-gon needs a vertex *inside*
//! it, and adding vertices is what the crack-free contract forbids at this layer
//! (design §5) — so Delaunay refinement is not the missing step here, it is a
//! different layer's decision.

use super::P2;
use nacre_predicates::{incircle, orient2d};
use std::collections::{HashMap, HashSet};

/// Flip every free edge that fails the incircle test, until none does.
///
/// `constrained` holds the ring edges as `(min, max)` pairs — the polygon's boundary
/// and its holes' — and they are the edges a flip may not take. Triangles must arrive
/// counter-clockwise and leave the same way.
pub(super) fn refine(uv: &[P2], tris: &mut [[usize; 3]], constrained: &HashSet<(usize, usize)>) {
    // Which triangle contains each *directed* edge. In a triangulation every directed
    // edge belongs to exactly one triangle, so this is a function — and the neighbour
    // across `a → b` is simply whoever owns `b → a`.
    let mut owner: HashMap<(usize, usize), usize> = HashMap::new();
    for (i, t) in tris.iter().enumerate() {
        for k in 0..3 {
            owner.insert((t[k], t[(k + 1) % 3]), i);
        }
    }

    let key = |a: usize, b: usize| (a.min(b), a.max(b));
    let mut queued: HashSet<(usize, usize)> = HashSet::new();
    let mut stack: Vec<(usize, usize)> = Vec::new();
    for t in tris.iter() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if !constrained.contains(&key(a, b)) && queued.insert(key(a, b)) {
                stack.push((a, b));
            }
        }
    }

    // Lawson terminates, but a bound costs one counter and turns "we believe it does"
    // into "it did" — a mis-signed predicate would otherwise spin here forever.
    let mut budget = 64 * tris.len() + 64;
    while let Some((a, b)) = stack.pop() {
        queued.remove(&key(a, b));
        if budget == 0 {
            return;
        }
        budget -= 1;
        let (Some(&t1), Some(&t2)) = (owner.get(&(a, b)), owner.get(&(b, a))) else {
            continue; // a boundary edge, or one this loop has already replaced
        };
        let third = |t: [usize; 3], x: usize, y: usize| {
            t.into_iter()
                .find(|&v| v != x && v != y)
                .expect("a triangle")
        };
        let (c, d) = (third(tris[t1], a, b), third(tris[t2], a, b));

        // The quad `a → d → b → c` must be convex, or the flip would fold a triangle
        // over. Asking whether the *replacements* are properly wound is the same
        // question and needs no separate case analysis.
        if orient2d(uv[a], uv[d], uv[c]) <= 0.0 || orient2d(uv[d], uv[b], uv[c]) <= 0.0 {
            continue;
        }
        // `(a, b, c)` is counter-clockwise because `t1` is, so this asks exactly what
        // the Delaunay condition asks: is the opposite vertex inside the circumcircle?
        if incircle(uv[a], uv[b], uv[c], uv[d]) <= 0.0 {
            continue;
        }

        for t in [t1, t2] {
            for k in 0..3 {
                owner.remove(&(tris[t][k], tris[t][(k + 1) % 3]));
            }
        }
        tris[t1] = [a, d, c];
        tris[t2] = [d, b, c];
        for t in [t1, t2] {
            for k in 0..3 {
                owner.insert((tris[t][k], tris[t][(k + 1) % 3]), t);
            }
        }
        // The quad's four sides may have stopped being Delaunay now that the diagonal
        // moved; the new diagonal `c–d` cannot, since it was just chosen.
        for (x, y) in [(a, d), (d, b), (b, c), (c, a)] {
            if !constrained.contains(&key(x, y)) && queued.insert(key(x, y)) {
                stack.push((x, y));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The smallest angle anywhere in the mesh, in degrees — the quantity constrained
    /// Delaunay maximises, and therefore the one worth measuring.
    pub(in crate::polygon) fn min_angle(uv: &[P2], tris: &[[usize; 3]]) -> f64 {
        let mut worst = 180.0f64;
        for t in tris {
            for k in 0..3 {
                let (o, p, q) = (uv[t[k]], uv[t[(k + 1) % 3]], uv[t[(k + 2) % 3]]);
                let (u, v) = ([p[0] - o[0], p[1] - o[1]], [q[0] - o[0], q[1] - o[1]]);
                let dot = u[0] * v[0] + u[1] * v[1];
                let cross = u[0] * v[1] - u[1] * v[0];
                worst = worst.min(cross.abs().atan2(dot).to_degrees());
            }
        }
        worst
    }

    /// A long strip with vertices down both sides: fan it from one corner and the
    /// triangles are splinters, flip it and they are near-right-angled halves of the
    /// squares they sit in. Delaunay maximises the smallest angle, and here that is a
    /// large, visible number rather than a technicality.
    ///
    /// **The obvious fixture — a regular n-gon — would have proved nothing**, and that
    /// is worth knowing rather than discovering twice: its vertices are *cocircular*,
    /// so every triangulation of them is Delaunay, `incircle` answers exactly zero at
    /// every candidate, and not one flip fires. A cylinder cap is precisely that shape,
    /// which is why the flip pass improves the caps' *average* triangle (their sampled
    /// coordinates are only approximately cocircular) but cannot touch their worst one.
    #[test]
    fn a_fan_of_splinters_becomes_a_delaunay_strip() {
        let k = 12;
        let mut uv: Vec<P2> = (0..=k).map(|i| [i as f64, 0.0]).collect();
        uv.extend((0..=k).rev().map(|i| [i as f64, 1.0]));
        let n = uv.len();
        // The worst legal triangulation: every triangle hangs off vertex 0.
        let mut tris: Vec<[usize; 3]> = (1..n - 1).map(|i| [0, i, i + 1]).collect();
        let constrained: HashSet<(usize, usize)> = (0..n)
            .map(|i| (i.min((i + 1) % n), i.max((i + 1) % n)))
            .collect();

        let before = min_angle(&uv, &tris);
        let count = tris.len();
        refine(&uv, &mut tris, &constrained);
        let after = min_angle(&uv, &tris);

        assert_eq!(
            tris.len(),
            count,
            "a flip neither adds nor removes a triangle"
        );
        assert!(
            tris.iter()
                .all(|t| orient2d(uv[t[0]], uv[t[1]], uv[t[2]]) > 0.0),
            "every triangle is still counter-clockwise"
        );
        assert!(
            after > 40.0 && after > 5.0 * before,
            "min angle {before:.2}° → {after:.2}°"
        );
    }

    /// The regular polygon, kept as a *negative* result: cocircular points leave the
    /// Delaunay condition with nothing to say, so the mesh must come back untouched.
    /// Nothing is wrong when that happens — it is what the theorem says.
    #[test]
    fn cocircular_points_admit_every_triangulation() {
        // Every integer point on `x² + y² = 25²`, counter-clockwise. Sampling a circle
        // with `cos`/`sin` would *not* do — the rounding leaves the points only nearly
        // cocircular, which is exactly why a cylinder cap does see some flips.
        let uv: Vec<P2> = [
            [25, 0],
            [24, 7],
            [20, 15],
            [15, 20],
            [7, 24],
            [0, 25],
            [-7, 24],
            [-15, 20],
            [-20, 15],
            [-24, 7],
            [-25, 0],
            [-24, -7],
            [-20, -15],
            [-15, -20],
            [-7, -24],
            [0, -25],
            [7, -24],
            [15, -20],
            [20, -15],
            [24, -7],
        ]
        .into_iter()
        .map(|p: [i32; 2]| [p[0] as f64, p[1] as f64])
        .collect();
        let n = uv.len();
        let before: Vec<[usize; 3]> = (1..n - 1).map(|i| [0, i, i + 1]).collect();
        let mut tris = before.clone();
        let constrained: HashSet<(usize, usize)> = (0..n)
            .map(|i| (i.min((i + 1) % n), i.max((i + 1) % n)))
            .collect();
        refine(&uv, &mut tris, &constrained);
        assert_eq!(
            tris, before,
            "a flip fired on points the incircle test cannot separate"
        );
    }

    /// Flipping must not touch the boundary, or the mesh would stop matching the face.
    #[test]
    fn constrained_edges_survive() {
        let uv: Vec<P2> = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [0.0, 1.0]];
        let mut tris = vec![[0, 1, 2], [0, 2, 3]];
        let constrained: HashSet<(usize, usize)> =
            [(0, 1), (1, 2), (2, 3), (0, 3)].into_iter().collect();
        refine(&uv, &mut tris, &constrained);
        let edges: HashSet<(usize, usize)> = tris
            .iter()
            .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
            .collect();
        for e in &constrained {
            assert!(edges.contains(e), "constrained edge {e:?} was flipped away");
        }
    }
}
