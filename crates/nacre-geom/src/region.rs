//! Area and centroid of a planar region bounded by 3-D polygons.
//!
//! Pure geometry — a list of points in, two numbers out; nothing here knows about topology. It
//! lives in `nacre-geom` for the same reason the ring predicates do: two crates above need it
//! (`nacre-props` reports it, `nacre-ops` places a sketch frame with it), and `nacre-props` is a
//! *dev*-dependency of `nacre-ops`, so a shared helper cannot live there.

use nacre_math::{Point3, Vector3};

/// Area and area-weighted centroid of a planar region: one outer ring minus its holes.
///
/// The rings are closed implicitly and must be **simple and properly nested** (the profile
/// contract, which `Profile2d::check` enforces on authored input). Returns `None` if the region
/// encloses nothing — a ring of fewer than three points, or an outer ring whose area the holes
/// cancel exactly.
///
/// **Not a vertex average.** Each triangle of a fan from the first vertex is weighted by its
/// *signed* area about the ring's own normal, so a concave outline comes out right where a plain
/// mean of the corners would not. The distinction matters beyond accuracy: the area centroid is a
/// property of the *region*, so inserting extra vertices along a straight edge does not move it,
/// while a vertex average shifts with every such change.
pub fn planar_region_area_centroid(outer: &[Point3], holes: &[&[Point3]]) -> Option<(f64, Point3)> {
    let (mut area, c_out) = ring_area_centroid(outer)?;
    // Area and first moment are both additive, so a hole subtracts each. Accumulate the moment
    // about `c_out` to keep the numbers small and independent of the model's placement.
    let mut moment = Vector3::zero();
    for hole in holes {
        let (a_in, c_in) = ring_area_centroid(hole)?;
        area -= a_in;
        moment -= (c_in - c_out) * a_in;
    }
    if area == 0.0 {
        return None;
    }
    Some((area, c_out + moment * (1.0 / area)))
}

/// Area and area-weighted centroid of one simple planar ring, via a signed triangle fan from its
/// first vertex.
fn ring_area_centroid(points: &[Point3]) -> Option<(f64, Point3)> {
    if points.len() < 3 {
        return None;
    }
    let base = points[0];
    // Area vector A_vec = ½ Σ (vᵢ − v₀) × (vᵢ₊₁ − v₀); |A_vec| is the true area.
    let mut area_vec = Vector3::zero();
    for pair in points[1..].windows(2) {
        area_vec += (pair[0] - base).cross(pair[1] - base);
    }
    if area_vec == Vector3::zero() {
        return None;
    }

    let mut weighted = Vector3::zero();
    let mut weight = 0.0;
    for pair in points[1..].windows(2) {
        // ★ Weighted by the **un-normalized** area vector. `normalize` divides by a square root,
        // so every `signed` carried that rounding; the centroid is a *ratio* of these weights, so
        // the common `|area_vec|` factor cancels and the normalization only ever added error.
        let signed = (pair[0] - base).cross(pair[1] - base).dot(area_vec);
        // ★ The `1/3` is **not** applied here. `1.0 / 3.0` is not an exact `f64`, so folding it
        // into every triangle rounds once per triangle and the terms of a symmetric ring stop
        // cancelling — a rectangle centred on the origin came out at `−1.1e-16` instead of `0`,
        // which then moved a sketch placed in that frame by one ulp. Factored out, the loop is
        // exact for exactly-representable rings and one division does all the rounding.
        weighted += ((pair[0] - base) + (pair[1] - base)) * signed;
        weight += signed;
    }
    if weight == 0.0 {
        return None;
    }
    let d = 3.0 * weight;
    let [wx, wy, wz] = weighted.as_array();
    Some((
        0.5 * area_vec.norm(),
        base + Vector3::from_array([wx / d, wy / d, wz / d]),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ★★★ **A ring symmetric about the origin has its centroid exactly at the origin.**
    ///
    /// It did not. `ring_area_centroid` folded `1.0 / 3.0` — not an exact `f64` — into every
    /// triangle of its fan, and weighted each by a dot with the **normalized** area vector, so a
    /// square root's rounding rode along too. The terms of a symmetric ring then failed to cancel
    /// and the answer came out one ulp off zero.
    ///
    /// That is not cosmetic. `face_frame` uses this point as a sketch frame's **origin**, so a
    /// profile placed on such a face lands one ulp away from where the same profile landed on the
    /// face below it — and padding a footprint twice left two faces of area `2.2e-16` in the
    /// result, which `validate` did not report.
    ///
    /// The `1/3` is factored out and the weights use the un-normalized area vector: the centroid
    /// is a *ratio* of those weights, so the common `|area_vec|` cancels and normalizing only ever
    /// added error.
    ///
    /// ★★★★ **Exactness here is not general and this test does not claim it is.** Removing two
    /// roundings makes the answer exact whenever what remains is exactly representable — as it is
    /// for the fixture that exposed the defect. `±0.1 × ±7.7` still lands one ulp off, because the
    /// weighted sum itself rounds. Closing that needs the centroid derived in `Rat` from the face's
    /// own exact vertices, not computed in `f64` and lifted.
    #[test]
    fn a_symmetric_ring_centroids_exactly_on_its_centre() {
        let p = |x: f64, y: f64, z: f64| Point3::from_array([x, y, z]);
        // (0.6, 1.0, 2.1) is the footprint and height that left the zero-area faces.
        for (hx, hy, z) in [(0.6, 1.0, 2.1), (1.5, 1.5, 0.0), (0.25, 4.0, -3.5)] {
            let ring = [p(-hx, -hy, z), p(hx, -hy, z), p(hx, hy, z), p(-hx, hy, z)];
            let (area, c) = planar_region_area_centroid(&ring, &[]).expect("a real rectangle");
            assert_eq!(
                c.as_array(),
                [0.0, 0.0, z],
                "rectangle ±{hx}×±{hy} at z={z} centroids off-centre (area {area})"
            );
        }
    }

    /// The same, shifted off the origin: the centroid must be exactly the rectangle's own centre,
    /// not merely close to it. Pins that the fix is not an accident of zeros cancelling.
    #[test]
    fn a_shifted_rectangle_centroids_exactly_on_its_centre() {
        let p = |x: f64, y: f64| Point3::from_array([x, y, 0.0]);
        let (x0, x1, y0, y1) = (1.25, 4.75, -2.5, 0.5);
        let ring = [p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)];
        let (_, c) = planar_region_area_centroid(&ring, &[]).expect("a real rectangle");
        assert_eq!(c.as_array(), [(x0 + x1) / 2.0, (y0 + y1) / 2.0, 0.0]);
    }

    fn pts(v: &[[f64; 3]]) -> Vec<Point3> {
        v.iter().map(|p| Point3::from_array(*p)).collect()
    }

    fn close(a: Point3, b: Point3) -> bool {
        (a - b).norm() < 1e-12
    }

    #[test]
    fn a_square_centres_on_its_middle() {
        let sq = pts(&[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ]);
        let (area, c) = planar_region_area_centroid(&sq, &[]).unwrap();
        assert!((area - 16.0).abs() < 1e-12);
        assert!(close(c, Point3::from_array([2.0, 2.0, 0.0])));
    }

    /// The case a vertex average gets wrong: an L has more corners on the thin arm, so their mean
    /// is pulled towards it while the area centroid stays with the material.
    #[test]
    fn a_reflex_outline_is_not_the_vertex_average() {
        // L = [0,4]×[0,4] minus the [2,4]×[2,4] corner. Area 12, centroid (5/3, 5/3).
        let l = pts(&[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [2.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ]);
        let (area, c) = planar_region_area_centroid(&l, &[]).unwrap();
        assert!((area - 12.0).abs() < 1e-12);
        assert!(
            close(c, Point3::from_array([5.0 / 3.0, 5.0 / 3.0, 0.0])),
            "{c:?}"
        );
        // The mean of the six corners is (2, 2) — a different point.
        let mean = Point3::centroid(&l).unwrap();
        assert!((mean - c).norm() > 0.3, "the two rules must disagree here");
    }

    /// **The property the sketch frame depends on**: the area centroid is a property of the
    /// region, so subdividing an edge leaves it where it was. A vertex average does not.
    #[test]
    fn extra_collinear_vertices_do_not_move_the_centroid() {
        let plain = pts(&[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [2.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ]);
        let subdivided = pts(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], // mid-run on a straight edge
            [3.0, 0.0, 0.0], // and another
            [4.0, 0.0, 0.0],
            [4.0, 2.0, 0.0],
            [2.0, 2.0, 0.0],
            [2.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ]);
        let (a0, c0) = planar_region_area_centroid(&plain, &[]).unwrap();
        let (a1, c1) = planar_region_area_centroid(&subdivided, &[]).unwrap();
        assert!((a0 - a1).abs() < 1e-12);
        assert!((c0 - c1).norm() < 1e-12, "{c0:?} vs {c1:?}");
        // The rule this replaces moves by a large fraction of the face instead.
        let m0 = Point3::centroid(&plain).unwrap();
        let m1 = Point3::centroid(&subdivided).unwrap();
        assert!(
            (m0 - m1).norm() > 0.1,
            "vertex average must move: {m0:?} vs {m1:?}"
        );
    }

    #[test]
    fn an_off_centre_hole_pulls_the_centroid_away() {
        // [0,10]² with a [1,3]² hole near the low corner: area 100 − 4 = 96.
        let outer = pts(&[
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
        ]);
        let hole = pts(&[
            [1.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 3.0, 0.0],
            [1.0, 3.0, 0.0],
        ]);
        let (area, c) = planar_region_area_centroid(&outer, &[&hole]).unwrap();
        assert!((area - 96.0).abs() < 1e-12);
        // (100·5 − 4·2)/96 on each axis.
        let want = (100.0 * 5.0 - 4.0 * 2.0) / 96.0;
        assert!(close(c, Point3::from_array([want, want, 0.0])), "{c:?}");
    }

    #[test]
    fn a_region_that_encloses_nothing_is_none() {
        let two = pts(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert!(planar_region_area_centroid(&two, &[]).is_none());
        let collinear = pts(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        assert!(planar_region_area_centroid(&collinear, &[]).is_none());
    }
}
