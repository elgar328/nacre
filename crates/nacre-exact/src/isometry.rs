use super::*;
/// A coordinate axis — the fixed axis of an axis-aligned rotation. Rotations turn about
/// `X`/`Y`/`Z` only; an arbitrary rational axis would take Rodrigues' formula.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    /// The coordinate index this axis names (`X → 0`, `Y → 1`, `Z → 2`) — the one a reflection in
    /// a plane perpendicular to it negates.
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }

    /// The two in-plane coordinate indices (the third is the fixed rotation axis).
    /// The order gives a right-handed (CCW-about-the-axis) rotation.
    pub fn plane(self) -> (usize, usize) {
        match self {
            Axis::X => (1, 2), // rotate y,z
            Axis::Y => (2, 0), // rotate z,x
            Axis::Z => (0, 1), // rotate x,y
        }
    }
}

/// Whether **any** rotation about an axis-parallel line fixes the plane
/// `a·x + b·y + c·z + d = 0` as a set: the normal rides the axis (its other two
/// components exactly zero), so the normal coordinate of every point is untouched
/// whatever the pivot or angle. One of the three set-invariance atoms — the others
/// are [`translation_fixes_plane`] and [`mirror_fixes_plane`]; `Isometry::fixes_plane`
/// composes the first two, and a recorded motion chain is checked node by node.
pub fn axis_rotation_fixes_plane(axis: Axis, coeffs: &[Rat; 4]) -> bool {
    let i = axis.index();
    let zero = Rat::from_int(0);
    coeffs[i] != zero && (0..3).all(|k| k == i || coeffs[k] == zero)
}

/// Whether the translation `offset` fixes the plane as a set: its component along the
/// normal is exactly zero (`n · offset == 0`, checked — an overflow answers `false`,
/// a conservative miss).
pub fn translation_fixes_plane(offset: &[Rat; 3], coeffs: &[Rat; 4]) -> bool {
    let zero = Rat::from_int(0);
    let dot = coeffs[..3]
        .iter()
        .zip(offset)
        .try_fold(zero, |acc, (&a, &t)| a.checked_mul(t)?.checked_add(acc));
    dot == Some(zero)
}

/// Whether the reflection in the coordinate plane `axis = offset` fixes the plane as a
/// set: either the normal is perpendicular to the mirror axis (the plane contains the
/// mirrored direction, so the set maps to itself), or the plane **is** the mirror plane
/// (`n ∥ axis` and `n_axis · offset + d == 0`, checked — overflow answers `false`).
/// Set-fixed only: the second case flips the normal's sense, which a canonical name
/// does not carry.
pub fn mirror_fixes_plane(axis: Axis, offset: Rat, coeffs: &[Rat; 4]) -> bool {
    let i = axis.index();
    let zero = Rat::from_int(0);
    if coeffs[i] == zero {
        return true;
    }
    (0..3).all(|k| k == i || coeffs[k] == zero)
        && coeffs[i]
            .checked_mul(offset)
            .and_then(|p| p.checked_add(coeffs[3]))
            == Some(zero)
}

/// An axis-aligned rigid rotation: turn about `axis` (the line through the rational
/// `pivot`) by the rational `angle`. Exact for the 90°-family (`try_exact_cos_sin`);
/// otherwise the realized coordinate is irrational (cos/sin) and carries tol.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation {
    pub axis: Axis,
    pub pivot: [Rat; 3],
    pub angle: Angle,
}

/// A rigid-body isometry: a rotation (optional) then a translation.
/// The exact rational data is the **definition**; the `apply_*`/`offset_f64`
/// realizers give the f64 cache. Math-type independent — operates on plain
/// `[f64; 3]`, so `nacre-exact` never depends on `nacre-math`; the caller
/// (`nacre-ops`) applies it to `Point3`/`Plane`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Isometry {
    /// Applied first: an axis-aligned rotation, or `None` (pure translation).
    pub rotate: Option<Rotation>,
    /// Applied second: an exact rational translation.
    pub translate: [Rat; 3],
}

impl Isometry {
    /// A pure translation by the rational vector `translate`.
    pub fn translation(translate: [Rat; 3]) -> Self {
        Isometry {
            rotate: None,
            translate,
        }
    }

    /// A pure rotation (no translation).
    pub fn rotation(rotate: Rotation) -> Self {
        Isometry {
            rotate: Some(rotate),
            translate: [Rat::from_int(0); 3],
        }
    }

    /// A rotation followed by a translation.
    pub fn rigid(rotate: Rotation, translate: [Rat; 3]) -> Self {
        Isometry {
            rotate: Some(rotate),
            translate,
        }
    }

    /// The translation realized in f64.
    pub fn offset_f64(&self) -> [f64; 3] {
        [
            self.translate[0].to_f64(),
            self.translate[1].to_f64(),
            self.translate[2].to_f64(),
        ]
    }

    /// Whether the isometry realizes exactly: no rotation, or a 90°-family rotation
    /// (`try_exact_cos_sin` gives rational cos/sin, so an f64-representable point
    /// stays exact — tol 0). A non-90° rotation realizes to irrational f64 (tol > 0).
    pub fn is_exact(&self) -> bool {
        match self.rotate {
            None => true,
            Some(r) => r.angle.try_exact_cos_sin().is_some(),
        }
    }

    /// Whether this isometry maps the plane `a·x + b·y + c·z + d = 0` onto itself
    /// **as a set** — exactly, on the rational definition.
    ///
    /// A rotation about an axis parallel to the plane's normal permutes the plane
    /// within itself whatever the pivot (the normal coordinate is untouched), and a
    /// translation moves the plane iff its component along the normal is nonzero.
    /// So the condition is: the rotation, if any, has its axis parallel to
    /// `(a, b, c)`, and `(a, b, c) · translate == 0`. `d` plays no part — every
    /// plane sharing the normal answers alike. Checked arithmetic; an overflow
    /// answers `false`, a conservative miss (the plane is then carried on the
    /// recorded path, slower but never wrong).
    ///
    /// A reflection is not an `Isometry`; [`mirror_fixes_plane`] answers for the mirror
    /// motion, and a recorded chain is checked node by node with the same three atoms.
    pub fn fixes_plane(&self, coeffs: &[Rat; 4]) -> bool {
        self.rotate
            .is_none_or(|r| axis_rotation_fixes_plane(r.axis, coeffs))
            && translation_fixes_plane(&self.translate, coeffs)
    }

    /// Apply the full isometry (rotate about the axis point, then translate) to a
    /// point realized in f64.
    pub fn apply_point(&self, p: [f64; 3]) -> [f64; 3] {
        let mut q = p;
        if let Some(r) = self.rotate {
            let (i, j) = r.axis.plane();
            let (px, py) = (r.pivot[i].to_f64(), r.pivot[j].to_f64());
            let (c, s) = r.angle.cos_sin_f64();
            let (dx, dy) = (p[i] - px, p[j] - py);
            q[i] = px + dx * c - dy * s;
            q[j] = py + dx * s + dy * c;
        }
        let off = self.offset_f64();
        [q[0] + off[0], q[1] + off[1], q[2] + off[2]]
    }

    /// Apply only the rotation (no axis point, no translation) to a direction.
    pub fn apply_dir(&self, d: [f64; 3]) -> [f64; 3] {
        let mut q = d;
        if let Some(r) = self.rotate {
            let (i, j) = r.axis.plane();
            let (c, s) = r.angle.cos_sin_f64();
            let (dx, dy) = (d[i], d[j]);
            q[i] = dx * c - dy * s;
            q[j] = dx * s + dy * c;
        }
        q
    }

    /// **This isometry applied to a rational plane**, exactly — the definition-level twin of
    /// [`apply_point`](Isometry::apply_point).
    ///
    /// A plane is not a bag of points, so it does not go through `apply_point`: under
    /// `x ↦ R(x − p) + p + t` the plane `n·x + d = 0` becomes
    ///
    /// ```text
    /// n' = R·n            d' = d + n·p − n'·(p + t)
    /// ```
    ///
    /// (For a pure translation that collapses to the familiar `d' = d − n·t`.)
    ///
    /// `None` unless the rotation is one the rationals can state — the 90°-family, where
    /// [`Angle::try_exact_cos_sin`] gives `cos`/`sin` in `{0, ±1}` — or on `i128` overflow.
    /// That is the same condition [`Isometry::is_exact`] reports, so a caller that already
    /// checked it will not be surprised here.
    ///
    /// The result is canonicalized ([`canonical_plane_coeffs`]), so a plane reached by two
    /// different routes lands on the **same array**.
    pub fn plane_coeffs(&self, c: [Rat; 4]) -> Option<[Rat; 4]> {
        let n = [c[0], c[1], c[2]];
        let (n2, pivot) = match self.rotate {
            None => (n, [Rat::from_int(0); 3]),
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                let mut m = n;
                m[i] = n[i].checked_mul(cos)?.checked_sub(n[j].checked_mul(sin)?)?;
                m[j] = n[i].checked_mul(sin)?.checked_add(n[j].checked_mul(cos)?)?;
                (m, r.pivot)
            }
        };
        // d' = d + n·p − n'·(p + t)
        let mut d = c[3];
        for k in 0..3 {
            d = d.checked_add(n[k].checked_mul(pivot[k])?)?;
            let q = pivot[k].checked_add(self.translate[k])?;
            d = d.checked_sub(n2[k].checked_mul(q)?)?;
        }
        canonical_plane_coeffs([n2[0], n2[1], n2[2], d])
    }

    /// **This isometry applied to a rational point**, exactly — the rational twin
    /// [`apply_point`](Isometry::apply_point) never had, and the map
    /// [`plane_coeffs`](Isometry::plane_coeffs) is described against.
    ///
    /// `x ↦ R(x − p) + p + t`, with `R` read from [`Angle::try_exact_cos_sin`]. `None` under the
    /// same conditions as `plane_coeffs` — a rotation the rationals cannot state, or `i128`
    /// overflow — so a caller that can move one description exactly can move the other.
    ///
    /// ★ **Why a plane needs it even though a plane is not a bag of points.** A plane's *truth* is
    /// three points on it (`Model::surface_points`): its coefficients are a product of two point
    /// differences and overflow `i128` far sooner than the points do. Moving such a plane means
    /// moving its points, and carrying them through f64 would put a rounded coordinate back into a
    /// definition — the thing storing points was meant to stop.
    pub fn point_rat(&self, x: [Rat; 3]) -> Option<[Rat; 3]> {
        let turned = match self.rotate {
            None => x,
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                // Pivot-relative, exactly as the realization does it: `u = x − p`, turn, shift back.
                let u = x[i].checked_sub(r.pivot[i])?;
                let v = x[j].checked_sub(r.pivot[j])?;
                let mut m = x;
                m[i] = r.pivot[i]
                    .checked_add(u.checked_mul(cos)?.checked_sub(v.checked_mul(sin)?)?)?;
                m[j] = r.pivot[j]
                    .checked_add(u.checked_mul(sin)?.checked_add(v.checked_mul(cos)?)?)?;
                m
            }
        };
        Some([
            turned[0].checked_add(self.translate[0])?,
            turned[1].checked_add(self.translate[1])?,
            turned[2].checked_add(self.translate[2])?,
        ])
    }

    /// **This isometry's rotation applied to a rational direction**, exactly — the rational twin
    /// of [`apply_dir`](Isometry::apply_dir). A direction is a difference of points, so the pivot
    /// and the translation cancel and only the turn remains. `None` under the same conditions as
    /// [`point_rat`](Isometry::point_rat) — a rotation the rationals cannot state, or `i128`
    /// overflow.
    pub fn dir_rat(&self, d: [Rat; 3]) -> Option<[Rat; 3]> {
        match self.rotate {
            None => Some(d),
            Some(r) => {
                let (cos, sin) = r.angle.try_exact_cos_sin()?;
                let (i, j) = r.axis.plane();
                let mut m = d;
                m[i] = d[i].checked_mul(cos)?.checked_sub(d[j].checked_mul(sin)?)?;
                m[j] = d[i].checked_mul(sin)?.checked_add(d[j].checked_mul(cos)?)?;
                Some(m)
            }
        }
    }
}

/// **A reflection in `axis = offset` applied to a rational point**, exactly — `x_a ↦ 2·offset − x_a`.
///
/// The point twin of [`mirror_plane_coeffs`], and the same one-line map: a reflection is its own
/// inverse and touches one coordinate. `None` only on `i128` overflow.
pub fn mirror_point_rat(x: [Rat; 3], axis: Axis, offset: Rat) -> Option<[Rat; 3]> {
    let a = axis.index();
    let mut out = x;
    out[a] = offset.checked_mul(Rat::from_int(2))?.checked_sub(x[a])?;
    Some(out)
}

/// **A reflection in `axis = offset` applied to a rational plane**, exactly.
///
/// `x_a ↦ 2·offset − x_a` negates that component of the normal and shifts the offset:
/// `n'_a = −n_a`, `d' = d + 2·offset·n_a`. The reflection is its own inverse, which is why the
/// map and its transpose-inverse coincide and no case analysis is needed.
///
/// ★ **The determinant is `−1`.** A caller that relies on orientation being preserved has to
/// account for that itself; this function states where the plane goes, nothing more.
pub fn mirror_plane_coeffs(c: [Rat; 4], axis: Axis, offset: Rat) -> Option<[Rat; 4]> {
    let a = axis.index();
    let mut out = c;
    out[a] = Rat::from_int(0).checked_sub(c[a])?;
    out[3] = c[3].checked_add(offset.checked_mul(Rat::from_int(2))?.checked_mul(c[a])?)?;
    canonical_plane_coeffs(out)
}
