//! **Solids built through the product's own operations** — test-only (`test-util`).
//!
//! ★ A fixture here states nothing the product does not: it is a sequence of calls an application
//! makes (a datum plane, a circle sketched on it, an extrude). That is the difference from a
//! constructor that writes cells directly — such a door skips the realization funnel, and a
//! population built through it can hold what no operation produces (a tilted cap lifted from
//! realized `f64`, whose rim is an ellipse while its edge cache says circle).
//!
//! ★ **A box is a floor and a height**, the way an application states one (kit's `cuboid`
//! extrudes a rectangle by its size): a datum plane, a rectangle on it, an extrude. Two corners
//! are not a statement the product takes — the extrude reads its distance as the shortest decimal
//! that round-trips (design 「구성 산술은 사용자가 쓴 십진수로 한다」), so a box whose two heights
//! differ by no such decimal cannot be built with its far cap where the corners say.

use super::*;

/// A cylinder built by [`cylinder`] or [`cylinder_with_seam`], with its faces **named** — the
/// extrude's `[base, top, sides…]` order is a fact of `OpOutput::Extrude`, and a caller reading
/// a face by index would be guessing it.
#[derive(Clone, Copy, Debug)]
pub struct CylinderSolid {
    pub solid: Handle<Solid>,
    /// The cap on the sketch plane, facing against the axis.
    pub base: Handle<Face>,
    /// The far cap, facing along the axis.
    pub top: Handle<Face>,
    pub lateral: Handle<Face>,
}

/// A closed cylinder from `base` along `axis` for `height`, of `radius`: the plane through `base`
/// with normal `axis` stated as a datum ([`SketchPlane::from_origin_normal`]), a whole circle
/// about its origin, extruded. The seam sits where the product puts it — the frame's `+x̂`.
///
/// Panics on a statement the product refuses (a zero axis, a value outside the decimal window):
/// a fixture's input is the test's own.
pub fn cylinder(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    let plane = SketchPlane::from_origin_normal(base, axis).expect("a nonzero cylinder axis");
    extruded_circle(m, plane, radius, height)
}

/// [`cylinder`] with the seam written out — for a test whose proposition involves **where the
/// seam is**. The sketch frame is `+x̂ = seam`, `+ŷ = axis × seam`
/// ([`SketchPlane::from_axes`]), so the circle's seam lands on the `seam` side of the axis and
/// the frame's normal is `axis`. `seam` must be perpendicular to `axis`; `axis.any_perpendicular()`
/// is the usual choice. Exact when the two and their cross product are the decimals written (an
/// axis-aligned pair always is).
pub fn cylinder_with_seam(
    m: &mut Model,
    base: Point3,
    axis: Vector3,
    seam: Vector3,
    radius: f64,
    height: f64,
) -> CylinderSolid {
    extruded_circle(
        m,
        SketchPlane::from_axes(base, seam, axis.cross(seam)),
        radius,
        height,
    )
}

fn extruded_circle(m: &mut Model, plane: SketchPlane, radius: f64, height: f64) -> CylinderSolid {
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating the cylinder's base plane: {other:?}"),
    };
    let ring = Ring2d::circle(Point2::from_array([0.0, 0.0]), radius).expect("a stated radius");
    let profile = Profile2d::from_normalized_rings(ring, vec![]);
    match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: height,
        },
    ) {
        Ok(OpOutput::Extrude { solid, faces }) => {
            let [base, top, lateral] = faces[..] else {
                panic!("a circle extrudes to three faces, got {}", faces.len())
            };
            CylinderSolid {
                solid,
                base,
                top,
                lateral,
            }
        }
        other => panic!("extruding the cylinder's circle: {other:?}"),
    }
}

/// A boss with a window through its side: the cylinder `r = 1` on the `z` axis over `z ∈ [0, 4]`,
/// seam on the `seam_x`·`x̂` side (`±1`), less the box `x ∈ [0.5, 2]`, `|y| < half_width`,
/// `z ∈ [1.5, 2.5]`. The lateral face is absent from the line `x = 1, y = 0` over
/// `z ∈ [1.5, 2.5]`, and the window is an inner loop either way; `seam_x = 1` puts the rims' own
/// vertices (`θ = 0`) at the window's angle, `−1` opposite it.
///
/// The window's corners are rational exactly when `1 − half_width²` is a square (`0.6` → `x = 0.8`);
/// `0.3` gives `x = √0.91`.
pub fn windowed_boss(m: &mut Model, half_width: f64, seam_x: f64) -> Handle<Solid> {
    let boss = cylinder_with_seam(
        m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([seam_x, 0.0, 0.0]),
        1.0,
        4.0,
    )
    .solid;
    let window = cuboid(
        m,
        Point3::from_array([0.5, -half_width, 1.5]),
        Point3::from_array([2.0, half_width, 2.5]),
    );
    m.rebuild_adjacency();
    let out = boolean(m, BoolKind::Cut, boss, window).expect("the window");
    assert_eq!(out.len(), 1, "a window leaves one body");
    m.rebuild_adjacency();
    out[0]
}

/// **A cylinder whose two rims are both cut** — the lateral face with no seam edge. The cylinder
/// `r = 1` on the `z` axis over `z ∈ [0, 4]` (seam toward `+x̂`), combined by `kind` with the box
/// `x ∈ [0.5, 3]`, `y ∈ [y_min, 2]` over each end — `z ∈ [−1, 0.5]` and `z ∈ [3.5, 5]`. Each box
/// bites one rim, so the middle of the lateral wraps the axis between two **cut** rims; the region
/// walk emits it with no whole rim, and the face comes out bounded by those two rims alone — one
/// the outer loop, the other an inner one. `y_min = 0` stands the boxes' walls on the plane through
/// the axis, so each rim's `θ = 0` point — the chart's excluded point — is a box corner's pierce.
///
/// Volumes for `y_min = −2`: Fuse `4π + 30 − s`, Cut `4π − s`, with `s = acos 0.5 − 0.5·√0.75` the
/// area of the disk past `x = 0.5` — each bite is that segment `0.5` tall, and there are two.
///
/// Asserted here, because every test built on it measures this shape: one lateral face, no
/// self-adjacent edge on it, and exactly two loops.
pub fn cut_at_both_rims(m: &mut Model, kind: BoolKind, y_min: f64) -> Handle<Solid> {
    let c = cylinder(
        m,
        Point3::from_array([0.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        1.0,
        4.0,
    )
    .solid;
    let lo = cuboid(
        m,
        Point3::from_array([0.5, y_min, -1.0]),
        Point3::from_array([3.0, 2.0, 0.5]),
    );
    m.rebuild_adjacency();
    let a = boolean(m, kind, c, lo).expect("the lower bite");
    assert_eq!(a.len(), 1, "the lower bite leaves one body");
    m.rebuild_adjacency();
    let hi = cuboid(
        m,
        Point3::from_array([0.5, y_min, 3.5]),
        Point3::from_array([3.0, 2.0, 5.0]),
    );
    m.rebuild_adjacency();
    let out = boolean(m, kind, a[0], hi).expect("the upper bite");
    assert_eq!(out.len(), 1, "the upper bite leaves one body");
    m.rebuild_adjacency();
    let laterals: Vec<&Face> = m
        .shell(m.solid(out[0]).outer)
        .faces
        .iter()
        .map(|&f| m.face(f))
        .filter(|f| matches!(m.surface(f.surface), Surface::Cylinder { .. }))
        .collect();
    assert_eq!(laterals.len(), 1, "one lateral face");
    let face = laterals[0];
    assert_eq!(
        face.inner.len(),
        1,
        "two loops: one rim outer, the other inner"
    );
    assert!(
        std::iter::once(&face.outer)
            .chain(&face.inner)
            .flat_map(|l| &l.half_edges)
            .all(|he| {
                let [a, b] = m.edge(he.edge).surfaces;
                a != b
            }),
        "no seam edge"
    );
    out[0]
}

/// An axis-aligned box from `min` to `max` — [`cuboid_on`] with the floor at `min.z` and the
/// height `max.z − min.z`, **exactly**: the height is the difference of the two corners'
/// decimals, and it must be a decimal an `f64` carries (its shortest decimal), or the far cap
/// would not stand at `max.z`. Panics otherwise, saying so — state the floor and the height
/// instead ([`cuboid_on`]); a fixture's input is the test's own.
pub fn cuboid(m: &mut Model, min: Point3, max: Point3) -> Handle<Solid> {
    let [x0, y0, z0] = min.as_array();
    let [x1, y1, z1] = max.as_array();
    // ★ Reversed corners are refused rather than read: a lower `max.z` would become a negative
    // height — a box below the floor — and a reversed `x`/`y` a clockwise rectangle the extrude
    // quietly winds the other way.
    assert!(
        x1 > x0 && y1 > y0 && z1 > z0,
        "cuboid needs max > min on every axis: {min:?} .. {max:?}"
    );
    let dec = |v: f64| Rat::from_decimal(v).expect("a cuboid corner inside the decimal window");
    let diff = dec(z1)
        .checked_sub(dec(z0))
        .expect("a cuboid height inside `Rat`");
    let height = diff.to_f64();
    assert_eq!(
        Rat::from_decimal(height),
        Some(diff),
        "the height {z0} .. {z1} is no decimal an f64 carries — state the floor and the height \
         (`cuboid_on`)"
    );
    cuboid_on(m, [x0, y0], [x1, y1], z0, height)
}

/// An axis-aligned box over the rectangle `xy_min .. xy_max`, standing on the plane `z` and
/// reaching `height` along `+z` — or, for a negative `height`, hanging below the plane (the
/// extrude's negative distance). The plane is a datum through `(0, 0, z)`, so two boxes on one
/// plane in opposite directions **share** it — one statement of `z`, which a stack needs and two
/// heights summed in `f64` do not give.
///
/// The extrude's `[floor cap, far cap, walls…]` order is the fixture's: `faces[0]` is the cap on
/// the plane `z`, `faces[1]` the far one — for [`cuboid`] the bottom and the top. Asserted here,
/// with the box's six faces, twelve edges and eight corners, because tests read the caps by
/// position. The walls' order is the extrude's and nothing reads it.
pub fn cuboid_on(
    m: &mut Model,
    xy_min: [f64; 2],
    xy_max: [f64; 2],
    z: f64,
    height: f64,
) -> Handle<Solid> {
    let ([x0, y0], [x1, y1]) = (xy_min, xy_max);
    assert!(
        x1 > x0 && y1 > y0 && height != 0.0,
        "cuboid_on needs a rectangle and a nonzero height: {xy_min:?} .. {xy_max:?}, {height}"
    );
    let plane = SketchPlane::from_axes(
        Point3::from_array([0.0, 0.0, z]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    );
    let frame = match apply(
        m,
        &Operation::DatumPlane {
            def: DatumDef::Stated(plane),
        },
    ) {
        Ok(OpOutput::DatumPlane { frame, .. }) => frame,
        other => panic!("stating the box's floor plane: {other:?}"),
    };
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let profile = Profile2d::polygon(vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)])
        .expect("a rectangle inside the decimal window");
    let (solid, faces) = match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist: height,
        },
    ) {
        Ok(OpOutput::Extrude { solid, faces }) => (solid, faces),
        other => panic!("extruding the box's rectangle: {other:?}"),
    };
    // What the box is, and where its caps sit — the choices a positional reader leans on.
    let dec = |v: f64| Rat::from_decimal(v).expect("decimal");
    let far = dec(z)
        .checked_add(dec(height))
        .expect("a box height inside `Rat`")
        .to_f64();
    assert_eq!(faces.len(), 6, "a rectangle extrudes to six faces");
    let corners_of = |f: Handle<Face>| -> Vec<[u64; 3]> {
        m.face(f)
            .outer
            .half_edges
            .iter()
            .flat_map(|he| m.edge(he.edge).vertices)
            .map(|v| m.vertex_point(v).as_array().map(f64::to_bits))
            .collect()
    };
    for (f, at) in [(faces[0], z), (faces[1], far)] {
        assert!(
            corners_of(f).iter().all(|c| c[2] == at.to_bits()),
            "faces[0] is the cap on the plane, faces[1] the far one"
        );
    }
    let mut got: Vec<[u64; 3]> = faces.iter().flat_map(|&f| corners_of(f)).collect();
    got.sort_unstable();
    got.dedup();
    let mut want: Vec<[u64; 3]> = [x0, x1]
        .into_iter()
        .flat_map(|x| {
            [y0, y1]
                .into_iter()
                .flat_map(move |y| [z, far].map(|zz| [x, y, zz]))
        })
        .map(|c| c.map(f64::to_bits))
        .collect();
    want.sort_unstable();
    assert_eq!(
        got, want,
        "the box's eight corners are the stated ones, realized"
    );
    let edges: std::collections::HashSet<_> = faces
        .iter()
        .flat_map(|&f| m.face(f).outer.half_edges.iter().map(|he| he.edge))
        .collect();
    assert_eq!(edges.len(), 12, "a box has twelve edges");
    solid
}

/// What a pad or pocket built: every solid the boolean left, and the tool's far cap — its plane,
/// which way it faced, and its corners' mean — where the boss top or the pocket floor lies.
#[derive(Clone, Debug)]
pub struct Feature {
    pub solids: Vec<Handle<Solid>>,
    pub cap: Handle<Surface>,
    cap_normal: Vector3,
    cap_at: Point3,
}

/// Why [`pad`] or [`pocket`] built nothing.
#[derive(Debug, PartialEq)]
pub enum FeatureError {
    /// The extrude or the boolean refused.
    Op(OpError),
    /// The pad's fuse came back in pieces: two one-shell solids do that only when they never
    /// touched, so the footprint missed the face. The boolean answered right; the pad's premise
    /// broke — and a test that places footprints through [`face_plane`] reads this as «the frame
    /// it was told is not the one the tool was built in».
    Missed,
}

/// A boss on a planar `face`, the way an application builds one: the face's sketch frame
/// ([`face_sketch_frame`]), the profile extruded `dist` outward, and a `Fuse` with the face's solid
/// first.
pub fn pad(
    m: &mut Model,
    face: Handle<Face>,
    profile: Profile2d,
    dist: f64,
) -> Result<Feature, FeatureError> {
    feature(m, face, profile, dist, BoolKind::Fuse)
}

/// A pocket in a planar `face`: [`pad`]'s tool swept `dist` inward (a negative extrude in the same
/// frame) and cut away. Deeper than the body cuts through.
pub fn pocket(
    m: &mut Model,
    face: Handle<Face>,
    profile: Profile2d,
    dist: f64,
) -> Result<Feature, FeatureError> {
    feature(m, face, profile, -dist, BoolKind::Cut)
}

/// The `0.4` boss/pocket footprint on `[0.3, 0.7]²` of the unit cube's lid.
///
/// ★ **On a lid these are world coordinates.** The sketch origin is the world origin projected
/// onto the face's plane, and the axes are the arbitrary-axis convention's, which for `n = ẑ` are
/// `u = +x̂`, `v = +ŷ` — so a frame point `(a, b)` is world `(a, b, 1)`, the identity.
pub fn lid_square() -> Profile2d {
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    Profile2d::polygon(vec![p(0.3, 0.7), p(0.3, 0.3), p(0.7, 0.3), p(0.7, 0.7)])
        .expect("a square inside the decimal window")
}

/// The unit cube with a [`lid_square`] pocket `0.5` deep in its top face: the void is
/// `[0.3, 0.7]² × [0.5, 1]` and the solid measures `1 − 0.16·0.5 = 0.92`.
///
/// ★ Its lid carries an inner loop, and that is what its consumers measure — so it is asserted
/// here. A footprint the frame places elsewhere (on the lid's corner, say) builds a solid every
/// test still accepts, and a boolean scored only against another kernel stays green on it.
pub fn pocketed_cube() -> (Model, Handle<Solid>) {
    let mut m = Model::new();
    let p = |x: f64, y: f64| Point2::from_array([x, y]);
    let square = Profile2d::polygon(vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)])
        .expect("the unit square");
    let frame = SketchFrame::world(&m, Axis::Z);
    let top = match apply(
        &mut m,
        &Operation::Extrude {
            frame,
            profile: square,
            dist: 1.0,
        },
    ) {
        Ok(OpOutput::Extrude { faces, .. }) => faces[1],
        other => panic!("extruding the unit cube: {other:?}"),
    };
    let solid = pocket(&mut m, top, lid_square(), 0.5)
        .expect("the pocket")
        .solid();
    let holed: Vec<usize> = m
        .shell(m.solid(solid).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).inner.len())
        .filter(|&n| n > 0)
        .collect();
    assert_eq!(
        holed,
        [1],
        "the pocketed cube's lid carries the pocket as its one hole"
    );
    (m, solid)
}

impl From<OpError> for FeatureError {
    fn from(e: OpError) -> Self {
        Self::Op(e)
    }
}

fn feature(
    m: &mut Model,
    face: Handle<Face>,
    profile: Profile2d,
    dist: f64,
    kind: BoolKind,
) -> Result<Feature, FeatureError> {
    let frame = face_sketch_frame(m, face).map_err(FeatureError::Op)?;
    let body = *m
        .live_solids()
        .iter()
        .find(|&&s| m.shell(m.solid(s).outer).faces.contains(&face))
        .expect("face_sketch_frame found the face on a live solid's outer shell");
    let (tool, faces) = match apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .map_err(FeatureError::Op)?
    {
        OpOutput::Extrude { solid, faces } => (solid, faces),
        other => unreachable!("an extrude answered {other:?}"),
    };
    let cap = m.face(faces[1]).surface;
    let (cap_normal, cap_at) = facing_and_mean(m, faces[1]);
    let solids = match apply(
        m,
        &Operation::Boolean {
            kind,
            a: body,
            b: tool,
        },
    )
    .map_err(FeatureError::Op)?
    {
        OpOutput::Boolean { solids } => solids,
        other => unreachable!("a boolean answered {other:?}"),
    };
    if matches!(kind, BoolKind::Fuse) && solids.len() > 1 {
        return Err(FeatureError::Missed);
    }
    Ok(Feature {
        solids,
        cap,
        cap_normal,
        cap_at,
    })
}

/// A planar face's outward unit normal and its outer corners' mean, from the caches — a test's
/// way of pointing at a face, not a judgment.
fn facing_and_mean(m: &Model, f: Handle<Face>) -> (Vector3, Point3) {
    let face = m.face(f);
    let nacre_geom::Surface::Plane(pl) = m.surface_cache(face.surface) else {
        panic!("a feature's cap is planar")
    };
    let pts: Vec<[f64; 3]> = face
        .outer
        .half_edges
        .iter()
        .map(|&he| m.vertex_point(m.he_start(he)).as_array())
        .collect();
    let n = pts.len() as f64;
    let mean = std::array::from_fn(|k| pts.iter().map(|p| p[k]).sum::<f64>() / n);
    (
        pl.normal() * f64::from(face.orientation.sign()),
        Point3::from_array(mean),
    )
}

impl Feature {
    /// The one solid the feature left. Panics when it left none or several — a test that expects
    /// a severing cut reads [`Feature::solids`].
    pub fn solid(&self) -> Handle<Solid> {
        match self.solids.as_slice() {
            [s] => *s,
            other => panic!("the feature left {} solids, not one", other.len()),
        }
    }

    /// The result face where the tool's far cap went — the boss top or the pocket floor — picked
    /// the way a test picks a face: on the cap's plane, facing the cap's way, nearest the cap's
    /// corners' mean. **A selection, not a judgment.** "On the cap's plane" is the cap's handle
    /// *or* the cap's plane read off the caches, because a result face carries its plane class's
    /// representative: where the body already held that plane under another handle (a frame-node
    /// statement on a tilted face), the cap survives on that one. Two caps on one plane (two
    /// bosses of one height) are told apart by where they stand. Panics when nothing qualifies.
    pub fn cap_face(&self, m: &Model) -> Handle<Face> {
        let planar = |f: Handle<Face>| {
            matches!(
                m.surface_cache(m.face(f).surface),
                nacre_geom::Surface::Plane(_)
            )
        };
        let near = |f: Handle<Face>| -> Option<f64> {
            let (n, at) = facing_and_mean(m, f);
            let off = (self.cap_at - at).dot(self.cap_normal).abs();
            let on = m.face(f).surface == self.cap
                || (n.dot(self.cap_normal) > 1.0 - 1e-9 && off < 1e-9);
            on.then(|| (at - self.cap_at).norm())
        };
        self.solids
            .iter()
            .flat_map(|&s| {
                let solid = m.solid(s);
                std::iter::once(solid.outer)
                    .chain(solid.cavities.iter().copied())
                    .flat_map(|sh| m.shell(sh).faces.iter().copied())
                    .collect::<Vec<_>>()
            })
            .filter(|&f| planar(f))
            .filter_map(|f| near(f).map(|d| (f, d)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(f, _)| f)
            .expect("a result face where the tool's far cap went")
    }
}
