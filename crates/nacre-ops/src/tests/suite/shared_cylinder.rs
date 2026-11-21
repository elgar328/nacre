//! **One cylinder surface under two solids** — the scoreboard of booleans whose two operands both
//! carry lateral faces on one cylinder: under one handle, and under two.
//!
//! Each family is built in one frame from one circle, so the two operands' laterals intern to one
//! `Handle<Surface>` (asserted per family before anything else): two identical prisms, a circle
//! swept both ways, a plate bored by the circle and the circle's pin, a short and a tall prism, a
//! half disk beside the whole one (a ruling station on the shared class), a notched plate over the
//! half pin below the notch (the material side flips across the station), a bored block above the
//! bore's pin (it flips across a rim), the two halves of a disk, a filleted corner swept
//! both ways and the fillet's pin (tangent walls on the shared class), two plates rounded on
//! opposite quarters of one circle, two identical cylinders on a stated datum, two chord segments
//! of the disk meeting on one ruling (each lateral ending there from its own side), two half disks
//! whose statements of the diameter's rulings overlap in part and one inside the other with both
//! ends inside (the rulings cut at each other's ends), a bored block below the bore's pin (the
//! rim joined the other way round), two blocks bitten by one circle into cusps that meet on its
//! tangent ruling — alone, mirrored and offset along it, joined round the bore, joined on their
//! floor, and beside a plane through the axis — a cusp and the block bitten beyond its floor (the
//! union the cylinder's outside whole round that ruling), and two blocks bitten into the walls of a
//! wedge that is empty below the ruling, apart and joined above the bore. Every family is run as
//! Fuse, both Cuts and Common.
//!
//! The class carries both operands' lateral faces, and the chart reads its cells for both — the
//! inside chamber off the arrangement's labels, and the outside as the inside with every solid
//! whose lateral face is there flipped (`cyl_chart::read_cell`).
//!
//! The same cylinder stated under **two** handles is one class too
//! (`planes::plane_classes`): every family is run again with its two operands' cylinders under two
//! different statements ([`statement_pairs`]) — the seam turned, the origin moved along the axis,
//! the axis reversed ([`Says`]) — and answers as it does under one handle, as do twins a motion or
//! a frame restates, one of them with the axis written at another length. A motion that folds back
//! to no motion at all is one statement.

use super::scoreboard::{OPS, Want, check, counted, refused, run};
use super::*;
use crate::{Edge2d, Ring2d};
use nacre_exact::Rat;

type Build = fn(&mut Model, Says) -> (Handle<Solid>, Handle<Solid>);

/// A twin pair: the first operand as written, the second stating the first's cylinder again.
type Twins = fn(&mut Model) -> (Handle<Solid>, Handle<Solid>);

fn p(x: f64, y: f64) -> nacre_math::Point2 {
    nacre_math::Point2::from_array([x, y])
}

fn profile(edges: Vec<Stated>) -> Profile2d {
    stated(edges).unwrap().remove(0)
}

/// **How a family's solids state their cylinders** — one solid, four statements of it. Every sweep
/// of a family writes the same points of the world; what differs is the sketch frame the circle is
/// written in, and so the cylinder's statement: its seam direction, the place of its origin on the
/// axis, or the way its axis points.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Says {
    /// The world frame on `z = 0`: axis `+ẑ`, seam `+x̂`, origin on `z = 0`.
    Plain,
    /// The frame `+u = ŷ` on `z = 0`, the profile written a quarter turn back: seam `+ŷ`.
    Seam,
    /// A datum at the far end, swept back: the origin at the other end of the sweep.
    Origin,
    /// The plane `z = 0` facing `−ẑ` with `+u = x̂`, the profile written mirrored in `y`: axis
    /// `−ẑ`, seam and origin as [`Says::Plain`].
    Reversed,
}

/// The four statements.
const STATEMENTS: [Says; 4] = [Says::Plain, Says::Seam, Says::Origin, Says::Reversed];

/// **Every ordered pair of two different statements** — the first operand's, then the second's.
/// Ordered, because a class reads in its lowest row's statement, the first operand's: `(x, y)`
/// restates `y` in `x` and `(y, x)` restates `x` in `y`, while `Op::CutBa` swaps the two solids of
/// one pair, not the statements they are built in.
fn statement_pairs() -> impl Iterator<Item = (Says, Says)> {
    STATEMENTS.into_iter().flat_map(|x| {
        STATEMENTS
            .into_iter()
            .filter(move |&y| y != x)
            .map(move |y| (x, y))
    })
}

/// `profile` under the map `at` of the frame's coordinates — a turn (`mirrored: false`) or a
/// mirror. A whole circle keeps its seam on its frame's `+u` side: its one vertex moves with the
/// centre and is not turned about it, so the circle's world seam turns with the frame. A mirror
/// reverses a ring's sense, so the ring is written the other way round.
fn mapped(profile: &Profile2d, at: impl Fn([Rat; 2]) -> [Rat; 2], mirrored: bool) -> Profile2d {
    let ring = |r: &Ring2d| {
        let n = r.len();
        let arc = |e: &Edge2d| match *e {
            Edge2d::Line => Edge2d::Line,
            Edge2d::Arc { center, r2, ccw } => Edge2d::Arc {
                center: at(center),
                r2,
                ccw,
            },
        };
        if n == 1 {
            let Edge2d::Arc { center, .. } = r.edges()[0] else {
                unreachable!("a one-vertex ring is a circle")
            };
            let (v, c) = (r.vertices()[0], at(center));
            let seam = [
                c[0].checked_add(v[0].checked_sub(center[0]).unwrap())
                    .unwrap(),
                c[1].checked_add(v[1].checked_sub(center[1]).unwrap())
                    .unwrap(),
            ];
            return Ring2d::normalized(vec![seam], vec![arc(&r.edges()[0])]);
        }
        let vs: Vec<[Rat; 2]> = r.vertices().iter().map(|&v| at(v)).collect();
        if !mirrored {
            return Ring2d::normalized(vs, r.edges().iter().map(arc).collect());
        }
        // The ring the other way round from its first vertex: step `j` is the mirrored step
        // `n − j − 1` walked backwards, which turns the way the original did.
        let order = |j: usize| (n - j) % n;
        Ring2d::normalized(
            (0..n).map(|j| vs[order(j)]).collect(),
            (0..n).map(|j| arc(&r.edges()[n - j - 1])).collect(),
        )
    };
    Profile2d::from_normalized_rings(
        ring(profile.outer()),
        profile.holes().iter().map(ring).collect(),
    )
}

/// `profile` written on the world's `z = 0` and swept `dist` along `+ẑ`, as `says` states it.
fn sweep(m: &mut Model, says: Says, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let neg = |x: Rat| Rat::from_int(0).checked_sub(x).unwrap();
    let (frame, profile, dist) = match says {
        Says::Plain => (SketchFrame::world(m, Axis::Z), profile, dist),
        Says::Seam => {
            let plane = m.seeded_plane(Axis::Z);
            let frame = SketchFrame::named(
                m,
                plane,
                Point3::origin(),
                Vector3::from_array([0.0, 1.0, 0.0]),
            )
            .expect("the frame");
            // `+u = ŷ`, `+v = −x̂`: the world's `(x, y)` is the frame's `(y, −x)`.
            (frame, mapped(&profile, |q| [q[1], neg(q[0])], false), dist)
        }
        Says::Origin => (datum(m, [0.0, 0.0, dist], [0.0, 0.0, 1.0]), profile, -dist),
        Says::Reversed => {
            let frame = frame_on(m, [0.0, 0.0, 0.0], ([1.0, 0.0, 0.0], [0.0, -1.0, 0.0]));
            // `+u = x̂`, `+v = −ŷ`: the world's `(x, y)` is the frame's `(x, −y)`.
            (frame, mapped(&profile, |q| [q[0], neg(q[1])], true), -dist)
        }
    };
    extrude_in(m, frame, profile, dist)
}

/// The unit cylinder from `base` up `height` along `+ẑ` on a stated datum, as `says` states it.
fn upright(m: &mut Model, says: Says, base: [f64; 3], height: f64) -> Handle<Solid> {
    let (z, x, y) = (
        Vector3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([1.0, 0.0, 0.0]),
        Vector3::from_array([0.0, 1.0, 0.0]),
    );
    let at = Point3::from_array(base);
    let top = Point3::from_array([base[0], base[1], base[2] + height]);
    match says {
        Says::Plain => crate::fixtures::cylinder(m, at, z, 1.0, height),
        Says::Seam => crate::fixtures::cylinder_with_seam(m, at, z, y, 1.0, height),
        Says::Origin => crate::fixtures::cylinder(m, top, z, 1.0, -height),
        Says::Reversed => crate::fixtures::cylinder_with_seam(m, at, z * -1.0, x, 1.0, -height),
    }
    .solid
}

/// The unit circle about the origin.
fn disk() -> Profile2d {
    circle_profile([0.0, 0.0], 1.0)
}

/// The unit circle's half on `y ≥ 0`, closed by its diameter.
fn upper_half() -> Profile2d {
    profile(vec![
        line(p(-1.0, 0.0), p(1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(1.0, 0.0), 2),
    ])
}

/// The unit circle's half on `y ≤ 0`, closed by its diameter.
fn lower_half() -> Profile2d {
    profile(vec![
        line(p(1.0, 0.0), p(-1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(-1.0, 0.0), 2),
    ])
}

/// `[0, 4] × [0, 2]` with its `(4, 0)` corner rounded by the unit circle about `(3, 1)`.
fn filleted_right() -> Profile2d {
    profile(vec![
        line(p(0.0, 0.0), p(3.0, 0.0)),
        arc_turns(p(3.0, 1.0), p(3.0, 0.0), 1),
        line(p(4.0, 1.0), p(4.0, 2.0)),
        line(p(4.0, 2.0), p(0.0, 2.0)),
        line(p(0.0, 2.0), p(0.0, 0.0)),
    ])
}

/// `[2, 6] × [0, 2]` with its `(2, 0)` corner rounded by the same circle — the opposite quarter.
fn filleted_left() -> Profile2d {
    profile(vec![
        line(p(3.0, 0.0), p(6.0, 0.0)),
        line(p(6.0, 0.0), p(6.0, 2.0)),
        line(p(6.0, 2.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(2.0, 1.0)),
        arc_turns(p(3.0, 1.0), p(2.0, 1.0), 1),
    ])
}

fn identical(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, s, disk(), 2.0), sweep(m, s, disk(), 2.0))
}

fn swept_both_ways(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, s, disk(), 2.0), sweep(m, s, disk(), -2.0))
}

/// A `4 × 4 × 2` plate bored by the unit circle, and the circle's pin of the plate's height.
fn bore_and_pin(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let square =
        Profile2d::polygon(vec![p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0)]).unwrap();
    let plate = sweep(m, s, square, 2.0);
    let tool = sweep(m, s, disk(), 2.0);
    let bored = crate::boolean(m, BoolKind::Cut, plate, tool).expect("the bore")[0];
    (bored, sweep(m, s, disk(), 2.0))
}

fn short_and_tall(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, s, disk(), 2.0), sweep(m, s, disk(), 3.0))
}

fn half_and_whole(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, s, upper_half(), 2.0), sweep(m, s, disk(), 2.0))
}

fn short_half_and_tall_whole(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (sweep(m, s, upper_half(), 2.0), sweep(m, s, disk(), 3.0))
}

/// `[−2, 2] × [0, 2]` notched by the upper half disk, and the lower half disk's pin: the shared
/// lateral is the notch's on `y > 0` (material outside) and the pin's on `y < 0` (material inside).
fn notch_and_half_pin(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let notched = profile(vec![
        line(p(-2.0, 0.0), p(-1.0, 0.0)),
        arc_turns(p(0.0, 0.0), p(-1.0, 0.0), -2),
        line(p(1.0, 0.0), p(2.0, 0.0)),
        line(p(2.0, 0.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(-2.0, 2.0)),
        line(p(-2.0, 2.0), p(-2.0, 0.0)),
    ]);
    (sweep(m, s, notched, 2.0), sweep(m, s, lower_half(), 2.0))
}

/// A `4 × 4 × 2` block on `z ∈ [2, 4]` bored through by the unit circle swept from `z = 0`, and
/// the circle's pin on `z ∈ [0, 2]` below it: the shared lateral keeps the material outside above
/// the rim at `z = 2` and inside below it, and four faces meet on that circle.
fn bore_above_pin(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, 2.0], [2.0, 2.0, 4.0]);
    let tool = sweep(m, s, disk(), 4.0);
    let bored = crate::boolean(m, BoolKind::Cut, block, tool).expect("the bore")[0];
    (bored, sweep(m, s, disk(), 2.0))
}

/// The unit disk's two halves: one lateral across the two stations of the diameter, and the two
/// flats on it facing each other.
fn two_halves(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, s, upper_half(), 2.0),
        sweep(m, s, lower_half(), 2.0),
    )
}

fn fillet_swept_both_ways(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, s, filleted_right(), 2.0),
        sweep(m, s, filleted_right(), -2.0),
    )
}

/// The filleted plate and a pin on its fillet's circle, taller than the plate: the plate's two
/// walls are tangent to the pin's lateral.
fn fillet_and_its_pin(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, s, filleted_right(), 2.0),
        sweep(m, s, circle_profile([3.0, 1.0], 1.0), 3.0),
    )
}

fn opposite_quarters(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        sweep(m, s, filleted_right(), 2.0),
        sweep(m, s, filleted_left(), 2.0),
    )
}

/// The unit disk's segment beyond the chord `x + y = 1` and its segment beyond `x − y = 1`: the
/// two chords meet on the ruling `(1, 0)`, where one lateral face ends from `+θ` and the other's
/// from `−θ`, each against its own chord.
fn two_chord_segments(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let upper = profile(vec![
        arc_turns(p(0.0, 0.0), p(1.0, 0.0), 1),
        line(p(0.0, 1.0), p(1.0, 0.0)),
    ]);
    let lower = profile(vec![
        arc_turns(p(0.0, 0.0), p(0.0, -1.0), 1),
        line(p(1.0, 0.0), p(0.0, -1.0)),
    ]);
    (sweep(m, s, upper, 2.0), sweep(m, s, lower, 2.0))
}

/// Two identical cylinders on a stated datum plane — the other frame road.
fn identical_on_a_datum(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let on_datum = |m: &mut Model| upright(m, s, [0.0, 0.0, 0.0], 2.0);
    (on_datum(m), on_datum(m))
}

/// `solid` with everything outside `z ∈ [lo, hi]` cut away by boxes wider than the unit disk — the
/// cut keeps the lateral's cylinder handle.
fn trimmed(m: &mut Model, mut solid: Handle<Solid>, lo: f64, hi: f64) -> Handle<Solid> {
    let below = crate::fixtures::cuboid(m, [-2.0, -2.0, lo - 2.0], [2.0, 2.0, lo]);
    solid = crate::boolean(m, BoolKind::Cut, solid, below).expect("the trim below")[0];
    let above = crate::fixtures::cuboid(m, [-2.0, -2.0, hi], [2.0, 2.0, hi + 2.0]);
    crate::boolean(m, BoolKind::Cut, solid, above).expect("the trim above")[0]
}

/// The upper half disk on `z ∈ [0, 2]` and on `z ∈ [1, 3]`: the two statements of each ruling of
/// the diameter overlap on `[1, 2]`.
fn halves_overlapping(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let low = sweep(m, s, upper_half(), 2.0);
    let high = sweep(m, s, upper_half(), 3.0);
    (low, trimmed(m, high, 1.0, 3.0))
}

/// The upper half disk on `z ∈ [0, 3]` and on `z ∈ [1, 2]`: one statement of each ruling lies
/// inside the other, both its ends inside.
fn half_around_half(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let tall = sweep(m, s, upper_half(), 3.0);
    let inner = sweep(m, s, upper_half(), 3.0);
    (tall, trimmed(m, inner, 1.0, 2.0))
}

/// A `4 × 4 × 2` block on `z ∈ [−4, −2]` bored through by the unit circle swept down from `z = 0`,
/// and the circle's pin on `z ∈ [−2, 0]` above it — [`bore_above_pin`] the other way up, so each
/// face meets the rim at `z = −2` from the other side.
fn bore_below_pin(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let block = crate::fixtures::cuboid(m, [-2.0, -2.0, -4.0], [2.0, 2.0, -2.0]);
    let tool = sweep(m, s, disk(), -4.0);
    let bored = crate::boolean(m, BoolKind::Cut, block, tool).expect("the bore")[0];
    (bored, sweep(m, s, disk(), -2.0))
}

/// The block from `min` to `max` less the unit cylinder about the axis line `x = 1, y = 1` on
/// `z ∈ [lo, lo + len]` — one statement wherever it is written, so the bites of one family share
/// its handle.
fn bitten(
    m: &mut Model,
    s: Says,
    min: [f64; 3],
    max: [f64; 3],
    lo: f64,
    len: f64,
) -> Handle<Solid> {
    bitten_about(m, s, 1.0, min, max, lo, len)
}

/// [`bitten`] about the axis line `x = 1, y = axis_y`.
fn bitten_about(
    m: &mut Model,
    s: Says,
    axis_y: f64,
    min: [f64; 3],
    max: [f64; 3],
    lo: f64,
    len: f64,
) -> Handle<Solid> {
    let block = crate::fixtures::cuboid(m, min, max);
    let tool = upright(m, s, [1.0, axis_y, lo], len);
    crate::boolean(m, BoolKind::Cut, block, tool).expect("the bite")[0]
}

/// The unit cube and the next one along `x`, each less the circle about `(1, 1)`: two cusps whose
/// floors `y = 0` and laterals end on the tangent ruling `(1, 0)` from either side, meeting there
/// alone.
fn two_cusps(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten(m, s, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], -1.0, 3.0),
        bitten(m, s, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0], -1.0, 3.0),
    )
}

/// [`two_cusps`] mirrored across `y = 0` — the floor's axis side the other way — and two high, the
/// second raised by one: the cusps meet on the tangent ruling `(1, 0)` over `z ∈ [1, 2]` only, and
/// each runs on along the line alone beyond it.
fn two_cusps_under_the_axis_offset_along_it(
    m: &mut Model,
    s: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten_about(m, s, -1.0, [0.0, -1.0, 0.0], [1.0, 0.0, 2.0], -1.0, 5.0),
        bitten_about(m, s, -1.0, [1.0, -1.0, 1.0], [2.0, 0.0, 3.0], -1.0, 5.0),
    )
}

/// The blocks `[−1, 1] × [0, 3]` and `[1, 3] × [0, 3]`, each less the circle about `(1, 1)`: the
/// cusps meet on the tangent ruling `(1, 0)` and the blocks share the face `x = 1` above the
/// circle — one body round a bore tangent to its floor, whose boundary touches itself along the
/// line.
fn two_cusps_joined_round_the_bore(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten(m, s, [-1.0, 0.0, 0.0], [1.0, 3.0, 1.0], -1.0, 3.0),
        bitten(m, s, [1.0, 0.0, 0.0], [3.0, 3.0, 1.0], -1.0, 3.0),
    )
}

/// The block `[−1, 1] × [0, 3]` less the circle, and the block `[−1, 3] × [0, 3] × [0, 2]` with
/// `[−1, 1] × [0, 1]` notched out under its top and the circle bitten on `z ∈ [−1, 1]`: the cusps
/// meet on the tangent ruling `(1, 0)` below `z = 1`, and the two floors `y = 0` meet again along
/// `z = 1` beside it — one body round a bore tangent to its floor, the floor one face of its plane.
fn two_cusps_joined_on_their_floor(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    let first = bitten(m, s, [-1.0, 0.0, 0.0], [1.0, 3.0, 1.0], -1.0, 3.0);
    let block = crate::fixtures::cuboid(m, [-1.0, 0.0, 0.0], [3.0, 3.0, 2.0]);
    let notch = crate::fixtures::cuboid(m, [-2.0, -1.0, -1.0], [1.0, 4.0, 1.0]);
    let stepped = crate::boolean(m, BoolKind::Cut, block, notch).expect("the notch")[0];
    let tool = upright(m, s, [1.0, 1.0, -1.0], 2.0);
    let second = crate::boolean(m, BoolKind::Cut, stepped, tool).expect("the bite")[0];
    (first, second)
}

/// The first cusp of [`two_cusps`] and the block `[1, 2] × [−1, 1]` less the same circle: the
/// second's face `x = 1` runs through the axis and ends on the tangent ruling `(1, 0)` from below,
/// where its lateral ends too.
fn a_cusp_beside_a_plane_through_the_axis(
    m: &mut Model,
    s: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten(m, s, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], -1.0, 3.0),
        bitten(m, s, [1.0, -1.0, 0.0], [2.0, 1.0, 1.0], -1.0, 3.0),
    )
}

/// The prism over `pts`, one tall, less the circle about `(1, 1)` — the bite [`bitten`] takes, so a
/// family's prisms and blocks share its handle.
fn bitten_prism(m: &mut Model, s: Says, pts: &[[f64; 2]]) -> Handle<Solid> {
    let outline =
        Profile2d::polygon(pts.iter().map(|q| p(q[0], q[1])).collect()).expect("a polygon");
    let prism = sweep(m, s, outline, 1.0);
    let tool = upright(m, s, [1.0, 1.0, -1.0], 3.0);
    crate::boolean(m, BoolKind::Cut, prism, tool).expect("the bite")[0]
}

/// The second cusp of [`two_cusps`], its floor `y = 0` ending on the tangent ruling `(1, 0)`, and
/// the L `[−1, 3] × [−1, 0] ∪ [−1, 1] × [0, 1]` less the same circle: the first cusp and everything
/// beyond the floor's plane, its face `y = 0` on `x ≥ 1` ending on the line where its lateral does.
/// Around the line the union is the cylinder's outside whole.
fn a_cusp_and_the_block_bitten_beyond_its_floor(
    m: &mut Model,
    s: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten(m, s, [1.0, 0.0, 0.0], [2.0, 1.0, 1.0], -1.0, 3.0),
        bitten_prism(
            m,
            s,
            &[
                [-1.0, -1.0],
                [3.0, -1.0],
                [3.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [-1.0, 1.0],
            ],
        ),
    )
}

/// The triangles `x + y ≥ 1` and `x − y ≤ 1` of `[−1, 3] × [−2, 1]` less the circle: each lateral
/// ends on the tangent ruling `(1, 0)` on its own wall, which runs on below the line, and between
/// the walls there the right-angled wedge `y ≤ −|x − 1|` is empty. Above the line the walls lie in
/// the disk up to `(0, 1)` and `(2, 1)`, where they meet the circle again on the triangles' top
/// corners.
fn two_blocks_bitten_into_a_wedges_walls(m: &mut Model, s: Says) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten_prism(m, s, &[[0.0, 1.0], [3.0, -2.0], [3.0, 1.0]]),
        bitten_prism(m, s, &[[-1.0, 1.0], [-1.0, -2.0], [2.0, 1.0]]),
    )
}

/// [`two_blocks_bitten_into_a_wedges_walls`] up to `y = 3`: above the disk the two blocks overlap,
/// and their union is one body round the bore touching itself along the line.
fn two_blocks_bitten_into_a_wedges_walls_joined_above_the_bore(
    m: &mut Model,
    s: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bitten_prism(m, s, &[[3.0, -2.0], [3.0, 3.0], [-1.0, 3.0], [-1.0, 2.0]]),
        bitten_prism(m, s, &[[-1.0, -2.0], [3.0, 2.0], [3.0, 3.0], [-1.0, 3.0]]),
    )
}

const FAMILIES: [(&str, Build); 25] = [
    ("two identical prisms", identical),
    ("one circle swept both ways", swept_both_ways),
    ("a bored plate and the bore's pin", bore_and_pin),
    ("a short and a tall prism", short_and_tall),
    ("a half disk and the whole", half_and_whole),
    (
        "a short half disk and a tall whole",
        short_half_and_tall_whole,
    ),
    ("a notched plate and the half pin", notch_and_half_pin),
    ("a bored block above the bore's pin", bore_above_pin),
    ("the two halves of a disk", two_halves),
    ("a fillet swept both ways", fillet_swept_both_ways),
    ("a filleted plate and its fillet's pin", fillet_and_its_pin),
    ("two plates rounded on opposite quarters", opposite_quarters),
    ("two identical cylinders on a datum", identical_on_a_datum),
    ("two chord segments meeting on a ruling", two_chord_segments),
    ("two half disks on [0, 2] and [1, 3]", halves_overlapping),
    (
        "a half disk on [0, 3] around one on [1, 2]",
        half_around_half,
    ),
    ("a bored block below the bore's pin", bore_below_pin),
    ("two cusps meeting on their tangent ruling", two_cusps),
    (
        "two cusps under the axis, offset along it",
        two_cusps_under_the_axis_offset_along_it,
    ),
    (
        "two cusped blocks joined round the bore",
        two_cusps_joined_round_the_bore,
    ),
    (
        "two cusped blocks joined on their floor",
        two_cusps_joined_on_their_floor,
    ),
    (
        "a cusp beside a plane through the axis",
        a_cusp_beside_a_plane_through_the_axis,
    ),
    (
        "a cusp and the block bitten beyond its floor",
        a_cusp_and_the_block_bitten_beyond_its_floor,
    ),
    (
        "two blocks bitten into a wedge's walls",
        two_blocks_bitten_into_a_wedges_walls,
    ),
    (
        "two blocks bitten into a wedge's walls, joined above the bore",
        two_blocks_bitten_into_a_wedges_walls_joined_above_the_bore,
    ),
];

/// The cylinder surfaces of a solid's lateral faces.
fn lateral_surfaces(m: &Model, s: Handle<Solid>) -> Vec<Handle<Surface>> {
    let mut out: Vec<Handle<Surface>> = m
        .shell(m.solid(s).outer)
        .faces
        .iter()
        .map(|&f| m.face(f).surface)
        .filter(|&h| matches!(m.surface(h), Surface::Cylinder { .. }))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// **Every family shares one cylinder handle between its operands** — the population this
/// scoreboard is about, checked before its answers are read.
#[test]
fn each_familys_operands_share_one_cylinder_handle() {
    for (name, build) in FAMILIES {
        let mut m = Model::new();
        let (a, b) = build(&mut m, Says::Plain);
        let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
        assert!(
            la.iter().any(|h| lb.contains(h)),
            "{name}: the operands' laterals are {la:?} and {lb:?}"
        );
    }
}

const PI: f64 = std::f64::consts::PI;

/// A unit disk's segment cut off by a chord a quarter turn long, swept 2: `(π/4 − 1/2)·2`.
const SEGMENT: f64 = PI / 2.0 - 1.0;

/// Each family's answers as Fuse, the first less the second, the second less the first, Common —
/// the volumes set out from the shapes (a fillet's quarter is `1 − π/4` short of its square's
/// corner).
fn wants(name: &str) -> [Want; 4] {
    let plate = 14.0 + PI / 2.0; // `filleted_right`, swept 2
    let none = || counted(&[], &[]);
    match name {
        "two identical prisms" | "two identical cylinders on a datum" => [
            counted(&[2.0 * PI], &[3]),
            none(),
            none(),
            counted(&[2.0 * PI], &[3]),
        ],
        "one circle swept both ways" => [
            counted(&[4.0 * PI], &[3]),
            counted(&[2.0 * PI], &[3]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        "a bored plate and the bore's pin" => [
            counted(&[32.0], &[6]),
            counted(&[32.0 - 2.0 * PI], &[7]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        "a short and a tall prism" => [
            counted(&[3.0 * PI], &[3]),
            none(),
            counted(&[PI], &[3]),
            counted(&[2.0 * PI], &[3]),
        ],
        "a half disk and the whole" => [
            counted(&[2.0 * PI], &[3]),
            none(),
            counted(&[PI], &[4]),
            counted(&[PI], &[4]),
        ],
        // The tall less the short is the lower half under a whole slab: one lateral face round
        // both, the flat, the overhang's half disk, the two caps.
        "a short half disk and a tall whole" => [
            counted(&[3.0 * PI], &[3]),
            none(),
            counted(&[2.0 * PI], &[5]),
            counted(&[PI], &[4]),
        ],
        // The notch and the pin meet on the two rulings alone: two bodies — the pin's four faces,
        // and the notched plate's eight (its bottom two pieces either side of the notch).
        "a notched plate and the half pin" => [
            counted(&[PI, 16.0 - PI], &[4, 8]),
            counted(&[16.0 - PI], &[8]),
            counted(&[PI], &[4]),
            none(),
        ],
        // The pin meets the block on the rim circle alone: two bodies — the pin's three faces, the
        // bored block's seven.
        "a bored block above the bore's pin" | "a bored block below the bore's pin" => [
            counted(&[2.0 * PI, 32.0 - 2.0 * PI], &[3, 7]),
            counted(&[32.0 - 2.0 * PI], &[7]),
            counted(&[2.0 * PI], &[3]),
            none(),
        ],
        // The two segments meet on the ruling `(1, 0)` alone: two bodies of four faces each.
        "two chord segments meeting on a ruling" => [
            counted(&[SEGMENT, SEGMENT], &[4, 4]),
            counted(&[SEGMENT], &[4]),
            counted(&[SEGMENT], &[4]),
            none(),
        ],
        "the two halves of a disk" => [
            counted(&[2.0 * PI], &[3]),
            counted(&[PI], &[4]),
            counted(&[PI], &[4]),
            none(),
        ],
        "a fillet swept both ways" => [
            counted(&[2.0 * plate], &[7]),
            counted(&[plate], &[7]),
            counted(&[plate], &[7]),
            none(),
        ],
        // The pin's circle is also tangent to the plate's top wall, at `(3, 2)`: the plate less the
        // pin is two lumps meeting only along that line — the corner of five faces and the rest of
        // six, eleven as OCCT builds them. In the union the plate's top less the pin's disk is two
        // faces meeting at that point, the corner beyond the disk and the rest.
        "a filleted plate and its fillet's pin" => [
            counted(&[plate + PI], &[9]),
            counted(&[2.0 - PI / 2.0, 12.0 - PI], &[5, 6]),
            counted(&[PI], &[3]),
            counted(&[2.0 * PI], &[3]),
        ],
        "two plates rounded on opposite quarters" => [
            counted(&[24.0], &[6]),
            counted(&[10.0 - PI / 2.0], &[7]),
            counted(&[10.0 - PI / 2.0], &[7]),
            counted(&[4.0 + PI], &[6]),
        ],
        // Each half disk is four faces — the lateral, the flat, two caps — and so is their union:
        // the two statements of a ruling fold where they overlap.
        "two half disks on [0, 2] and [1, 3]" => [
            counted(&[1.5 * PI], &[4]),
            counted(&[PI / 2.0], &[4]),
            counted(&[PI / 2.0], &[4]),
            counted(&[PI / 2.0], &[4]),
        ],
        // The tall less the inner is two half disks, below and above it.
        "a half disk on [0, 3] around one on [1, 2]" => [
            counted(&[1.5 * PI], &[4]),
            counted(&[PI / 2.0, PI / 2.0], &[4, 4]),
            none(),
            counted(&[PI / 2.0], &[4]),
        ],
        // The cusps meet on the ruling `(1, 0)` alone: two bodies of five faces, as OCCT builds
        // them — each floor `y = 0` and each lateral ends on the line, which bounds a cusp on
        // either side and no material between them.
        "two cusps meeting on their tangent ruling" => [
            counted(&[1.0 - PI / 4.0, 1.0 - PI / 4.0], &[5, 5]),
            counted(&[1.0 - PI / 4.0], &[5]),
            counted(&[1.0 - PI / 4.0], &[5]),
            none(),
        ],
        // Two bodies of five faces again, the line an edge of both over the stretch they share and
        // of one beyond it.
        "two cusps under the axis, offset along it" => [
            counted(&[2.0 - PI / 2.0, 2.0 - PI / 2.0], &[5, 5]),
            counted(&[2.0 - PI / 2.0], &[5]),
            counted(&[2.0 - PI / 2.0], &[5]),
            none(),
        ],
        // The union is a block bored by a circle tangent to its floor: one body whose boundary
        // touches itself along `(1, 0)`. The two floors stay two faces there, so the line is an
        // edge four faces of one solid use.
        "two cusped blocks joined round the bore" => [
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[6.0 - PI / 2.0], &[7]),
            counted(&[6.0 - PI / 2.0], &[7]),
            none(),
        ],
        // One body whose boundary touches itself along `(1, 0)`, its floor one face across the
        // line: the laterals' edge lies inside it, and the self-touch verdict names the line.
        "two cusped blocks joined on their floor" => [
            refused(RejectReason::SelfTouchingResult),
            counted(&[6.0 - PI / 2.0], &[7]),
            counted(&[18.0 - PI / 2.0], &[9]),
            none(),
        ],
        // The two touch along the tangent ruling alone; OCCT builds the union as these two bodies.
        "a cusp beside a plane through the axis" => [
            counted(&[1.0 - PI / 4.0, 2.0 - PI / 4.0], &[5, 6]),
            counted(&[1.0 - PI / 4.0], &[5]),
            counted(&[2.0 - PI / 4.0], &[6]),
            none(),
        ],
        // Round the tangent ruling the union is the cylinder's outside whole: one body whose lateral
        // runs across the line as one face — the line's star with both sectors beyond the floor
        // kept (`Chart::slit_at`). OCCT 7.9 gives the same volumes, the union in twelve faces and
        // the second less the first in nine.
        "a cusp and the block bitten beyond its floor" => [
            counted(&[7.0 - PI / 2.0], &[9]),
            counted(&[1.0 - PI / 4.0], &[5]),
            counted(&[6.0 - PI / 4.0], &[8]),
            none(),
        ],
        // Each lateral ends on the tangent ruling at its own wall, and between the walls below the
        // line the wedge is empty: two bodies touching along the line, the chamber between the
        // station's two chords not kept. OCCT 7.9: the same volume, twelve faces.
        "two blocks bitten into a wedge's walls" => [
            counted(&[4.0 - PI / 4.0, 4.0 - PI / 4.0], &[6, 6]),
            counted(&[4.0 - PI / 4.0], &[6]),
            counted(&[4.0 - PI / 4.0], &[6]),
            none(),
        ],
        // Joined above the bore, the union is one body touching itself along the line, which four
        // of its faces use (OCCT 7.9 builds it, in sixteen faces); the cuts and the common as OCCT
        // builds them.
        "two blocks bitten into a wedge's walls, joined above the bore" => [
            refused(RejectReason::NonManifoldResultEdge),
            counted(&[4.5 - PI / 4.0], &[6]),
            counted(&[4.5 - PI / 4.0], &[6]),
            counted(&[7.0 - PI / 2.0], &[8]),
        ],
        other => panic!("no answers stated for {other}"),
    }
}

/// **A cylinder two solids share is one class read for both.** Every family, as each of the four
/// booleans, comes out as [`check`] reads its answer.
#[test]
fn a_cylinder_shared_by_two_solids_reads_both_solids() {
    for (name, build) in FAMILIES {
        for (op, want) in OPS.into_iter().zip(wants(name)) {
            let mut m = Model::new();
            let (a, b) = build(&mut m, Says::Plain);
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}

fn extrude_in(m: &mut Model, frame: SketchFrame, profile: Profile2d, dist: f64) -> Handle<Solid> {
    let OpOutput::Extrude { solid, .. } = apply(
        m,
        &Operation::Extrude {
            frame,
            profile,
            dist,
        },
    )
    .expect("the profile extrudes") else {
        unreachable!()
    };
    solid
}

/// The frame of a datum plane through `origin` with normal `normal`.
fn datum(m: &mut Model, origin: [f64; 3], normal: [f64; 3]) -> SketchFrame {
    let plane = SketchPlane::from_origin_normal(
        nacre_math::Point3::from_array(origin),
        nacre_math::Vector3::from_array(normal),
    )
    .expect("a plane");
    datum_frame(m, plane)
}

/// The frame of a datum plane through `origin` with the axes `(x, y)`.
fn frame_on(m: &mut Model, origin: [f64; 3], (x, y): ([f64; 3], [f64; 3])) -> SketchFrame {
    datum_frame(
        m,
        SketchPlane::from_axes(
            Point3::from_array(origin),
            Vector3::from_array(x),
            Vector3::from_array(y),
        ),
    )
}

/// Axes of a plane across the axis `(3, 4, 0)` with the rational lengths `1` and `5`: the frame's
/// `+u = ẑ` and `+v = (4, −3, 0)/5` are exact, so a circle's cylinder is stated in the world along
/// the frame's unit normal `(3, 4, 0)/5`.
const WORLD_345: ([f64; 3], [f64; 3]) = ([0.0, 0.0, 1.0], [4.0, -3.0, 0.0]);

/// Axes of a plane across the axis `(3, 4, 0)` whose `+u = (−4, 3, 5)` has the length `√50`: the
/// frame is a node, and its chain's fold carries a circle's axis out as the primitive direction
/// `(3, 4, 0)` — [`WORLD_345`]'s axis written five times as long. Only a circle about the frame's
/// origin has a rational world statement.
const NODE_345: ([f64; 3], [f64; 3]) = ([-4.0, 3.0, 5.0], [20.0, -15.0, 25.0]);

fn moved(m: &mut Model, solid: Handle<Solid>, isometry: nacre_exact::Isometry) -> Handle<Solid> {
    let OpOutput::Transform { solid } =
        apply(m, &Operation::Transform { solid, isometry }).expect("the motion")
    else {
        unreachable!()
    };
    solid
}

fn lift(dz: i128) -> nacre_exact::Isometry {
    let z = nacre_exact::Rat::from_int;
    nacre_exact::Isometry::translation([z(0), z(0), z(dz)])
}

/// The unit circle swept `2` up from `z = 0`, and its twin swept `2` from the same plane stated
/// facing down — the same cylinder, its axis written `−ẑ`.
fn twin_with_the_axis_reversed(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let down = datum(m, [0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
    (a, extrude_in(m, down, disk(), 2.0))
}

/// The twin sketched with its `+u` along `ŷ` — another seam direction.
fn twin_with_another_seam(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let plane = m.seeded_plane(Axis::Z);
    let frame = SketchFrame::named(
        m,
        plane,
        nacre_math::Point3::from_array([0.0, 0.0, 0.0]),
        nacre_math::Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("the frame");
    (a, extrude_in(m, frame, disk(), 2.0))
}

/// A prism on `z ∈ [0, 4]` and a twin sketched on the datum `z = 2` — another origin on the line.
fn twin_on_another_origin(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 4.0);
    let up = datum(m, [0.0, 0.0, 2.0], [0.0, 0.0, 1.0]);
    (a, extrude_in(m, up, disk(), 2.0))
}

/// A prism on `z ∈ [0, 4]` and a twin moved up `2` — another origin, carried by a motion.
fn twin_moved_up(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 4.0);
    let b = sweep(m, Says::Plain, disk(), 2.0);
    (a, moved(m, b, lift(2)))
}

/// The twin turned a quarter about the axis — a motion that turns the seam.
fn twin_turned_a_quarter(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let b = sweep(m, Says::Plain, disk(), 2.0);
    (a, moved(m, b, super::rot_iso(Axis::Z, 90)))
}

/// The twin moved up `1` and back down.
fn twin_moved_there_and_back(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let b = sweep(m, Says::Plain, disk(), 2.0);
    let b = moved(m, b, lift(1));
    (a, moved(m, b, lift(-1)))
}

/// The twin turned a quarter and then three quarters about the axis — a whole turn.
fn twin_turned_a_whole_turn(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let b = sweep(m, Says::Plain, disk(), 2.0);
    let b = moved(m, b, super::rot_iso(Axis::Z, 90));
    (a, moved(m, b, super::rot_iso(Axis::Z, 270)))
}

/// The upper half disk and its twin sketched with `+u` along `ŷ` — the half on `x ≤ 0`, its
/// corners ranked about another seam.
fn half_twin_with_another_seam(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, upper_half(), 2.0);
    let plane = m.seeded_plane(Axis::Z);
    let frame = SketchFrame::named(
        m,
        plane,
        Point3::origin(),
        Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("the frame");
    (a, extrude_in(m, frame, upper_half(), 2.0))
}

/// The upper half disk on `z ∈ [0, 2]` and its twin swept down from the datum `z = 1` facing
/// `−ẑ` — the same half on `z ∈ [−1, 1]`, its axis written the other way.
fn half_twin_with_the_axis_reversed(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, upper_half(), 2.0);
    let down = datum(m, [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]);
    (a, extrude_in(m, down, upper_half(), 2.0))
}

/// The whole disk on `z ∈ [0, 2]` and the half on `x ≤ 0` over `z ∈ [−1, 2]`, sketched with `+u`
/// along `ŷ` and moved down `1`.
fn whole_and_half_twin_with_another_seam(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, disk(), 2.0);
    let plane = m.seeded_plane(Axis::Z);
    let frame = SketchFrame::named(
        m,
        plane,
        Point3::origin(),
        Vector3::from_array([0.0, 1.0, 0.0]),
    )
    .expect("the frame");
    let b = extrude_in(m, frame, upper_half(), 3.0);
    (a, moved(m, b, lift(-1)))
}

/// The upper half disk and its twin turned a quarter about the axis — the half on `x ≤ 0`, its
/// seam turned by a motion.
fn half_twin_turned_a_quarter(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let a = sweep(m, Says::Plain, upper_half(), 2.0);
    let b = sweep(m, Says::Plain, upper_half(), 2.0);
    (a, moved(m, b, super::rot_iso(Axis::Z, 90)))
}

/// The unit circle about the axis `(3, 4, 0)` through the origin, swept `4` from there on the world
/// road ([`WORLD_345`]), and its twin swept `4` from `(1.8, 2.4, 0)` — `3` along the axis — on the
/// frame node ([`NODE_345`]): the same cylinder, its axis written five times as long, its seam and
/// origin elsewhere ([`the_scaled_twin_writes_one_axis_at_two_lengths`]).
fn twin_with_the_axis_scaled(m: &mut Model) -> (Handle<Solid>, Handle<Solid>) {
    let world = frame_on(m, [0.0, 0.0, 0.0], WORLD_345);
    let a = extrude_in(m, world, disk(), 4.0);
    let node = frame_on(m, [1.8, 2.4, 0.0], NODE_345);
    (a, extrude_in(m, node, disk(), 4.0))
}

/// **The scaled twin writes one axis at two lengths** — the population the row «the axis scaled» of
/// [`one_cylinder_under_two_handles_is_one_class`] is about: one surface, its axis `(3, 4, 0)/5` in
/// the first statement and `(3, 4, 0)` in the second, and the seams apart (a frame that writes the
/// axis at another length cannot write the seam the same way: the world road's `+u` has a rational
/// length, the node's has not).
#[test]
fn the_scaled_twin_writes_one_axis_at_two_lengths() {
    let mut m = Model::new();
    let (a, b) = twin_with_the_axis_scaled(&mut m);
    let def = |s| {
        let [h] = lateral_surfaces(&m, s)[..] else {
            panic!("one cylinder a solid")
        };
        m.world_cylinder_def(h).expect("a world statement")
    };
    let (da, db) = (def(a), def(b));
    let q = |n: i128, d: i128| Rat::new(n, d).expect("a rational");
    assert!(da.same_surface(&db), "one surface");
    assert_eq!(da.dir(), [q(3, 5), q(4, 5), q(0, 1)]);
    assert_eq!(db.dir(), [q(3, 1), q(4, 1), q(0, 1)]);
    assert!(
        !nacre_exact::parallel_rat(&da.ref_dir(), &db.ref_dir()),
        "{:?} {:?}",
        da.ref_dir(),
        db.ref_dir()
    );
}

/// **One cylinder under two handles is one class** — each twin states the operand's cylinder once
/// more, differing in one part of the statement — but the twin whose axis is written at another
/// length, which differs in its seam and origin too — (the last four with corners on the cylinder,
/// which rank from the seam of the handle that names them), lands on a handle of its own, and the
/// pair answers as the two solids do: the volumes set out from the shapes.
#[test]
fn one_cylinder_under_two_handles_is_one_class() {
    let none = || counted(&[], &[]);
    let whole = |v: f64| counted(&[v * PI], &[3]);
    // Two half disks a quarter turn apart: three quarters of the disk, and a quarter each way.
    let quarter_turned = || {
        [
            counted(&[1.5 * PI], &[5]),
            counted(&[0.5 * PI], &[5]),
            counted(&[0.5 * PI], &[5]),
            counted(&[0.5 * PI], &[5]),
        ]
    };
    let twins: [(&str, Twins, [Want; 4]); 10] = [
        (
            "the axis reversed",
            twin_with_the_axis_reversed,
            [whole(4.0), whole(2.0), whole(2.0), none()],
        ),
        (
            "another seam",
            twin_with_another_seam,
            [whole(2.0), none(), none(), whole(2.0)],
        ),
        (
            "another origin",
            twin_on_another_origin,
            [whole(4.0), whole(2.0), none(), whole(2.0)],
        ),
        (
            "moved up",
            twin_moved_up,
            [whole(4.0), whole(2.0), none(), whole(2.0)],
        ),
        (
            "turned a quarter",
            twin_turned_a_quarter,
            [whole(2.0), none(), none(), whole(2.0)],
        ),
        // `[0, 4]` and `[3, 7]` along the axis.
        (
            "the axis scaled",
            twin_with_the_axis_scaled,
            [whole(7.0), whole(3.0), whole(3.0), whole(1.0)],
        ),
        (
            "a half, another seam",
            half_twin_with_another_seam,
            quarter_turned(),
        ),
        // The same half on `z ∈ [0, 2]` and `[−1, 1]`: three high together, one high each way.
        (
            "a half, the axis reversed",
            half_twin_with_the_axis_reversed,
            [
                counted(&[1.5 * PI], &[4]),
                counted(&[0.5 * PI], &[4]),
                counted(&[0.5 * PI], &[4]),
                counted(&[0.5 * PI], &[4]),
            ],
        ),
        // The disk on `[0, 2]` and the half on `[−1, 2]`: the half's foot below the disk.
        (
            "a whole and a half, another seam",
            whole_and_half_twin_with_another_seam,
            [
                counted(&[2.5 * PI], &[5]),
                counted(&[PI], &[4]),
                counted(&[0.5 * PI], &[4]),
                counted(&[PI], &[4]),
            ],
        ),
        (
            "a half turned a quarter",
            half_twin_turned_a_quarter,
            quarter_turned(),
        ),
    ];
    for (name, build, wants) in twins {
        for (op, want) in OPS.into_iter().zip(wants) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
            assert_eq!((la.len(), lb.len()), (1, 1), "{name}");
            assert_ne!(la, lb, "{name}: two handles");
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}

/// **A motion that folds back to none is one statement** — the twin moved there and back, or
/// turned a whole turn, interns to the operand's own handle, and the pair builds as two identical
/// prisms do.
#[test]
fn a_motion_folded_back_states_the_cylinder_once() {
    let twins: [(&str, Twins); 2] = [
        ("moved there and back", twin_moved_there_and_back),
        ("turned a whole turn", twin_turned_a_whole_turn),
    ];
    for (name, build) in twins {
        for (op, want) in OPS.into_iter().zip(wants("two identical prisms")) {
            let mut m = Model::new();
            let (a, b) = build(&mut m);
            assert_eq!(
                lateral_surfaces(&m, a),
                lateral_surfaces(&m, b),
                "{name}: one handle"
            );
            let got = run(&mut m, op, a, b);
            check(&format!("{name} {op:?}"), &m, got, &want);
        }
    }
}

/// A family's first operand as `first` states it and its second as `second` does — the same two
/// solids, their cylinders under statements of their own.
fn restated(
    m: &mut Model,
    build: Build,
    (first, second): (Says, Says),
) -> (Handle<Solid>, Handle<Solid>) {
    let (a, _) = build(m, first);
    let (_, b) = build(m, second);
    (a, b)
}

/// **Every restated family is the same two solids under two handles** — the population the restated
/// scoreboard is about, checked before its answers are read. Under each statement both operands fill
/// the same volume, box and centroid as they do written plainly — the box reads an arc as its whole
/// circle and cannot tell a half from its mirror (`nacre_props::bounds`) — and each one's statement
/// differs from the plain one in the part the statement names and in no other; under every pair of
/// two statements ([`statement_pairs`]) the two operands share no cylinder handle. Each family is
/// built once per statement, all in one model — a statement interns to one handle whatever else the
/// model holds.
#[test]
fn each_restated_family_states_the_same_solids_under_two_handles() {
    let at = |says: Says| {
        STATEMENTS
            .iter()
            .position(|&s| s == says)
            .expect("a statement")
    };
    for (name, build) in FAMILIES {
        let mut m = Model::new();
        let built: Vec<[Handle<Solid>; 2]> = STATEMENTS
            .iter()
            .map(|&says| {
                let (a, b) = build(&mut m, says);
                [a, b]
            })
            .collect();
        let props = |s| {
            let volume = nacre_props::mass_props(&m, s).expect("props").volume;
            let (lo, hi) = nacre_props::bounds(&m, s).expect("bounds");
            let centroid = nacre_props::centroid(&m, s).expect("centroid");
            (volume, lo.as_array(), hi.as_array(), centroid.as_array())
        };
        let def = |s, label: &str| {
            let [h] = lateral_surfaces(&m, s)[..] else {
                panic!("{label}: one cylinder a solid")
            };
            m.world_cylinder_def(h).expect("a world statement")
        };
        for says in STATEMENTS.into_iter().skip(1) {
            for side in 0..2 {
                let label = format!("{name} {says:?} operand {side}");
                let (plain, solid) = (built[0][side], built[at(says)][side]);
                let ((pv, plo, phi, pc), (bv, blo, bhi, bc)) = (props(plain), props(solid));
                assert!((pv - bv).abs() < 1e-9, "{label}: {pv} and {bv}");
                assert_eq!((plo, phi), (blo, bhi), "{label}");
                assert!(
                    (0..3).all(|k| (pc[k] - bc[k]).abs() < 1e-9),
                    "{label}: centroid {pc:?} and {bc:?}"
                );
                let (p, r) = (def(plain, &label), def(solid, &label));
                let flipped = p.dir().map(|x| Rat::from_int(0).checked_sub(x).unwrap());
                let same = [
                    p.origin() == r.origin(),
                    p.dir() == r.dir(),
                    p.ref_dir() == r.ref_dir(),
                ];
                match says {
                    Says::Plain => unreachable!(),
                    Says::Seam => assert_eq!(same, [true, true, false], "{label}"),
                    Says::Origin => assert_eq!(same, [false, true, true], "{label}"),
                    Says::Reversed => {
                        assert_eq!(same, [true, false, true], "{label}");
                        assert_eq!(r.dir(), flipped, "{label}");
                    }
                }
            }
        }
        for (x, y) in statement_pairs() {
            let (a, b) = (built[at(x)][0], built[at(y)][1]);
            let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
            assert!(
                lb.iter().all(|h| !la.contains(h)),
                "{name} {x:?}/{y:?}: {la:?} {lb:?}"
            );
        }
    }
}

/// **A cylinder stated under two handles answers as it does under one** — every family, its first
/// operand stated `first` way and its second each other way ([`statement_pairs`]), as each of the
/// four booleans: the answers of the one-handle scoreboard ([`wants`]). Every row is run before any
/// is judged, so a reading that goes wrong shows every row it moves. One test per first statement,
/// so that the suite runs the four side by side (10 to 16 s each).
fn answers_as_under_one(first: Says) {
    let mut wrong: Vec<String> = Vec::new();
    for (name, build) in FAMILIES {
        for pair in statement_pairs().filter(|&(x, _)| x == first) {
            for (op, want) in OPS.into_iter().zip(wants(name)) {
                let label = format!("{name} {:?}/{:?} {op:?}", pair.0, pair.1);
                let row = std::panic::catch_unwind(|| {
                    let mut m = Model::new();
                    let (a, b) = restated(&mut m, build, pair);
                    let got = run(&mut m, op, a, b);
                    check(&label, &m, got, &want);
                });
                if let Err(e) = row {
                    let said = e
                        .downcast_ref::<String>()
                        .cloned()
                        .unwrap_or_else(|| label.clone());
                    wrong.push(said);
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} rows:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// [`answers_as_under_one`] with the first operand stated plainly.
#[test]
fn a_cylinder_restated_under_a_second_handle_answers_as_under_one() {
    answers_as_under_one(Says::Plain);
}

/// [`answers_as_under_one`] with the first operand's seam turned.
#[test]
fn the_first_operand_on_a_turned_seam_answers_as_under_one() {
    answers_as_under_one(Says::Seam);
}

/// [`answers_as_under_one`] with the first operand's origin at the far end of the sweep.
#[test]
fn the_first_operand_on_a_moved_origin_answers_as_under_one() {
    answers_as_under_one(Says::Origin);
}

/// [`answers_as_under_one`] with the first operand's axis reversed.
#[test]
fn the_first_operand_on_a_reversed_axis_answers_as_under_one() {
    answers_as_under_one(Says::Reversed);
}

/// A block `over` on `z ∈ [split, top]` bored through by the unit circle about `centre`, fused to
/// the block `under` on `z ∈ [0, split]` beneath it. The bore's lateral is whole — it touches no
/// wall — so where a wall of `under` meets another solid's lateral on the same cylinder, this
/// solid's own faces say nothing about it: the gate's reading of the other solid's row decides.
fn bored_over(
    m: &mut Model,
    says: Says,
    centre: [f64; 2],
    (split, top): (f64, f64),
    over: [[f64; 2]; 2],
    under: [[f64; 2]; 2],
) -> Handle<Solid> {
    let upper = crate::fixtures::cuboid(
        m,
        [over[0][0], over[0][1], split],
        [over[1][0], over[1][1], top],
    );
    let tool = sweep(m, says, circle_profile(centre, 1.0), top);
    let bored = crate::boolean(m, BoolKind::Cut, upper, tool).expect("the bore")[0];
    let lower = crate::fixtures::cuboid(
        m,
        [under[0][0], under[0][1], 0.0],
        [under[1][0], under[1][1], split],
    );
    crate::boolean(m, BoolKind::Fuse, bored, lower).expect("the block under the bore")[0]
}

/// A bored plate on `z ∈ [1, 2]` over the block `[−2, 2] × [0, 2]`, and the unit pin five high:
/// the block's wall `y = 0` cuts the pin's lateral over `z ∈ [0, 1]`, heights the pin's row holds
/// only once its range is moved into the class's statement.
fn a_wall_under_the_bore_and_the_pin(
    m: &mut Model,
    a: Says,
    b: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    (
        bored_over(
            m,
            a,
            [0.0, 0.0],
            (1.0, 2.0),
            [[-2.0, -2.0], [2.0, 2.0]],
            [[-2.0, 0.0], [2.0, 2.0]],
        ),
        sweep(m, b, disk(), 5.0),
    )
}

/// A bored plate on `z ∈ [1, 2]` over the block `[0.28, 2] × [0.6, 2]`, and the half disk on
/// `x ≥ 0`: the block's walls meet the circle at `(0.28, 0.96)` and `(0.8, 0.6)`, inside the half's
/// arc and off its complement.
fn a_corner_under_the_bore_and_the_half(
    m: &mut Model,
    a: Says,
    b: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    let right_half = profile(vec![
        line(p(0.0, 1.0), p(0.0, -1.0)),
        arc_turns(p(0.0, 0.0), p(0.0, -1.0), 2),
    ]);
    (
        bored_over(
            m,
            a,
            [0.0, 0.0],
            (1.0, 2.0),
            [[-2.0, -2.0], [2.0, 2.0]],
            [[0.8, 0.6], [2.0, 2.0]],
        ),
        sweep(m, b, right_half, 2.0),
    )
}

/// A plate on `z ∈ [3, 4]` bored about `(3, 1)` over the slab `[0, 6] × [−1, 0]`, and the lower
/// half of the bore's pin on `z ∈ [0, 3]`, lying on the slab: the slab's top `y = 0` is tangent to
/// the pin's lateral along a line inside the half's arc, which the tangency row reads off the
/// pin's face on both sides.
fn a_half_pin_lying_on_the_slab_under_the_bore(
    m: &mut Model,
    a: Says,
    b: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    let lower_half_pin = profile(vec![
        line(p(4.0, 1.0), p(2.0, 1.0)),
        arc_turns(p(3.0, 1.0), p(2.0, 1.0), 2),
    ]);
    (
        bored_over(
            m,
            a,
            [3.0, 1.0],
            (3.0, 4.0),
            [[0.0, -1.0], [6.0, 3.0]],
            [[0.0, -1.0], [6.0, 0.0]],
        ),
        sweep(m, b, lower_half_pin, 4.0),
    )
}

/// A plate on `z ∈ [2, 3]` bored by the unit circle, a slab under it and under that, on
/// `z ∈ [0, 1]`, a pin of radius `0.5` about `(1.2, 0)`; and the half disk on `x ≥ 0`, one high.
/// The small pin's lateral crosses the unit circle near `(1, 0)`, inside the half's arc and off
/// its complement — two lateral faces of different solids that meet, which the pair rule reads
/// off the half's arc about the class's axis.
fn a_small_pin_under_the_bore_and_the_half(
    m: &mut Model,
    a: Says,
    b: Says,
) -> (Handle<Solid>, Handle<Solid>) {
    let right_half = profile(vec![
        line(p(0.0, 1.0), p(0.0, -1.0)),
        arc_turns(p(0.0, 0.0), p(0.0, -1.0), 2),
    ]);
    let upper = crate::fixtures::cuboid(m, [-2.0, -2.0, 2.0], [2.0, 2.0, 3.0]);
    let tool = sweep(m, a, disk(), 3.0);
    let bored = crate::boolean(m, BoolKind::Cut, upper, tool).expect("the bore")[0];
    let slab = crate::fixtures::cuboid(m, [-2.0, -2.0, 1.0], [2.0, 2.0, 2.0]);
    let plate = crate::boolean(m, BoolKind::Fuse, bored, slab).expect("the slab")[0];
    let pin = crate::fixtures::cylinder(
        m,
        Point3::from_array([1.2, 0.0, 0.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
        0.5,
        1.0,
    )
    .solid;
    let first = crate::boolean(m, BoolKind::Fuse, plate, pin).expect("the pin under the slab")[0];
    (first, sweep(m, b, right_half, 1.0))
}

/// What a boolean came out as, for comparing two runs: the refusal's name, or each body's volume
/// (to nine places) and face count, ascending.
fn outcome(m: &Model, got: &super::scoreboard::BoolResult) -> String {
    match got {
        Err(BoolError::Rejected { reason, .. }) => format!("refused {reason:?}"),
        Err(e) => format!("{e:?}"),
        Ok(solids) => {
            let bad = nacre_validate::validate(m);
            assert!(bad.is_empty(), "{bad:?}");
            let mut bodies: Vec<(i64, usize)> = solids
                .iter()
                .map(|&s| {
                    let v = nacre_props::mass_props(m, s).expect("props").volume;
                    (
                        (v * 1e9).round() as i64,
                        m.shell(m.solid(s).outer).faces.len(),
                    )
                })
                .collect();
            bodies.sort_unstable();
            format!("{bodies:?}")
        }
    }
}

/// **What the gate reads off a restated row answers as under one handle** — four pairs whose
/// first solid's lateral on the shared cylinder is a bore's, whole and away from the contact, so
/// where its other faces meet the second solid's lateral the gate's reading of the second's row is
/// what the answer rests on: the arc about the class's axis (a second cylinder crossing it, the
/// pair rule's refusal), the range along the axis and the stretches on a tangent line (a slab the
/// half pin lies on, the self-touch verdict), and a corner edge on the lateral. Each pair is run
/// with both cylinders stated plainly (one handle) and again under every ordered pair of two
/// different statements ([`statement_pairs`]); the outcomes are the same. Planted, an arc left about its own axis builds the crossing pair
/// into a solid the mesher refuses, and a range left in its own parameter, an arc read about the
/// class's axis or a corner ranked from the class's seam lets the half pin's contact through. The
/// other two pairs are here for their answers alone, and no planted reading moves them: the wall
/// under the bore is a secant, which the gate records whatever the faces' ranges say, and the
/// corner edge's contact is found on the result.
#[test]
fn what_the_gate_reads_off_a_restated_row_answers_as_under_one() {
    type Pair = fn(&mut Model, Says, Says) -> (Handle<Solid>, Handle<Solid>);
    let pairs: [(&str, Pair); 4] = [
        (
            "a small pin under the bore and the half",
            a_small_pin_under_the_bore_and_the_half,
        ),
        (
            "a wall under the bore and the pin",
            a_wall_under_the_bore_and_the_pin,
        ),
        (
            "a corner under the bore and the half",
            a_corner_under_the_bore_and_the_half,
        ),
        (
            "a half pin lying on the slab under the bore",
            a_half_pin_lying_on_the_slab_under_the_bore,
        ),
    ];
    let mut wrong: Vec<String> = Vec::new();
    for (name, build) in pairs {
        for op in OPS {
            let plain = {
                let mut m = Model::new();
                let (a, b) = build(&mut m, Says::Plain, Says::Plain);
                let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
                assert!(lb.iter().all(|h| la.contains(h)), "{name}: one handle");
                let got = run(&mut m, op, a, b);
                outcome(&m, &got)
            };
            for (sa, sb) in statement_pairs() {
                let mut m = Model::new();
                let (a, b) = build(&mut m, sa, sb);
                let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
                assert!(lb.iter().all(|h| !la.contains(h)), "{name}: two handles");
                let got = run(&mut m, op, a, b);
                let got = outcome(&m, &got);
                if got != plain {
                    wrong.push(format!(
                        "{name} {op:?} {sa:?}/{sb:?}: {got}, one handle {plain}"
                    ));
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} rows:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// The plate `[−2, 2]²` of a frame across the axis `(3, 4, 0)` ([`WORLD_345`]), over `[0, 4]` along
/// the axis and bored there by the unit circle, fused to the block under the bore's tangent line
/// `u = −1` (the world's `z = −1`) over `[4, 7]` — every face on the world road.
fn a_bored_plate_on_a_block_along_the_axis_345(m: &mut Model) -> Handle<Solid> {
    let square = |u: [f64; 2], v: [f64; 2]| {
        Profile2d::polygon(vec![
            p(u[0], v[0]),
            p(u[1], v[0]),
            p(u[1], v[1]),
            p(u[0], v[1]),
        ])
        .expect("a rectangle")
    };
    let at_origin = frame_on(m, [0.0, 0.0, 0.0], WORLD_345);
    let plate = extrude_in(m, at_origin, square([-2.0, 2.0], [-2.0, 2.0]), 4.0);
    let at_origin = frame_on(m, [0.0, 0.0, 0.0], WORLD_345);
    let tool = extrude_in(m, at_origin, disk(), 4.0);
    let bored = crate::boolean(m, BoolKind::Cut, plate, tool).expect("the bore")[0];
    let at_four = frame_on(m, [2.4, 3.2, 0.0], WORLD_345);
    let block = extrude_in(m, at_four, square([-2.0, -1.0], [-2.0, 2.0]), 3.0);
    crate::boolean(m, BoolKind::Fuse, bored, block).expect("the block beyond the bore")[0]
}

/// The bore's unit pin over `[3, 7]` along the axis `(3, 4, 0)`: on the frame node, its axis five
/// times as long ([`NODE_345`]), or on the bore's own handle — the circle swept `7` from the origin
/// on the world road and kept on `[3, 7]`.
fn the_pin_along_the_axis_345(m: &mut Model, on_the_node: bool) -> Handle<Solid> {
    if on_the_node {
        let node = frame_on(m, [1.8, 2.4, 0.0], NODE_345);
        return extrude_in(m, node, disk(), 4.0);
    }
    let at_origin = frame_on(m, [0.0, 0.0, 0.0], WORLD_345);
    let pin = extrude_in(m, at_origin, disk(), 7.0);
    let at_three = frame_on(m, [1.8, 2.4, 0.0], WORLD_345);
    let keep = extrude_in(
        m,
        at_three,
        Profile2d::polygon(vec![p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0)])
            .expect("a square"),
        4.0,
    );
    crate::boolean(m, BoolKind::Common, pin, keep).expect("the pin's span")[0]
}

/// **A row whose axis is written at another length answers as under one handle** — the bored plate
/// on its block ([`a_bored_plate_on_a_block_along_the_axis_345`]) and the bore's pin
/// ([`the_pin_along_the_axis_345`]): the block's top is tangent to the pin along `[4, 7]`, past the
/// bore, so whether a contact is there the gate reads off the pin's range along the axis — on the
/// frame node, a range in the node's statement moved into the class's through the axis's length
/// (`CylFaceInfo::restated_in`). Each boolean is run with the pin on the bore's handle and on the
/// node, and the outcomes are the same. Planted, the length read as its sign alone narrows the pin's
/// range, and the union is refused as `NonManifoldVertex` where it is `SelfTouchingResult` (the pin
/// lies on the block along a line, inside one body).
#[test]
fn a_row_restated_through_the_axis_length_answers_as_under_one() {
    for op in OPS {
        let answer = |on_the_node: bool| {
            let mut m = Model::new();
            let a = a_bored_plate_on_a_block_along_the_axis_345(&mut m);
            let b = the_pin_along_the_axis_345(&mut m, on_the_node);
            let (la, lb) = (lateral_surfaces(&m, a), lateral_surfaces(&m, b));
            assert_eq!(la == lb, !on_the_node, "{op:?}: {la:?} {lb:?}");
            let got = run(&mut m, op, a, b);
            outcome(&m, &got)
        };
        assert_eq!(answer(true), answer(false), "{op:?}");
    }
}
