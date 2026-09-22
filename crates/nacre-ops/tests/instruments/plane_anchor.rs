//! **Which point states a plane is already a choice — and the judge does not read it.**
//!
//! A plane's f64 cache stores an *anchor*: `Plane::coefficients()` is `[raw, −raw·origin]`, so `d`
//! is whatever that dot product rounds to at whichever point the pusher happened to hold. And
//! `Model::push_plane` interns by `(name, motion)` — on a hit the newcomer's cache is **discarded**
//! and the first pusher's survives. So the moment a plane exists before the operation that would
//! have created it, that operation inherits someone else's `d`.
//!
//! A volume was measured going wrong when `extrude` was made to issue its base-cap surface up
//! front, and the tempting response is "hand coefficients, not a handle". But `build_prism`
//! still takes a surface handle and **pad/pocket passes one every day**
//! — the two roads are the same road: `Some(h)` picks its orientation from `n_h·(−N) > 0`, and the
//! `None` road gets `flipped = n_h·(−N) < 0` back from the intern and flips. So "pass points, not a
//! handle" is not an escape. It never was.
//!
//! This file measures what is actually at stake, because a datum-plane operation
//! makes "the plane already exists" ordinary rather than exotic.
//! Four questions, in the order that lets the cheap one end the enquiry:
//!
//! 1. **Do two anchors even disagree?** (`which_point_states_a_plane_changes_its_stored_d`)
//! 2. **Does the judge read the part that disagrees?** (`no_anchor_lets_a_tilted_plane_carry_…`)
//! 3. **Does the model move?** (`a_pre_pushed_plane_does_not_move_the_model`)
//! 4. **Do the two models' surface caches agree?**
//!    (`which_surface_caches_two_anchors_leave_disagreeing`)
//!
//! ★★ The fourth is the only one that **sees** the disagreement. The
//! first three all answer about things a boolean derives — topology, coordinates, volume — and
//! those are bit-identical whichever anchor states the plane. Without the fourth the cache
//! itself is never compared, and nothing in this repository witnesses the defect.
//!
//! The population is census's `wf` family — a fully tilted, exactly-orthonormal *decimal* frame.
//! It is not chosen for being exotic but for being the only kind that can show anything: see
//! [`an_axis_aligned_plane_is_anchor_blind`].

use nacre_exact::Axis;
use nacre_exact::Rat;
use nacre_geom::Plane;
use nacre_math::{Point2, Point3, Vector3};
use nacre_ops::SketchFrame;
use nacre_ops::{BoolKind, OpOutput, Operation, Profile2d, SketchPlane, apply};
use nacre_topo::Model;

use crate::fixtures::datum_frame;

/// **Which surface caches do two anchors leave disagreeing?** — the comparison this file never
/// made, and the one that sees what the other three miss.
///
/// The two models state **one plane with one truth** (`wf_points()`) and differ only in the f64
/// anchor the pre-statement plants; interning hands the datum the survivor's cache, so the
/// disagreement is carried, not created. Handles 0–2 are the seeded world planes (their anchor is
/// the origin either way) and **handle 3 is the `wf` datum itself** — which is why the answer is
/// `[3]` rather than something larger.
///
/// ☑ **Both rows are `[]` because the anchor is derived from the truth's first point.** Without
/// the derivation `worst_ring` is `[3]`, and removing it turns this red before anything else.
/// ★ The `sketch_origin` row is the control that ships with the fixture: it is `[]` **today**, so
/// the file shows what a lucky pass looks like beside the case that actually disagrees — the same
/// service `an_axis_aligned_plane_is_anchor_blind` does for the other questions.
/// ★ Counts are asserted **before** the element comparison: zipping two different lengths would
/// silently stop at the shorter one, which is the trap `a_pre_pushed_plane_does_not_move_the_model`
/// already records against itself.
#[test]
fn which_surface_caches_two_anchors_leave_disagreeing() {
    let plain = wf_model(None);
    let sp = wf_plane();
    let bits = |m: &Model| -> Vec<[u64; 10]> {
        (0..m.surface_count() as u32)
            .filter_map(|i| m.surface_handle_at(i))
            .map(|h| match m.surface_cache(h) {
                nacre_geom::Surface::Plane(p) => {
                    let o = p.origin().as_array();
                    let n = p.normal().as_array();
                    let c = p.coefficients();
                    [
                        o[0].to_bits(),
                        o[1].to_bits(),
                        o[2].to_bits(),
                        n[0].to_bits(),
                        n[1].to_bits(),
                        n[2].to_bits(),
                        c[0].to_bits(),
                        c[1].to_bits(),
                        c[2].to_bits(),
                        c[3].to_bits(),
                    ]
                }
                nacre_geom::Surface::Cylinder(cy) => {
                    let o = cy.axis().origin().as_array();
                    let d = cy.axis().direction().as_array();
                    let r = cy.ref_dir().as_array();
                    [
                        o[0].to_bits(),
                        o[1].to_bits(),
                        o[2].to_bits(),
                        d[0].to_bits(),
                        d[1].to_bits(),
                        d[2].to_bits(),
                        r[0].to_bits(),
                        r[1].to_bits(),
                        r[2].to_bits(),
                        cy.radius().to_bits(),
                    ]
                }
            })
            .collect()
    };
    for (what, anchor, want) in [
        ("sketch_origin", sp.origin(), &[][..]),
        (
            "worst_ring",
            wf_ring_point(WF_RING[1][0], WF_RING[1][1]),
            &[][..],
        ),
    ] {
        let stated = wf_model(Some(anchor));
        let (a, b) = (bits(&plain), bits(&stated));
        assert_eq!(
            a.len(),
            b.len(),
            "{what}: the two models hold different surface counts"
        );
        let diff: Vec<usize> = (0..a.len()).filter(|&i| a[i] != b[i]).collect();
        assert_eq!(
            diff.as_slice(),
            want,
            "{what}: which surface caches disagree has moved — see this test's doc before \
             re-pinning, the value is transitional"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// The `wf` fixture — census.rs's tilted-decimal-frame family, restated
// ─────────────────────────────────────────────────────────────────────────────

const WF_ORIGIN: [f64; 3] = [0.1234567890123456, 0.2345678901234567, 0.3456789012345678];
const WF_X: [f64; 3] = [0.6, 0.8, 0.0];
const WF_Y: [f64; 3] = [-0.48, 0.36, 0.8];
const WF_RING: [[f64; 2]; 4] = [
    [0.1111111111111111, 0.1234567890123456],
    [4.123456789012345, 0.2345678901234567],
    [3.9876543210987654, 3.1234567890123459],
    [0.2222222222222222, 2.765432109876543],
];
const WF_DIST: f64 = 2.5;

fn wf_plane() -> SketchPlane {
    SketchPlane::from_axes(
        Point3::from_array(WF_ORIGIN),
        Vector3::from_array(WF_X),
        Vector3::from_array(WF_Y),
    )
}

fn wf_profile() -> Profile2d {
    Profile2d::polygon(WF_RING.iter().map(|&p| Point2::from_array(p)).collect())
        .expect("the wf ring is a fair profile")
}

fn lift(v: [f64; 3]) -> [Rat; 3] {
    v.map(|c| Rat::from_decimal(c).expect("the wf constants are inside the decimal window"))
}

/// The three exact points `SketchPlane::from_axes` records: `[o, o + x, o + y]`, added in `Rat`.
///
/// Restated here rather than read through a new accessor — **stating them is the experiment.**
/// The datum-plane operation this file is measuring for will hand exactly these to `push_plane`.
fn wf_points() -> [[Rat; 3]; 3] {
    let (o, x, y) = (lift(WF_ORIGIN), lift(WF_X), lift(WF_Y));
    let add = |a: [Rat; 3], b: [Rat; 3]| {
        core::array::from_fn(|k| a[k].checked_add(b[k]).expect("no overflow on wf constants"))
    };
    [o, add(o, x), add(o, y)]
}

/// A profile point placed on the plane, in exact rationals, then realized once — the same road
/// `exact::prism_rings` takes, so these are the prism's actual base-ring points.
fn wf_ring_point(u: f64, v: f64) -> Point3 {
    let (o, x, y) = (lift(WF_ORIGIN), lift(WF_X), lift(WF_Y));
    let (uu, vv) = (
        Rat::from_decimal(u).expect("in window"),
        Rat::from_decimal(v).expect("in window"),
    );
    let term = |a: Rat, b: Rat| a.checked_mul(b).expect("no overflow on wf constants");
    let sum = |a: Rat, b: Rat| a.checked_add(b).expect("no overflow on wf constants");
    Point3::from_array(core::array::from_fn(|k| {
        sum(sum(o[k], term(uu, x[k])), term(vv, y[k])).to_f64()
    }))
}

fn ulps(a: f64, b: f64) -> i64 {
    (a.to_bits() as i64 - b.to_bits() as i64).abs()
}

/// `coefficients()[3]` for a plane anchored at `p` with the sweep-facing normal the kernel uses.
///
/// ★ **`−normal`, always.** Pushing `+normal` would change `flipped` as well, and then a moved
/// result could not be attributed to `d` alone. With `−normal` on both sides the intern reports
/// `flipped == false` and the anchor is the only variable.
fn d_at(anchor: Point3, sp: &SketchPlane) -> f64 {
    Plane::from_point_normal(anchor, -sp.normal())
        .expect("the wf normal is nonzero")
        .coefficients()[3]
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Do two anchors even disagree?
// ─────────────────────────────────────────────────────────────────────────────

/// **They do — by up to 22 ulps.**
///
/// The anchors compared here are not "the datum's" versus "the prism's". They are the *sketch
/// origin* and each of the *four ring points*, and `build_prism` anchors at the **oriented** ring's
/// first point — which reversing the ring changes. So which `d` a *producer* hands in depends on a
/// winding decision made inside the prism builder, with no datum operation in sight.
///
/// ⚠ **The stored `d` does not inherit this.** This test builds its planes with
/// `Plane::from_point_normal` and never pushes one, so it measures two anchors, two `d`s.
/// The *stored* `d` does not inherit the producer's choice: `push_plane_raw` derives the anchor
/// from the truth, so the model keeps one `d` per plane however many anchors are offered;
/// `which_surface_caches_two_anchors_leave_disagreeing` is where that is asserted.
#[test]
fn which_point_states_a_plane_changes_its_stored_d() {
    let sp = wf_plane();
    let origin_d = d_at(sp.origin(), &sp);
    let mut worst = 0;
    for (i, &[u, v]) in WF_RING.iter().enumerate() {
        let d = d_at(wf_ring_point(u, v), &sp);
        let n = ulps(d, origin_d);
        worst = worst.max(n);
        println!("stat anchor_d ring[{i}] d={d:?} ulps_from_sketch_origin={n}");
    }
    println!("stat anchor_d sketch_origin d={origin_d:?} worst_ring_ulps={worst}");
    assert!(
        worst > 0,
        "the premise of this file is that anchors disagree on a tilted plane; they did not, so \
         either the fixture stopped being tilted or `coefficients` stopped reading the anchor"
    );
}

/// **The control that says why the population had to be tilted.**
///
/// `d = −raw·origin`. When `raw` has a single nonzero component that is one exact multiply, so two
/// anchors on the plane — which by definition agree in that coordinate — produce **bit-identical**
/// `d`. An axis-aligned fixture, however far apart its anchors, measures nothing at all. (This is
/// the trap the plan for this work walked into three times: the obvious fixture is the blind one.)
#[test]
fn an_axis_aligned_plane_is_anchor_blind() {
    let sp = SketchPlane::world_xy().with_origin(Point3::from_array([0.0, 0.0, 0.5]));
    let near = d_at(Point3::from_array([0.0, 0.0, 0.5]), &sp);
    let far = d_at(Point3::from_array([50.0, -37.25, 0.5]), &sp);
    assert_eq!(
        near.to_bits(),
        far.to_bits(),
        "an axis-aligned plane's d is one exact multiply; anchors cannot disagree"
    );
    println!("stat anchor_d axis_aligned near={near:?} far={far:?} ulps=0");
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Does the judge read the part that disagrees?
// ─────────────────────────────────────────────────────────────────────────────

/// **No — and not by luck: `reconcile`'s net is already fully engaged on this population.**
///
/// `WorkingPlane::reconcile` carries a plane's f64 coefficients into the predicates only when
/// [`Plane::spans_exactly`] says they describe the very plane the witness triangle spans, tested in
/// exact expansion arithmetic. For a tilted decimal plane `d` is a rounded sum of products, so the
/// test fails **for every anchor** — the coefficients are never carried, and swapping anchors
/// cannot change which road the judgment takes.
///
/// The positive control is the point of the second half: "false everywhere" is only informative
/// once the instrument is shown to be able to say true.
#[test]
fn no_anchor_lets_a_tilted_plane_carry_its_coefficients() {
    let sp = wf_plane();
    let tri = [
        wf_ring_point(WF_RING[0][0], WF_RING[0][1]),
        wf_ring_point(WF_RING[1][0], WF_RING[1][1]),
        wf_ring_point(WF_RING[2][0], WF_RING[2][1]),
    ];
    let mut anchors: Vec<(String, Point3)> = vec![("sketch_origin".into(), sp.origin())];
    for (i, &[u, v]) in WF_RING.iter().enumerate() {
        anchors.push((format!("ring[{i}]"), wf_ring_point(u, v)));
    }
    for (name, a) in &anchors {
        let pl = Plane::from_point_normal(*a, -sp.normal()).expect("nonzero normal");
        let spans = pl.spans_exactly(tri);
        println!("stat anchor_carry {name} spans_exactly={spans}");
        assert!(
            !spans,
            "{name}: a tilted decimal plane is not supposed to carry its coefficients — if this \
             now passes, the anchor became load-bearing for the judging path and the enquiry in \
             this file must be redone"
        );
    }

    // ★ Positive control: the instrument can say `true`, so "false everywhere" above is a
    // measurement and not a dead probe.
    let flat = Plane::from_point_normal(
        Point3::from_array([0.0, 0.0, 1.0]),
        Vector3::from_array([0.0, 0.0, 1.0]),
    )
    .expect("nonzero normal");
    let flat_tri = [
        Point3::from_array([0.0, 0.0, 1.0]),
        Point3::from_array([2.0, 0.0, 1.0]),
        Point3::from_array([0.0, 3.0, 1.0]),
    ];
    assert!(
        flat.spans_exactly(flat_tri),
        "positive control: an axis-aligned plane does carry its coefficients"
    );
    println!("stat anchor_carry positive_control spans_exactly=true");
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Does the model move?
// ─────────────────────────────────────────────────────────────────────────────

/// Store sizes and the whole topology by index — everything but coordinates.
fn topo_sig(m: &Model) -> Vec<String> {
    let mut out = vec![format!(
        "len {} {} {} {} {} {}",
        m.vertex_count(),
        m.edge_count(),
        m.face_count(),
        m.shell_count(),
        m.solid_count(),
        m.surface_count()
    )];
    let mut i = 0u32;
    while let Some(h) = m.edge_handle_at(i) {
        i += 1;
        let e = m.edge(h);
        out.push(format!(
            "e{} {},{} {},{}",
            h.index(),
            e.surfaces[0].index(),
            e.surfaces[1].index(),
            e.vertices[0].index(),
            e.vertices[1].index()
        ));
    }
    let mut i = 0u32;
    while let Some(h) = m.face_handle_at(i) {
        i += 1;
        let f = m.face(h);
        let loops: Vec<String> = std::iter::once(&f.outer)
            .chain(f.inner.iter())
            .map(|lp| {
                lp.half_edges
                    .iter()
                    .map(|he| format!("{}{}", he.edge.index(), if he.forward { "+" } else { "-" }))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect();
        out.push(format!(
            "f{} s{} {:?} {}",
            h.index(),
            f.surface.index(),
            f.orientation,
            loops.join(" | ")
        ));
    }
    out.push(
        m.live_solids()
            .iter()
            .map(|h| h.index().to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    out
}

/// Build the `wf` prism and fuse a second body onto it, optionally stating the sketch plane
/// **before** the extrude — which is exactly what a datum-plane operation will do.
///
/// ★ The boolean is not decoration. A prism's vertices come from its own rational swept ring, so a
/// plain extrude's geometry cannot depend on the stored plane at all; `three_planes` — the only
/// consumer of a plane's f64 cache — runs in the arrangement. A fixture that stops at the extrude
/// measures nothing.
fn wf_model(pre_state_at: Option<Point3>) -> Model {
    let mut m = Model::new();
    let sp = wf_plane();
    // ★★ **Both arms now state the plane** — an extrude names its plane, so "without a
    // pre-existing plane" is not expressible. What still varies, and what this
    // file measures, is **which point anchors the cache**: `None` lets the datum anchor at the
    // caller's sketch origin, `Some(p)` plants the same plane at `p` first so the datum interns
    // onto it and inherits that anchor. Two production roads, one of them deliberately worse.
    if let Some(anchor) = pre_state_at {
        let (_h, flipped) = m.push_plane(
            Plane::from_point_normal(anchor, -sp.normal()).expect("nonzero normal"),
            wf_points(),
            None,
        );
        assert!(
            !flipped,
            "the pre-statement must face the way the base cap does, or the experiment conflates \
             the anchor with the orientation"
        );
    }
    let __frame0 = datum_frame(&mut m, sp);
    let OpOutput::Extrude { solid: a, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __frame0,
            profile: wf_profile(),
            dist: WF_DIST,
        },
    )
    .expect("the wf base prism") else {
        unreachable!()
    };
    // A world-axis block through the tilted prism: transverse contact, so the arrangement has to
    // intersect planes from both frames and read their f64 caches.
    let __w0 = SketchFrame::world(&m, Axis::Z);
    let OpOutput::Extrude { solid: b, .. } = apply(
        &mut m,
        &Operation::Extrude {
            frame: __w0,
            profile: Profile2d::polygon(vec![
                Point2::from_array([0.5, 0.5]),
                Point2::from_array([2.5, 0.5]),
                Point2::from_array([2.5, 2.5]),
                Point2::from_array([0.5, 2.5]),
            ])
            .expect("a square"),
            dist: 3.0,
        },
    )
    .expect("the cutting block") else {
        unreachable!()
    };
    apply(
        &mut m,
        &Operation::Boolean {
            kind: BoolKind::Cut,
            a,
            b,
        },
    )
    .expect("the wf cut");
    m.rebuild_adjacency();
    m
}

/// ★★ **The gate: which point anchors a plane's cache does not move the model.**
///
/// An extrude *names* its plane, so "the plane did not exist yet" is not a thing that can happen
/// — both arms below state it. What varies is the anchor: the datum's own (the caller's sketch
/// origin) against one planted deliberately at the worst ring point. Both are production roads.
///
/// Topology must match exactly; coordinates and volume are held to the project's ε (model size ×
/// `2⁻⁴⁰`), and the deviation is *printed* rather than merely passed, because a gate that only says
/// yes cannot show drift accumulating.
///
/// If this ever fails, the repair is not "pass points instead of a handle" — that road is the same
/// road (see the module doc). It is to derive the cache's anchor from the definition, which is
/// what ships.
///
/// ⚠★★★ **But not from the foot of perpendicular** — "the
/// minimum-norm point, therefore the smallest rounding available" — which was measured out:
/// the foot minimizes the rounding of `d`, not the residual at the face's own points, and it puts
/// the anchor near the *world origin* while every consumer of an anchor (conditioning, STEP's
/// required point, tess's chart) wants it near the **face**. Wired, it moved 32 census result
/// rows, closed an exact-coefficient road, and overflowed `i128` on 8 planes because the
/// projection squares the coefficients. What shipped is the truth's **first point**: 0 rows moved,
/// 61 caches changed. See `which_surface_caches_two_anchors_leave_disagreeing` below.
#[test]
fn a_pre_pushed_plane_does_not_move_the_model() {
    let plain = wf_model(None);
    let sp = wf_plane();

    // Two pre-statements. The first is what a datum operation will actually do; the second exists
    // so the gate has teeth.
    //
    // ★★ **The realistic anchor alone would be a lucky pass.** Measured: stating the plane at the
    // sketch origin *does* change the stored `Plane`'s anchor point (from the oriented ring's first
    // point to the origin) — but on this fixture the two round to the **same `d`**, so nothing
    // downstream could move whatever the answer was. A gate that only ever sees a zero-difference
    // input is not a gate. `ring[1]` is 22 ulps away in `d`, which is the largest disagreement this
    // plane admits, so the control measures the real ceiling rather than a convenient corner.
    for (what, anchor) in [
        ("sketch_origin", sp.origin()),
        ("worst_ring", wf_ring_point(WF_RING[1][0], WF_RING[1][1])),
    ] {
        let stated = wf_model(Some(anchor));

        let (sa, sb) = (topo_sig(&plain), topo_sig(&stated));
        for (x, y) in sa.iter().zip(&sb) {
            assert_eq!(x, y, "{what}: topology moved: {x:?} vs {y:?}");
        }
        assert_eq!(
            sa.len(),
            sb.len(),
            "{what}: topology moved (different cell counts)"
        );

        let size = plain
            .live_solids()
            .iter()
            .filter_map(|&s| nacre_props::bounds(&plain, s).ok())
            .map(|(lo, hi)| (hi - lo).as_array().iter().cloned().fold(0.0f64, f64::max))
            .fold(0.0f64, f64::max);
        let eps = size * f64::powi(2.0, -40);

        let mut worst_coord = 0.0f64;
        // `zip` stopped at the shorter side; the min keeps that exactly.
        let n = plain.vertex_count().min(stated.vertex_count()) as u32;
        for i in 0..n {
            let (ha, hb) = (
                plain.vertex_handle_at(i).expect("in range"),
                stated.vertex_handle_at(i).expect("in range"),
            );
            let (p, q) = (plain.vertex_point(ha), stated.vertex_point(hb));
            for k in 0..3 {
                worst_coord = worst_coord.max((p.as_array()[k] - q.as_array()[k]).abs());
            }
        }

        let vol = |m: &Model| -> f64 {
            m.live_solids()
                .iter()
                .map(|&s| {
                    nacre_props::mass_props(m, s)
                        .map(|p| p.volume)
                        .unwrap_or(0.0)
                })
                .sum()
        };
        let (va, vb) = (vol(&plain), vol(&stated));
        let rel = if va == 0.0 {
            0.0
        } else {
            (va - vb).abs() / va.abs()
        };

        println!(
            "stat pre_pushed_plane {what} d_ulps={} volume={va:.17e} vs {vb:.17e} \
             rel_dev={rel:.3e} max_coord_dev={worst_coord:.3e} eps={eps:.3e} bit_identical={}",
            ulps(d_at(anchor, &sp), d_at(plain_base_anchor(&plain), &sp)),
            va.to_bits() == vb.to_bits()
        );
        assert!(
            worst_coord <= eps,
            "{what}: coordinates moved past ε: {worst_coord:.3e} > {eps:.3e}"
        );
        assert!(
            (va - vb).abs() <= eps * size * size,
            "{what}: volume moved past ε: |{va} − {vb}| > {:.3e}",
            eps * size * size
        );
    }
}

/// The anchor `build_prism` chose for the base cap — read back so the report can say how far the
/// pre-statement actually moved `d`, instead of assuming which ring point won the winding.
fn plain_base_anchor(m: &Model) -> Point3 {
    let f = m.shell(m.solid(m.live_solids()[0]).outer).faces[0];
    match m.surface_cache(m.face(f).surface) {
        nacre_geom::Surface::Plane(p) => p.origin(),
        nacre_geom::Surface::Cylinder(_) => unreachable!("the wf base cap is planar"),
    }
}
