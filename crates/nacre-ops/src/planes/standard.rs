use super::*;
#[cfg(test)]
use nacre_judge::predicate::Notes;
/// Which sub-phase of [`crate::arrangement::plane_index_setup`] a [`Watch`] charges — spike instrumentation, and only
/// in a test build (see `crate::phase`).
pub(crate) enum Sub {
    TriPt3,
    Std,
    Collect,
    Edges,
    Classes,
    Dense,
}

/// Times the scope it is charged from, or does nothing at all in a release build.
pub(crate) struct Watch(#[cfg(test)] std::time::Instant);

impl Watch {
    pub(crate) fn new() -> Self {
        Watch(
            #[cfg(test)]
            std::time::Instant::now(),
        )
    }
    #[allow(unused_variables)]
    pub(crate) fn charge(self, which: Sub) {
        #[cfg(test)]
        {
            use crate::phase;
            let c = match which {
                Sub::TriPt3 => &phase::S_TRIPT3,
                Sub::Std => &phase::S_STD,
                Sub::Collect => &phase::S_COLLECT,
                Sub::Edges => &phase::S_EDGES,
                Sub::Classes => &phase::S_CLASSES,
                Sub::Dense => &phase::S_DENSE,
            };
            c.fetch_add(
                self.0.elapsed().as_nanos() as u64,
                std::sync::atomic::Ordering::Relaxed,
            );
        }
    }
}
/// **How precisely this operation's rotated definitions must be realized.**
///
/// The judges' error radius is `C · 2⁻ᵖʳᵉᶜ`, and `C` belongs to the model — it grows about one
/// bit per turn of rotation history and with the coordinate magnitudes. A fixed precision
/// therefore decides, silently, how long a model's history may be: at 256 bits a solid turned 245
/// times stops building, with a reject that names a symptom rather than the cause. So the
/// precision is read off the model instead.
///
/// The target is the **coincidence precision**: two things closer than this are treated as
/// coincident, and the kernel will only say so once it has *proved* the separation is below it.
/// Its default is derived rather than chosen —
///
/// - `output_precision = scale · 2⁻⁵²`, the finest distinction the `f64` coordinates this kernel
///   emits can carry. Below it nothing survives export, so distinguishing is meaningless.
/// - `coincidence_precision = output_precision · 2⁻¹²⁸`, two whole words further down. Erring low
///   only costs bits, while erring high merges features that were genuinely apart, so the
///   asymmetry says push it down; and a word is the natural unit because astro-float allocates
///   whole words anyway.
///
/// `scale` is the largest coordinate magnitude in either operand, taken over the whole table so
/// the result does not depend on traversal order (replay must reproduce it exactly).
///
/// The precision that reaches the target is then [`nacre_judge::judge_precision`]'s to compute;
/// [`JUDGE_PREC_CAP`] is where the kernel stops and says so instead, and [`CLIMB_HEADROOM`] is
/// what a single hard judgement may spend on top of it.
///
/// **Only the coincidence limit is a candidate for a setting.** Everything else here — the
/// precision, the cap, the headroom — is derived from it and from the model, because a bit count
/// means a different physical thing in every model ("256 bits" is `1e-76` for a solid turned once
/// and `1e+15` for one turned three hundred times).
/// The plane data of a row, `None` for a cylinder — the kind filter the plane-only sweeps
/// share.
#[inline]
pub(crate) fn plane_of(r: &FaceRow) -> Option<&FaceInfo> {
    match r {
        FaceRow::Plane(p) => Some(p),
        FaceRow::Cylinder(_) => None,
    }
}

// A cylinder row contributes no witness points to the standard — structurally right, not an
// omission: an axis-aligned-grade rational cylinder is the `rotated == false` case (its exact
// def needs no high-precision realization); the rotated-cylinder story is CIP's.
pub(crate) fn standard_for(rows: &[FaceRow]) -> Standard {
    // ★ **A face that was never moved contributes exactly nothing, so it is not asked.**
    //
    // Its `tri_pt3` are `WitnessPoint::at_nearest` of the plane's own rational points and the chain
    // is empty: an f64-representable point (base = `mantissa · 2^exp`, a power-of-two denominator
    // and a numerator within `Rat`'s 127 bits) realizes as an *exact* interval, and a decimal one
    // carries the ½-ulp bound it was stated with — nothing a replay could add. The realization has
    // no rotation error to report (`an_exact_point_demands_no_precision` in `nacre-judge`, and the
    // const assert at `TRIAL_PREC` that keeps it true).
    //
    // So the loop below used to spend a full high-precision replay per point to compute a zero —
    // measured, an axis-aligned 60-fin fold did that 24,120 times for 15.7ms and a `worst` of
    // exactly `Mag::ZERO`. The same shape was removed one level down when a stated zero replaced
    // `WitnessPoint::at`'s measurement for these points ("nine BigFloat operations to compute a
    // zero").
    //
    // `max` over the empty set is `Mag::ZERO`, which is the right answer for a model with no
    // rotation history — `precision_for` reads that as "nothing to size" and returns `TRIAL_PREC`.
    let worst = worst_trial(
        rows.iter()
            .filter_map(plane_of)
            .filter(|p| p.rotated)
            .flat_map(|p| p.tri_pt3.iter()),
    );
    // `scale`, by contrast, is every point's business: it is the model's size, and an unmoved face
    // is as far from the origin as any other.
    standard_from(
        rows.iter()
            .filter_map(plane_of)
            .flat_map(|p| p.tri_pt3.iter()),
        worst,
    )
}

/// **How deep a model may be before the operation is rejected instead.**
///
/// Not a resolution limit — the arithmetic is correct at any depth — but a **cost** limit, so it
/// is set from measured cost. A judgement's realization is quadratic-ish in the precision, and the
/// cap is placed where a single boolean's judging stays in the seconds rather than the minutes:
/// 4096 bits covers a rotation history of roughly four thousand turns (measured: `C` grows one bit
/// per turn), which is far past any real model, and a model that does exceed it is told *why*
/// rather than handed a wrong answer or an unbounded wait.
pub(crate) const JUDGE_PREC_CAP: usize = 4096;

/// **How thin a witness the kernel will still judge**, expressed as the bits a single judgement
/// may ask for *beyond* what the model itself needed.
///
/// This is a **separate budget from [`JUDGE_PREC_CAP`], and it has to be.** Sharing one absolute
/// ceiling would mean a deeply-turned model — already near the cap — leaves a hard judgement no
/// room at all, so the same sliver would be judged in a fresh model and abandoned in a turned one.
/// The model's depth and a judgement's difficulty are different quantities; only the second
/// belongs here.
///
/// It has a physical reading. A judgement's uncertainty is `(C / |cofactor|) · 2⁻ᵖʳᵉᶜ`, and the
/// model already chose `prec` so that `C · 2⁻ᵖʳᵉᶜ` clears the coincidence limit; what is left is
/// `log₂(1 / cofactor)` — the **thinness of the witness**, a needle triangle or three planes that
/// almost share a line. Two words says: a witness up to `2¹²⁸` (≈ 3·10³⁸) times more degenerate
/// than the model's own size is still judged to the end.
///
/// Two words, and not a measured number, for the same reason the coincidence limit is two words
/// below the output resolution: the error is asymmetric. Too small abandons a judgement that had
/// an answer; too large only spends bits. And measurement says there is nothing to tune — across
/// the rotation corpus and models turned 100 and 800 times, **no judgement asked for even one bit
/// beyond the model's own precision** (measured with the headroom forced to zero).
pub(crate) const CLIMB_HEADROOM: usize = 128;

/// A judging context over a hand-built table, for fixtures.
///
/// The standard is the derived default for a unit-scale model, and the collector is leaked so a
/// fixture is a one-liner — a handful of `Vec`s per test run, and nothing reads them. A fixture
/// that *does* want the evidence builds its own [`Notes`] and calls [`Judge::new`].
#[cfg(test)]
pub(crate) fn test_judge<W>(planes: &[W]) -> Judge<'_, W> {
    let notes: &'static Notes = Box::leak(Box::new(Notes::new()));
    Judge::new(
        planes,
        Standard {
            prec: 256,
            coincidence: Mag::pow2(-180),
            scale: Mag::of(1.0),
            cap: 256 + CLIMB_HEADROOM,
        },
        notes,
    )
}

/// [`standard_for`] over a bare set of definitions, with no plane table — **a test helper.**
///
/// It once served witness selection in `rotated_vertex`; that consumer is gone, and splitting
/// [`worst_trial`] out of [`standard_from`] is what surfaced it. Kept because a fixture that asks
/// "what precision does *this* point demand" wants exactly the two halves in order, and spelling
/// them out at every call site says less than the name does.
#[cfg(test)]
pub(crate) fn standard_for_points<'a>(
    pts: impl IntoIterator<Item = &'a WitnessPoint> + Clone,
) -> Standard {
    standard_from(pts.clone(), worst_trial(pts))
}

/// **The realization depth this set of definitions demands** — `max` over their trial bounds.
///
/// Split from [`standard_from`] because the two halves have nothing in common but the answer: this
/// one is **all of the cost** (a full high-precision replay per point), and the other is f64
/// arithmetic on already-known numbers. Keeping them apart is what lets a caller that already knows
/// this maximum skip straight to the second half.
///
/// **This is where a boolean spends most of what is left after the arrangement went parallel**
/// (measured: 76% of setup, and setup is 43% of the largest booleans once the trace is off the
/// critical path). Each point's trial realization is independent and they combine by **maximum**,
/// which is associative and exact — so evaluating them across cores cannot move the answer the way
/// a reassociated sum would.
fn worst_trial<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>) -> Mag {
    let pts: Vec<&WitnessPoint> = pts.into_iter().collect();
    let bounds = crate::par::map_range(pts.len(), |i| nacre_judge::trial_bound(pts[i]));
    bounds
        .into_iter()
        .fold(Mag::ZERO, |w, b| if w.lt(b) { b } else { w })
}

/// The standard for points whose worst trial bound is already known: `scale` off the f64
/// coordinates, the coincidence limit derived from it, and the precision that reaches it.
fn standard_from<'a>(pts: impl IntoIterator<Item = &'a WitnessPoint>, worst: Mag) -> Standard {
    let mut scale = 1.0f64;
    for p in pts {
        for c in p.coord() {
            scale = scale.max(c.abs());
        }
    }
    let scale = Mag::of(scale);
    let output_precision = scale.times(Mag::pow2(-52));
    let coincidence = output_precision.times(Mag::pow2(-128));
    let prec = nacre_judge::precision_for(worst, coincidence);
    Standard {
        prec,
        coincidence,
        scale,
        // Relative to this model's own depth — see [`CLIMB_HEADROOM`]. The absolute ceiling that
        // leaves is `JUDGE_PREC_CAP + CLIMB_HEADROOM`, since a model deeper than the first is
        // rejected before any judging starts.
        cap: prec + CLIMB_HEADROOM,
    }
}
