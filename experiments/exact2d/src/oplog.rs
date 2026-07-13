//! H5 — operation log: `Document = Vec<Op>` (the truth) + a derived model that
//! is a deterministic **replay** of the log (DNA ③ "model = replay", §7).
//!
//! This first piece establishes replay determinism with **raw-Handle** references
//! (an op refers to an earlier element by its handle). The provenance A/B — raw
//! Handle vs semantic lineage, for edit resilience / TNP — is the next piece; it
//! is what lets a *dimension* edit survive and a *topology* edit be detected
//! (§7). Semantic lineage lives here (not in H2/H3) because it references *ops*.

use crate::Rat;
use crate::share::{Handle, Point2, Segment, Sketch};

/// A logged operation. References to earlier elements are raw handles for now.
#[derive(Clone, Copy, Debug)]
pub enum Op {
    AddPoint {
        x: Rat,
        y: Rat,
    },
    AddSegment {
        a: Handle<Point2>,
        b: Handle<Point2>,
    },
}

/// Owns the log (truth) and the derived model (replay result / live state).
pub struct Document {
    log: Vec<Op>,
    model: Sketch,
}

impl Document {
    pub fn new() -> Self {
        Document {
            log: Vec::new(),
            model: Sketch::new(),
        }
    }

    pub fn log(&self) -> &[Op] {
        &self.log
    }

    pub fn model(&self) -> &Sketch {
        &self.model
    }

    /// Apply = record the op AND update the live model, returning the new handle.
    pub fn add_point(&mut self, x: Rat, y: Rat) -> Handle<Point2> {
        let h = self.model.add_point(x, y);
        self.log.push(Op::AddPoint { x, y });
        h
    }

    /// Apply a segment referencing existing point handles.
    pub fn add_segment(&mut self, a: Handle<Point2>, b: Handle<Point2>) -> Handle<Segment> {
        let h = self.model.add_segment(a, b);
        self.log.push(Op::AddSegment { a, b });
        h
    }

    /// Rebuild a fresh model by folding the log — DNA ③ "model = replay". With an
    /// append-only, deterministic-order model, raw-handle references resolve to
    /// the same elements every time.
    pub fn replay(&self) -> Sketch {
        replay_ops(&self.log)
    }
}

/// Fold a raw-reference op log into a model. Exposed so an *edited* log can be
/// replayed (H5 provenance A/B): inserting an op shifts every later index, so the
/// raw handles stored in `AddSegment` then resolve to the **wrong** points.
pub fn replay_ops(log: &[Op]) -> Sketch {
    let mut m = Sketch::new();
    for &op in log {
        match op {
            Op::AddPoint { x, y } => {
                m.add_point(x, y);
            }
            Op::AddSegment { a, b } => {
                m.add_segment(a, b);
            }
        }
    }
    m
}

/// A stable element id — assigned at creation, independent of position in the
/// store or the log. This is the minimal "semantic" reference: identity by name,
/// not by index. (The full lineage — "the edge op X generated from faces …" —
/// is the richer form; a stable id already survives insertion.)
pub type PointId = u64;

/// Operations that reference points by **stable id** instead of raw handle.
#[derive(Clone, Copy, Debug)]
pub enum SemOp {
    AddPoint { id: PointId, x: Rat, y: Rat },
    AddSegment { a: PointId, b: PointId },
}

/// A document whose ops reference by stable id.
pub struct SemDoc {
    pub log: Vec<SemOp>,
}

impl SemDoc {
    /// Replay: build the model and resolve id references through an id→handle
    /// map. A segment whose endpoint id is missing (its creator was removed —
    /// a topology change) is **detected** and returned in `broken`, rather than
    /// silently connecting the wrong point (§7 TNP: detect, don't guess).
    pub fn replay(&self) -> (Sketch, Vec<(PointId, PointId)>) {
        use std::collections::HashMap;
        let mut m = Sketch::new();
        let mut id_to_handle: HashMap<PointId, Handle<Point2>> = HashMap::new();
        let mut broken = Vec::new();
        for &op in &self.log {
            match op {
                SemOp::AddPoint { id, x, y } => {
                    let h = m.add_point(x, y);
                    id_to_handle.insert(id, h);
                }
                SemOp::AddSegment { a, b } => match (id_to_handle.get(&a), id_to_handle.get(&b)) {
                    (Some(&ha), Some(&hb)) => {
                        m.add_segment(ha, hb);
                    }
                    _ => broken.push((a, b)),
                },
            }
        }
        (m, broken)
    }
}

impl Default for Document {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(n: i128) -> Rat {
        Rat::from_int(n)
    }

    fn build_square(d: &mut Document) {
        let c = [
            d.add_point(r(0), r(0)),
            d.add_point(r(1), r(0)),
            d.add_point(r(1), r(1)),
            d.add_point(r(0), r(1)),
        ];
        for i in 0..4 {
            d.add_segment(c[i], c[(i + 1) % 4]);
        }
    }

    fn sketches_eq(a: &Sketch, b: &Sketch) -> bool {
        a.points.len() == b.points.len()
            && a.segments.len() == b.segments.len()
            && a.points.iter().zip(b.points.iter()).all(|(x, y)| x == y)
            && a.segments
                .iter()
                .zip(b.segments.iter())
                .all(|(x, y)| x == y)
    }

    /// H5: the log is the truth; replaying it is deterministic — two replays give
    /// identical models (down to handle indices), and both match the live model
    /// built incrementally by `apply`.
    #[test]
    fn replay_is_deterministic_and_matches_live() {
        let mut d = Document::new();
        build_square(&mut d);
        assert_eq!(d.log().len(), 8); // 4 points + 4 segments

        let r1 = d.replay();
        let r2 = d.replay();
        assert!(sketches_eq(&r1, &r2), "two replays must be identical");
        assert!(
            sketches_eq(d.model(), &r1),
            "replay must match the live model"
        );
        assert_eq!(r1.points.len(), 4);
        assert_eq!(r1.segments.len(), 4);
    }

    /// The segment ops carry raw handles that resolve identically on replay: the
    /// replayed segments reference the same point indices as the live model.
    #[test]
    fn raw_handle_refs_are_replay_stable() {
        let mut d = Document::new();
        build_square(&mut d);
        let replay = d.replay();
        for (live, re) in d.model().segments.iter().zip(replay.segments.iter()) {
            assert_eq!(live.a, re.a);
            assert_eq!(live.b, re.b);
        }
    }

    // ---- H5 provenance A/B: raw handle vs semantic (stable id) under an edit ----

    fn build_square_sem() -> SemDoc {
        let mut log = Vec::new();
        for (id, (x, y)) in [(0u64, (0, 0)), (1, (1, 0)), (2, (1, 1)), (3, (0, 1))] {
            log.push(SemOp::AddPoint {
                id,
                x: r(x),
                y: r(y),
            });
        }
        for i in 0..4u64 {
            log.push(SemOp::AddSegment {
                a: i,
                b: (i + 1) % 4,
            });
        }
        SemDoc { log }
    }

    /// Raw handle refs are positional: inserting an op shifts every later index,
    /// so an `AddSegment{Handle(0), Handle(1)}` silently connects the wrong point
    /// after the edit — the failure the semantic scheme must avoid.
    #[test]
    fn raw_ref_breaks_on_insert() {
        let mut d = Document::new();
        build_square(&mut d);
        // edit: insert a new point at the FRONT of the log → indices all shift.
        let mut edited: Vec<Op> = d.log().to_vec();
        edited.insert(0, Op::AddPoint { x: r(9), y: r(9) });
        let m = replay_ops(&edited);
        // the first segment stored raw Handle(0); index 0 is now the (9,9) point.
        let seg = *m.segments.iter().next().unwrap();
        assert_eq!(
            m.points.get(seg.a).xy,
            (r(9), r(9)),
            "raw ref silently resolves to the inserted point (broken)"
        );
    }

    /// Semantic (stable-id) refs survive the same insert: ids do not shift, so the
    /// segment still resolves to the intended points.
    #[test]
    fn semantic_ref_survives_insert() {
        let mut doc = build_square_sem();
        doc.log.insert(
            0,
            SemOp::AddPoint {
                id: 99,
                x: r(9),
                y: r(9),
            },
        );
        let (m, broken) = doc.replay();
        assert!(broken.is_empty(), "no broken refs after an insert");
        let seg = *m.segments.iter().next().unwrap();
        let (a, b) = (m.points.get(seg.a).xy, m.points.get(seg.b).xy);
        assert_eq!(
            (a, b),
            ((r(0), r(0)), (r(1), r(0))),
            "semantic ref resolves to the intended points"
        );
    }

    /// A topology change — removing the op that created a referenced point — is
    /// **detected** (the id no longer resolves), not silently mis-connected. This
    /// is what lets §7 ask the user instead of guessing.
    #[test]
    fn semantic_ref_detects_removed_target() {
        let mut doc = build_square_sem();
        doc.log
            .retain(|op| !matches!(op, SemOp::AddPoint { id: 1, .. }));
        let (_m, broken) = doc.replay();
        assert!(
            !broken.is_empty(),
            "removing a referenced point must be detected"
        );
        assert!(broken.iter().any(|&(a, b)| a == 1 || b == 1));
    }
}
