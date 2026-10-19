use super::*;

/// **What carries a ring edge, as reading a planar face's ring against a class line needs it** —
/// one view of the two ring forms the engine holds: the tracer's [`Wall`]s and a [`RingEdge`]'s
/// [`Carrier`]. A carrier cannot become a wall (an [`ArcCarrier`] keeps no plane), so the view is
/// its own type and each form converts into it here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LineCarrier {
    /// A straight edge on the face's plane and this wall class.
    Plane(usize),
    /// An arc of cylinder `cyl`, travelled counter-clockwise about its axis when `ccw`.
    Arc { cyl: usize, ccw: bool },
    /// A ruling of cylinder `cyl`: `side` names which of `plane`'s two, measured against `plane`.
    Ruling { cyl: usize, side: i8, plane: usize },
}

impl From<Wall> for LineCarrier {
    #[inline]
    fn from(w: Wall) -> Self {
        match w {
            Wall::Plane(p) => LineCarrier::Plane(p),
            Wall::Arc { cyl, ccw, .. } => LineCarrier::Arc { cyl, ccw },
            Wall::Ruling {
                cyl, side, plane, ..
            } => LineCarrier::Ruling { cyl, side, plane },
        }
    }
}

impl From<&Carrier> for LineCarrier {
    #[inline]
    fn from(c: &Carrier) -> Self {
        match c {
            Carrier::Plane { wall, .. } => LineCarrier::Plane(*wall),
            Carrier::Arc(a) => LineCarrier::Arc {
                cyl: a.cyl,
                ccw: a.ccw,
            },
            Carrier::Ruling(r) => LineCarrier::Ruling {
                cyl: r.cyl,
                side: r.side,
                plane: r.plane,
            },
        }
    }
}

/// **A planar face's ring read against the class line `face ∩ cut`** — the one rule for where
/// that line meets the ring and what it meets there, read by every road that walks a plane face's
/// ring along a line: the tracer (`trace_transversal_face`), the ray caster ([`every_ray`]) and the
/// self-touch chord test ([`segment_meets_face`]). Each brings its own cylinder table — the whole
/// one, or the narrower one a road reads by design — and what it does with a feature stays its
/// own: the order it handles runs and crossings in, the names it records, the words it refuses in.
///
/// The meet is a **line**, so a straight edge between two nodes on it lies on it, and an arc leaves
/// it to the side its tangent points ([`arc_departure_side`]). An arc whose ends cannot say how
/// often it crosses is asked ([`arc_crossings_of`]), and the same list names the crossings the walk
/// counted. A crossing is named by what carries its edge: a plane wall pins it on the line; on an
/// arc it is the `nth` root the arc strictly contains, on a ruling the root that ruling is
/// ([`crossing_on_ruling`]) — pierce points of the line's own pair either way.
///
/// **Read lazily.** A road walks the ring first and names a crossing only when it reaches it,
/// after the runs it handles before (a ray skips itself on an unnameable run end; the tracer records
/// a fourth plane through a run's node). Naming every crossing up front would change which refusal
/// comes first.
///
/// **An arc's crossings are solved once per edge**, so the count the walk read and the point a
/// road names cannot disagree. The walk asks an edge only where its ends cannot say, with the end
/// on the line where there is one (that end's own root is not counted); an edge the walk did not
/// ask (ends on opposite sides) is solved on naming, with no end on the line. Nothing is allocated
/// until an edge is asked — most rings a ray walks carry no arc, and a ray builds one reading per
/// line it casts.
pub(crate) struct ClassLine<'a, C> {
    jd: &'a Judge<'a, WorkingPlane>,
    cyls: &'a [crate::planes::WorkingCyl],
    /// The face's own plane class.
    face: usize,
    /// The class that cuts it along the line. Named apart from `face` because the primitives
    /// take the pair as `(cut, face)` and a swap would be silent.
    cut: usize,
    nodes: &'a [NodeId],
    carrier: C,
    arcs: std::cell::RefCell<Vec<Option<Option<Vec<NodeId>>>>>,
}

impl<'a, C: Fn(usize) -> LineCarrier> ClassLine<'a, C> {
    /// The ring `nodes` (edge `i` from node `i` to node `i + 1`, carried by `carrier(i)`) on plane
    /// `face`, against the line `face ∩ cut`, with the cylinder table `cyls`.
    pub(crate) fn new(
        jd: &'a Judge<'a, WorkingPlane>,
        cyls: &'a [crate::planes::WorkingCyl],
        face: usize,
        cut: usize,
        nodes: &'a [NodeId],
        carrier: C,
    ) -> Self {
        ClassLine {
            jd,
            cyls,
            face,
            cut,
            nodes,
            carrier,
            arcs: std::cell::RefCell::new(Vec::new()),
        }
    }

    /// Where the ring meets the line — [`ring_against_plane`] with this reading's two answers.
    #[inline]
    pub(crate) fn walk(&self) -> RingWalk {
        ring_against_plane(
            self.jd,
            self.cyls,
            self.nodes,
            self.cut,
            |i| match (self.carrier)(i) {
                LineCarrier::Arc { cyl, ccw } => {
                    arc_departure_side(self.jd, self.cyls, self.nodes[i], self.cut, cyl, ccw)
                        .map(EdgeMeet::Departs)
                }
                LineCarrier::Plane(_) | LineCarrier::Ruling { .. } => Some(EdgeMeet::On),
            },
            |i, on_end| {
                self.arc_roots(i, on_end)
                    .and_then(|roots| u8::try_from(roots.len()).ok())
            },
        )
    }

    /// The point where the walk's `Crossing { edge, nth }` lies on the line; `None` where it has no
    /// exact name — a ruling measured against another plane than the face's, a cylinder missing
    /// from the table, a root the arc's list does not hold.
    #[inline]
    pub(crate) fn crossing(&self, edge: usize, nth: u8) -> Option<PointOn> {
        match (self.carrier)(edge) {
            LineCarrier::Plane(wall) => Some(PointOn::Class(wall)),
            LineCarrier::Arc { .. } => self
                .arc_roots(edge, None)?
                .get(usize::from(nth))
                .copied()
                .map(PointOn::Pierce),
            // `side` reads only against the plane it was measured on — on a plane face's ring,
            // the face's own.
            LineCarrier::Ruling { cyl, side, plane } if plane == self.face => crossing_on_ruling(
                self.jd,
                &self.cyls.get(cyl)?.def,
                self.cut,
                self.face,
                cyl,
                side,
            )
            .ok()
            .map(PointOn::Pierce),
            LineCarrier::Ruling { .. } => None,
        }
    }

    /// Where edge `i` crosses the line strictly inside it — the first ask's answer, whatever end a
    /// later ask passes.
    #[inline]
    fn arc_roots(&self, i: usize, on_end: Option<usize>) -> Option<Vec<NodeId>> {
        if let Some(Some(got)) = self.arcs.borrow().get(i) {
            return got.clone();
        }
        let n = self.nodes.len();
        let got = match (self.carrier)(i) {
            LineCarrier::Arc { cyl, ccw } => arc_crossings_of(
                self.jd,
                self.cyls,
                self.cut,
                self.face,
                cyl,
                ccw,
                [self.nodes[i], self.nodes[(i + 1) % n]],
                on_end.map(|k| self.nodes[k]),
            ),
            LineCarrier::Plane(_) | LineCarrier::Ruling { .. } => Some(Vec::new()),
        };
        let mut memo = self.arcs.borrow_mut();
        if memo.is_empty() {
            memo.resize(n, None);
        }
        memo[i] = Some(got.clone());
        got
    }
}
