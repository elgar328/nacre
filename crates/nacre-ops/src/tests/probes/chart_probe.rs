//! **An instrument, so the file lives in the test tree and the module does not.**
//! [`crate::arrangement::cyl_chart`] mounts it with `#[path]` as `probe`, which is what keeps `super::`
//! here meaning that module — moving this file did not move what it belongs to.

use crate::ledger::Ledger;

/// One chart's shape.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Row {
    pub(crate) z_lines: usize,
    pub(crate) theta: usize,
    /// How many **rows** (lateral faces) the class carries.
    pub(crate) rows: usize,
    /// The θ order could not be formed, so this chart has **no** cells — recorded rather than
    /// silently dropped, because a partial chart compared as a whole one is how a counting
    /// claim goes blind.
    pub(crate) refused: bool,
    pub(crate) cells: usize,
    /// Intervals with no ruling over them: the cell is the whole circle.
    pub(crate) whole_circle: usize,
    /// Intervals crossed by an **odd** number of rulings — the only way to reach the
    /// one-cut-one-cell branch. A zero here means that branch is *unexercised*, not verified.
    pub(crate) odd_k: usize,
}

impl Row {
    pub(crate) fn empty() -> Self {
        Self::default()
    }
}

pub(crate) static ROWS: Ledger<Row> = Ledger::new();

/// **The chart's vertical answers have their own ledger**, deliberately not more fields on
/// [`Row`]: every field needs a reader or `dead_code` stops the build, which is how a table nobody
/// can read grows contrived invariants to feed it. One instrument per question.
pub(crate) mod rulings {
    use crate::ledger::Ledger;

    /// One chart's vertical answers.
    #[derive(Clone, Copy, Debug, Default)]
    pub(crate) struct Row {
        /// Alive-ruling observations across all of this chart's intervals — the denominator
        /// the two counts below are subsets of. Kept so a ratio is never read against a
        /// denominator that lives in another table.
        pub(crate) rulings: usize,
        /// Rulings whose wall changes the material at the lateral.
        pub(crate) wall_flips: usize,
        pub(crate) intervals_with_flip: usize,
        /// Intervals whose walk around the circle returns to where it started — `label_cells`'
        /// final verification, stated on the chart and asserted for every paired interval.
        pub(crate) closes: usize,
        /// Intervals where some wall is crossed an odd number of times (a panel's or a
        /// chain's wall with one ruling in the interval): the walk has an open end there and
        /// the closure is not asserted.
        pub(crate) unpaired: usize,
        /// Rulings whose every mark is a **graze** — the face stops at the line rather than
        /// crossing it, so whether it reaches the interval is `face_spans`' question and not a
        /// label's. Existence, not membership, on the vertical axis.
        pub(crate) grazing_rulings: usize,
    }

    pub(crate) static ROWS: Ledger<Row> = Ledger::new();

    pub(crate) fn push(r: Row) {
        ROWS.push(r);
    }
}

pub(crate) fn push(r: Row) {
    ROWS.push(r);
}

/// **The cells' two ends**: what the cell reader read, what the census's own walk of the chart
/// predicts, and how many faces the emitter put on the class. One row per chart with cells,
/// whether or not the emitter produced a face there.
///
/// Compared **within a row** (`cells` is copied in), never row-by-row against another
/// ledger: tests run in parallel and the ledgers interleave independently.
pub(crate) mod cell_ends {
    use crate::ledger::Ledger;

    #[derive(Clone, Debug, Default)]
    pub(crate) struct Row {
        pub(crate) cells: usize,
        /// Every cell's two ends, by what they said (an `Other` end is not counted — the census
        /// asserts no present cell has one).
        pub(crate) end_disk: usize,
        pub(crate) end_exact: usize,
        pub(crate) end_nocircle: usize,
        /// Ends between adjacent rim nodes no arc covers — the face is absent there.
        pub(crate) end_uncovered: usize,
        /// Cells whose trace and span disagree about the face being there.
        pub(crate) exist_disagree: usize,
        /// The reader refused a cell of this chart by name (`face_spans`' refusal, the
        /// two-cut-ends refusal, or an arrangement-table inconsistency) — the whole class is
        /// refused.
        pub(crate) read_refused: usize,
        /// Both-cut cells the trace says the face is not in — sectors dropped for existence.
        pub(crate) exist_marks_false: usize,
        /// Cells the reader would emit.
        pub(crate) emit: usize,
        /// Cut ends read — one per arc of an `End::Exact` run.
        pub(crate) arcs_read: usize,
        /// Σ(len − 1) over `End::Exact` runs: the arcs an end reads **beyond its first**.
        /// A sector that spans a rim node the chart has no line for reads the whole run.
        pub(crate) exact_run_arcs: usize,
        /// The emitter's lateral faces of this class (0 when it refused).
        pub(crate) emitted_faces: usize,
    }

    pub(crate) static ROWS: Ledger<Row> = Ledger::new();

    pub(crate) fn push(r: Row) {
        ROWS.push(r);
    }
}
