use std::sync::Mutex;

/// One chart's shape, and how today's emitted lateral faces cover its cells.
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
    /// one-cut-one-cell branch. ★ A zero here means that branch is *unexercised*, not verified.
    pub(crate) odd_k: usize,
}

impl Row {
    pub(crate) fn empty() -> Self {
        Self::default()
    }
}

pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

/// **Capability D's third rung has its own ledger**, deliberately not more fields on [`Row`].
/// That one is already fourteen wide, and every field of it needs a reader or `dead_code`
/// stops the build — which is how a table nobody can read grows contrived invariants to feed
/// it. One instrument per rung.
pub(crate) mod d2 {
    use std::sync::Mutex;

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
        /// Intervals whose walk around the circle returns to where it started, and those whose
        /// does not — `label_cells`' final verification, stated on the chart.
        pub(crate) closes: usize,
        pub(crate) does_not_close: usize,
        /// Intervals where some wall is crossed an odd number of times (a panel's or a
        /// chain's wall with one ruling in the interval): the walk has an open end there and
        /// the closure is not asserted.
        pub(crate) unpaired: usize,
        /// Rulings whose every mark is a **graze** — the face stops at the line rather than
        /// crossing it, so whether it reaches the interval is `face_spans`' question and not a
        /// label's. Existence, not membership, on the vertical axis.
        pub(crate) grazing_rulings: usize,
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    pub(crate) fn push(r: Row) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
}

pub(crate) fn push(r: Row) {
    ROWS.lock()
        .expect("the probe's lock is never held across a panic")
        .push(r);
}

/// **The fourth rung's ledger**: what the cell reader read, what the
/// census's own walk of the chart predicts, and how many faces the emitter put on the class.
/// One row per chart with cells, whether or not the emitter produced a face there.
///
/// ★ Compared **within a row** (`cells` is copied in), never row-by-row against another
/// ledger: tests run in parallel and the ledgers interleave independently.
/// **The vertical answer beside the horizontal one** — the shadow the cutover is measured
/// against (capability D's own pattern: a rung that ships no capability and *measures*).
///
/// The chart's two axes are symmetric in principle, but only the horizontal ones have ever
/// decided a face, so the sign bridge on the vertical side ([`crate::arrangement::
/// plus_theta_is_above`]) is unexercised. These count what it would say.
pub(crate) mod shadow {
    use std::sync::Mutex;

    /// `(agree, disagree, vertical_only, horizontal_only, both_silent, one_ruling)`
    pub(crate) static COUNTS: Mutex<(usize, usize, usize, usize, usize, usize)> =
        Mutex::new((0, 0, 0, 0, 0, 0));

    pub(crate) fn record(
        horizontal: Option<(bool, bool)>,
        vertical: &[(bool, bool)],
        one_ruling: bool,
    ) {
        let mut g = COUNTS
            .lock()
            .expect("the probe's lock is never held across a panic");
        g.5 += usize::from(one_ruling);
        // The vertical readings must agree with each other before they may agree with anyone.
        let v = match vertical {
            [] => None,
            [a] => Some(*a),
            [a, b] if a == b => Some(*a),
            _ => {
                g.1 += 1; // two walls of one cell disagreeing is a disagreement of its own
                return;
            }
        };
        match (horizontal, v) {
            (Some(h), Some(x)) if h == x => g.0 += 1,
            (Some(_), Some(_)) => g.1 += 1,
            (None, Some(_)) => g.2 += 1,
            (Some(_), None) => g.3 += 1,
            (None, None) => g.4 += 1,
        }
    }
}

/// The counts behind [`super::OtherWhy`], and the whole-circle disagreements in full: which
/// test, which line, and the per-arc bits that disagreed.
pub(crate) mod other {
    use super::super::OtherWhy;
    use std::sync::Mutex;

    pub(crate) static COUNTS: Mutex<Vec<(OtherWhy, usize)>> = Mutex::new(Vec::new());

    #[derive(Clone, Debug)]
    pub(crate) struct Whole {
        pub(crate) test: String,
        pub(crate) cyl: usize,
        pub(crate) t: f64,
        pub(crate) end: usize,
        pub(crate) above: bool,
        pub(crate) bits: Vec<(bool, bool)>,
    }

    pub(crate) static WHOLE: Mutex<Vec<Whole>> = Mutex::new(Vec::new());

    pub(crate) fn record(why: OtherWhy) {
        let mut c = COUNTS
            .lock()
            .expect("the probe's lock is never held across a panic");
        match c.iter_mut().find(|(w, _)| *w == why) {
            Some((_, n)) => *n += 1,
            None => c.push((why, 1)),
        }
    }

    pub(crate) fn record_whole(
        cyl: usize,
        t: nacre_scalar::Rat,
        end: usize,
        above: bool,
        bits: Vec<(bool, bool)>,
    ) {
        record(OtherWhy::WholeDisagree);
        WHOLE
            .lock()
            .expect("the probe's lock is never held across a panic")
            .push(Whole {
                test: std::thread::current().name().unwrap_or("?").to_string(),
                cyl,
                t: t.to_f64(),
                end,
                above,
                bits,
            });
    }
}

/// **The emitter's regions**, one row per class: the faces the walk made and the cells
/// they cover (the census asserts the cell→face assignment where it is made).
pub(crate) mod regions {
    use std::sync::Mutex;

    #[derive(Clone, Debug)]
    pub(crate) struct Row {
        pub(crate) test: String,
        pub(crate) cyl: usize,
        pub(crate) emitter_refused: bool,
        pub(crate) faces: usize,
        pub(crate) band_faces: usize,
        pub(crate) ring_faces: usize,
        pub(crate) emitted_cells: usize,
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    pub(crate) fn push(r: Row) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
}

pub(crate) mod d2b {
    use std::sync::Mutex;

    #[derive(Clone, Debug, Default)]
    pub(crate) struct Row {
        /// Every (z-line, station) pair asked for the station's canonical
        /// name on that line (`crossing_on_ruling`), and how many could not be named.
        pub(crate) station_pairs: usize,
        pub(crate) station_name_failures: usize,
        /// Of `end_other`, the ends of a single-cut sector (both walls one ruling).
        pub(crate) end_other_single_cut: usize,
        pub(crate) cells: usize,
        /// The emitter refused this boolean's lateral faces (some class's cells could not be
        /// read) — the census still records every chart of it.
        pub(crate) emitter_refused: bool,
        /// Rulings whose `end` had to be swapped alongside `z` — predicted 0.
        pub(crate) end_swapped: usize,
        /// Every cell's two ends, by what they said.
        pub(crate) end_disk: usize,
        pub(crate) end_exact: usize,
        pub(crate) end_other: usize,
        pub(crate) end_nocircle: usize,
        /// Ends between adjacent rim nodes no arc covers — the face is absent there.
        pub(crate) end_uncovered: usize,
        /// Present cells with an `Other` end — read from their other end alone. The
        /// population a θ-placement finer than [`Chart::arc_around`] would serve.
        pub(crate) other_present: usize,
        /// Present cells with a `NoCircle` end — the `circle_on_class` premise at one end.
        pub(crate) nocircle_present: usize,
        /// Cells whose two speaking ends disagreed — the two-end refusal
        /// (`CylinderGateUndecided`), on the chart; the emitter refuses the class.
        pub(crate) src2_disagree: usize,
        /// Cells with a face and no speaking end. Asserted 0 at the record; here as a count
        /// so the reporting test can say the assertion was live.
        pub(crate) src0_present: usize,
        pub(crate) exist_disagree: usize,
        /// The reader refused a cell of this chart by name (`face_spans`' refusal, the
        /// two-cut-ends refusal, or an arrangement-table inconsistency) — the whole class is
        /// refused.
        pub(crate) read_refused: usize,
        /// Both-cut cells the trace says the face is not in — the panel road's dropped sectors.
        pub(crate) exist_marks_false: usize,
        /// Cells the reader would emit, and cells it could not decide.
        pub(crate) emit: usize,
        pub(crate) emit_unknown: usize,
        /// Cut ends read, and those carrying no / more than one lateral mark of their own solid
        /// (the one-mark contract).
        pub(crate) arcs_read: usize,
        pub(crate) arcs_no_mark: usize,
        pub(crate) arcs_multi_mark: usize,
        /// The emitter's lateral faces of this class (0 when it refused).
        /// Σ(len − 1) over `End::Exact` runs: the arcs an end reads **beyond its first**.
        /// A sector that spans a rim node the chart has no line for reads the whole run.
        pub(crate) exact_run_arcs: usize,
        pub(crate) emitted_faces: usize,
    }

    pub(crate) static ROWS: Mutex<Vec<Row>> = Mutex::new(Vec::new());

    pub(crate) fn push(r: Row) {
        ROWS.lock()
            .expect("the probe's lock is never held across a panic")
            .push(r);
    }
}
