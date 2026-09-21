//! **Every silenced `dead_code` warning says why it is silenced.**
//!
//! `#[cfg_attr(not(test), allow(dead_code))]` is this crate's idiom for an item a *product* build
//! legitimately never reads — an instrument's field, a variant only the census builds
//! (`arrangement::result`'s note argues for it over a blanket allow). It is also the exact
//! spelling that hid a defect: `crossing_probe` and `cycle_probe` were whole instrumentation
//! modules with no `#[cfg(test)]` at their mounts, so a release library carried them, and this
//! attribute silenced the one warning a product build raised about it.
//!
//! ★ **The attribute is not the tell; the sentence above it is.** Of the nine sites that existed
//! when this was written, seven said why in the line above ("Read by the census only…", "Used by
//! `loop_winding`'s tests…") and the two that said nothing were the two defects. The population
//! separated perfectly, so that is what this pins: silencing the compiler is allowed, silencing it
//! without a word is not.
//!
//! ★★ **The window is four lines, and it was calibrated against the real population.** At three or
//! more it flags exactly the two defects; at two it also flags `cyl_chart::regions`' second field,
//! which shares the comment written above its sibling. Four leaves a line of margin over the
//! narrowest width that works.
//!
//! Being textual it can be talked around (any `//` line in the window satisfies it, including an
//! unrelated one). It is a tripwire on the ordinary way to get this wrong, not a proof — so
//! [`the_window_reads_what_it_should`] pins the rule against literal lines rather than against
//! files, which would go red whenever a file moves and mean nothing when it did.

/// The spelling that silences a product build's `dead_code` warning.
const SILENCER: &str = "cfg_attr(not(test), allow(dead_code))";

/// How far above the attribute a justifying comment may sit.
const WINDOW: usize = 4;

/// Does a comment stand within [`WINDOW`] lines above `lines[at]`?
fn justified(lines: &[&str], at: usize) -> bool {
    lines[at.saturating_sub(WINDOW)..at]
        .iter()
        .any(|l| l.trim_start().starts_with("//"))
}

fn sources() -> Vec<std::path::PathBuf> {
    let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut out = Vec::new();
    let mut dirs = vec![src];
    while let Some(dir) = dirs.pop() {
        for entry in std::fs::read_dir(&dir).expect("src") {
            let path = entry.expect("entry").path();
            if path.is_dir() {
                dirs.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

/// ★ No exemptions, and that is the point: the list is empty because the two that could not have
/// been written down were the defect, not because every awkward case was excused.
#[test]
fn a_silenced_warning_states_its_reason() {
    let mut offenders = Vec::new();
    for path in sources() {
        let text = std::fs::read_to_string(&path).expect("read");
        let lines: Vec<&str> = text.lines().collect();
        for (n, line) in lines.iter().enumerate() {
            if line.contains(SILENCER) && !justified(&lines, n) {
                offenders.push(format!("{}:{}", path.display(), n + 1));
            }
        }
    }
    assert_eq!(
        offenders,
        Vec::<String>::new(),
        "a product build's `dead_code` warning is silenced with no sentence saying why -- write \
         one, or ask whether the item belongs behind `#[cfg(test)]` instead"
    );
}

/// ★ **The rule is exercised on strings, not on the tree.** The whole workspace's population of
/// this attribute lives in this crate, and it is small; a scan that only ever sees a clean tree
/// has not been shown to say anything at all.
#[test]
fn the_window_reads_what_it_should() {
    let bare = ["", "#[derive(Clone, Copy, Debug)]", "#[cfg_attr(..)]"];
    assert!(
        !justified(&bare, 2),
        "a derive and a blank line are not a reason"
    );

    let doc = ["/// Read by the census only.", "#[cfg_attr(..)]"];
    assert!(
        justified(&doc, 1),
        "a doc comment directly above is a reason"
    );

    let line_comment = ["// production reads `chamber`.", "", "#[cfg_attr(..)]"];
    assert!(
        justified(&line_comment, 2),
        "a `//` comment counts, blank lines do not break it"
    );

    // The shape the window exists for: a run of fields sharing the sentence above the first.
    let run = [
        "/// Read by the census only -- an instrument's fields.",
        "#[cfg_attr(..)]",
        "pub(crate) cells: Vec<usize>,",
        "#[cfg_attr(..)]",
    ];
    assert!(
        justified(&run, 3),
        "a sibling field shares the sentence above the run"
    );

    // At the top of a file there is nothing above, and `saturating_sub` must not panic.
    assert!(
        !justified(&["#[cfg_attr(..)]"], 0),
        "the first line has no reason above it"
    );
}
