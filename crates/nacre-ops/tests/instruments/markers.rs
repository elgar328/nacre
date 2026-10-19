//! **No importance marks in the repository's prose — the sentence carries the weight.**
//!
//! overview's document rule 8: a risk is stated as a risk and a premise as a premise, not flagged
//! with a mark. The marks are the black star (U+2605), the white star (U+2606), the warning sign
//! (U+26A0) and the ballot box with check (U+2611). A check mark (U+2713, U+2714) is
//! content where it appears — the README's milestone legend, a cell of the table the deadname sweep
//! prints — and is not on the list. This file spells the marks as escapes, so it does not trip
//! itself.
//!
//! What it reads: every file under `crates/`, `docs/` and `tools/`, and the root's prose files by
//! name ([`ROOT_FILES`]). A name starting with a dot, and `__pycache__`, are skipped — what lives
//! there is local (a `.DS_Store`, a bytecode cache). Not the whole repository and not every root
//! file: what exists only on one machine (`target/`, an index folder, an ignored `CLAUDE.md`) would
//! make this red there alone. So a new prose file at the root has to be added to [`ROOT_FILES`]. A
//! file that does not read as UTF-8 fails the check instead of being skipped: a scan that quietly
//! reads less is how this one would fail.
//!
//! Being textual it is a tripwire, not a proof — so [`the_scan_reads_what_it_should`] pins the rule
//! against literal lines rather than against files.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The importance marks: black star, white star, warning sign, ballot box with check.
const MARKS: [char; 4] = ['\u{2605}', '\u{2606}', '\u{26A0}', '\u{2611}'];

/// The root's prose files, read by name (relative to the repository root).
const ROOT_FILES: [&str; 5] = [
    "Cargo.toml",
    "clippy.toml",
    "README.md",
    "CHANGELOG.md",
    ".cargo-husky/hooks/pre-commit",
];

/// Lines that may carry a mark, per file (relative to the repository root), with why.
const ALLOWED: &[(&str, usize, &str)] =
    &[("docs/overview.md", 1, "the rule names the marks it forbids")];

fn marked(line: &str) -> bool {
    line.chars().any(|c| MARKS.contains(&c))
}

/// Every file the check reads, as `(path relative to the repository root, path)`, sorted.
fn prose_files() -> Vec<(String, PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("the repository root, two above this crate");
    // Asked twice, of two kinds of name: the walk asks it of folders it might enter, this loop of
    // the files it hands back — a `.DS_Store` is a file, so the folder rule never sees it.
    let local = |name: &str| name.starts_with('.') || name == "__pycache__";
    let mut paths = Vec::new();
    for dir in ["crates", "docs", "tools"] {
        for path in crate::walk::files_under(&root.join(dir), |name| !local(name)) {
            if !path.file_name().and_then(|n| n.to_str()).is_some_and(local) {
                paths.push(path);
            }
        }
    }
    paths.extend(ROOT_FILES.iter().map(|f| root.join(f)));
    let mut out: Vec<(String, PathBuf)> = paths
        .into_iter()
        .map(|path| {
            let rel = path
                .strip_prefix(root)
                .expect("under the repository root")
                .to_string_lossy()
                .replace('\\', "/");
            (rel, path)
        })
        .collect();
    out.sort();
    out
}

#[test]
fn no_prose_carries_an_importance_mark() {
    let mut found: BTreeMap<String, usize> = BTreeMap::new();
    for (rel, path) in prose_files() {
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("{rel} does not read as UTF-8 text: {e}"));
        let n = text.lines().filter(|l| marked(l)).count();
        if n > 0 {
            found.insert(rel, n);
        }
    }
    let allowed: BTreeMap<String, usize> = ALLOWED
        .iter()
        .map(|(f, n, _)| (f.to_string(), *n))
        .collect();
    assert_eq!(
        found, allowed,
        "an importance mark outside the list -- say the weight in the sentence (a risk as a risk, \
         a premise as a premise) and drop the mark; the list holds only lines that must name a \
         mark (and drop an entry whose mark is gone)"
    );
}

/// **The rule is exercised on strings, not on the tree** — a scan that only ever sees a tree it
/// agrees with has not been shown to say anything. The lines are built from escapes, like
/// [`MARKS`], so this file holds no mark of its own.
#[test]
fn the_scan_reads_what_it_should() {
    for mark in MARKS {
        assert!(
            marked(&format!("/// {mark} **A claim.**")),
            "{mark:?} at a comment's head"
        );
        assert!(
            marked(&format!("text. {mark}{mark} More text")),
            "{mark:?} mid-line"
        );
    }
    assert!(
        marked("> \u{26A0}\u{FE0F} **Experimental**"),
        "a warning sign with its emoji selector"
    );
    assert!(!marked("/// **A claim.**"), "bold is not a mark");
    assert!(
        !marked("`M1 \u{2713}   M2 \u{2713}   M6 \u{25D4}   M7 \u{25CB}`"),
        "the README's milestone legend is content"
    );
    assert!(
        !marked("| 12 | `name` | x | \u{2714} | |"),
        "a check in a table cell is content"
    );
}
