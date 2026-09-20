//! **Which top-level module names which** — the crate's own dependency graph, read off the source.
//!
//! This crate is 35k lines of product code in eighteen top-level modules, and the question
//! "can this be split" is really "is the graph acyclic". Nothing in `fmt`, `clippy` or the suite
//! answers it, so the graph was invisible until it was measured: ten module pairs point at each
//! other, but only two of those are *behaviour* — the rest name a type that lives on the wrong
//! side. The difference is the whole finding, and it is only visible if something counts.
//!
//! ★ **The parser is what gets calibrated, not the files.** An earlier hand-rolled version of this
//! measurement was wrong three times — it counted `std::ops::Deref` as a reference to this crate's
//! `ops` module, then over-corrected and counted almost nothing, then broke on a lookbehind. So
//! [`targets`] is exercised against **literal lines** in [`the_parser_reads_what_it_should`]. Those
//! fixtures are strings on purpose: a fixture that named a file would go red the moment a module
//! moves, and that red would mean nothing while looking exactly like a real one.
//!
//! Being textual it can be talked around (a path inside a string literal counts, a module reached
//! through a `use super::*` glob does not). It is a tripwire on the ordinary way to name another
//! module, not a proof.

use std::collections::{BTreeMap, BTreeSet};

/// The top-level modules of this crate: every folder and every `.rs` file directly under `src`,
/// minus the crate root and the test tree.
fn modules() -> BTreeSet<String> {
    let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut out = BTreeSet::new();
    for entry in std::fs::read_dir(&src).expect("src") {
        let path = entry.expect("entry").path();
        let name = if path.is_dir() {
            path.file_name()
                .expect("dir name")
                .to_string_lossy()
                .into_owned()
        } else if path.extension().is_some_and(|e| e == "rs") {
            path.file_stem()
                .expect("stem")
                .to_string_lossy()
                .into_owned()
        } else {
            continue;
        };
        if name != "lib" && name != "tests" {
            out.insert(name);
        }
    }
    out
}

/// Every **other** top-level module this one line names, with duplicates kept — one entry per
/// spelling, because the count is what tells a type reference from a call site.
///
/// A path is a run of `ident (:: ident)*`. Its head decides what is being named: `crate::x::…`
/// names `x`, a bare `x::…` names `x`, and `std`/`core`/`alloc`/`super`/`self` name nothing here
/// (`std::ops::Deref` is not this crate's `ops`, and `super::` is inside a module, not across).
fn targets(line: &str, owner: &str, mods: &BTreeSet<String>) -> Vec<String> {
    let b = line.as_bytes();
    let ident_start = |c: u8| c.is_ascii_alphabetic() || c == b'_';
    let ident = |c: u8| c.is_ascii_alphanumeric() || c == b'_';
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < b.len() {
        if !ident_start(b[i]) {
            i += 1;
            continue;
        }
        let mut segs: Vec<&str> = Vec::new();
        let mut saw_sep = false;
        let mut j = i;
        loop {
            let s = j;
            while j < b.len() && ident(b[j]) {
                j += 1;
            }
            if j == s {
                break;
            }
            segs.push(&line[s..j]);
            if j + 1 < b.len() && b[j] == b':' && b[j + 1] == b':' {
                saw_sep = true;
                j += 2;
                if j < b.len() && ident_start(b[j]) {
                    continue;
                }
            }
            break;
        }
        if saw_sep && !segs.is_empty() {
            let head = segs[0];
            let named = if head == "crate" && segs.len() > 1 {
                segs[1]
            } else {
                head
            };
            let foreign = matches!(head, "std" | "core" | "alloc" | "super" | "self");
            if !foreign && named != owner && mods.contains(named) {
                out.push(named.to_string());
            }
        }
        i = j.max(i + 1);
    }
    out
}

/// Product source files of one top-level module, test trees excluded — a module that became a
/// folder is still one module, and a scan that stopped at the top level would go quiet about it.
fn sources_of(module: &str) -> Vec<std::path::PathBuf> {
    let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let base = src.join(module);
    let mut files = Vec::new();
    if base.is_dir() {
        let mut dirs = vec![base];
        while let Some(dir) = dirs.pop() {
            for entry in std::fs::read_dir(&dir).expect("module folder") {
                let path = entry.expect("entry").path();
                if path.is_dir() {
                    if path.file_name().is_some_and(|f| f != "tests") {
                        dirs.push(path);
                    }
                } else if path.extension().is_some_and(|e| e == "rs") {
                    files.push(path);
                }
            }
        }
    } else {
        files.push(src.join(format!("{module}.rs")));
    }
    files.sort();
    files
}

/// `owner -> named -> count` over the whole crate, comments dropped.
fn edges() -> BTreeMap<String, BTreeMap<String, usize>> {
    let mods = modules();
    let mut out: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    for owner in &mods {
        for path in sources_of(owner) {
            let text = std::fs::read_to_string(&path).expect("read");
            for line in text.lines() {
                let t = line.trim_start();
                if t.starts_with("//") || t.starts_with('*') || t.starts_with("/*") {
                    continue;
                }
                for named in targets(line, owner, &mods) {
                    *out.entry(owner.clone())
                        .or_default()
                        .entry(named)
                        .or_default() += 1;
                }
            }
        }
    }
    out
}

/// ★★ **The calibration, and the reason it is spelled as strings.** Each case is a line that was
/// actually in this crate when the graph was first measured, but it is pinned here as text: the
/// subject under test is the parser, and a fixture that pointed at a file would go red when that
/// file moves — a red that means nothing and looks like a real one.
///
/// The fourth case is the one that matters. `std::ops::Deref` shares a name with this crate's
/// `ops` module, and counting it put three phantom edges into the first measurement.
#[test]
fn the_parser_reads_what_it_should() {
    let mods: BTreeSet<String> = ["arrangement", "boolean", "exact", "ops", "combinatorics"]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    let cases: [(&str, &str, &[&str]); 5] = [
        (
            "            let crate::boolean::Wall::Plane(c) = *w else {",
            "arrangement",
            &["boolean"],
        ),
        (
            "    cut_rims: &crate::arrangement::CutRims,",
            "boolean",
            &["arrangement"],
        ),
        (
            "use crate::ops::{Profile2d, SketchPlane};",
            "exact",
            &["ops"],
        ),
        ("impl std::ops::Deref for Ring {", "boolean", &[]),
        // A glob re-export names its module even though no second segment follows.
        (
            "pub use combinatorics::*;",
            "arrangement",
            &["combinatorics"],
        ),
    ];
    for (line, owner, want) in cases {
        let got = targets(line, owner, &mods);
        assert_eq!(got, want, "parsing {line:?} as {owner}");
    }
}

/// The graph as it stands, printed for a reader. No assertion yet — this is the census the
/// restructure is measured against, and the gate that freezes it comes once the edges it is
/// meant to forbid are actually gone.
#[test]
#[ignore = "census: prints the module graph, asserts nothing"]
fn measure_module_graph() {
    let mods = modules();
    let e = edges();
    println!("\n{} top-level modules\n", mods.len());
    for owner in &mods {
        let row = e.get(owner);
        let mut cells: Vec<String> = Vec::new();
        for named in mods.iter() {
            if let Some(n) = row.and_then(|r| r.get(named)) {
                cells.push(format!("{named} {n}"));
            }
        }
        println!(
            "{owner:<16} -> {}",
            if cells.is_empty() {
                "-".into()
            } else {
                cells.join(" · ")
            }
        );
    }
    let mut cycles = Vec::new();
    for a in &mods {
        for b in &mods {
            if a >= b {
                continue;
            }
            let ab = e.get(a).and_then(|r| r.get(b)).copied().unwrap_or(0);
            let ba = e.get(b).and_then(|r| r.get(a)).copied().unwrap_or(0);
            if ab > 0 && ba > 0 {
                cycles.push(format!("{a} ⇄ {b}  ({ab} / {ba})"));
            }
        }
    }
    println!("\n{} cyclic pairs", cycles.len());
    for c in &cycles {
        println!("  {c}");
    }
}
