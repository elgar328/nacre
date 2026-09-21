//! **Which top-level module names which** — the crate's own dependency graph, read off the source.
//!
//! This crate is 35k lines of product code in eighteen top-level modules, and the question
//! "can this be split" is really "is the graph acyclic". Nothing in `fmt`, `clippy` or the suite
//! answers it, so the graph was invisible until it was measured: nine module pairs pointed at
//! each other, and only two of those were *behaviour* — the rest named a type, or took an
//! argument, that lived on the wrong side. The difference is the whole finding, and it is only
//! visible if something counts. Two pairs are left.
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

/// ★★ **The calibration, and why its module names are made up.**
///
/// The subject under test is [`targets`], and it does not know this crate: it takes the module
/// set as an argument. So the fixtures name `alpha`/`beta`/`gamma`, which exist nowhere. That is
/// not tidiness — it is the second thing that went wrong here.
///
/// The first was pinning fixtures to *files*: a fixture that points at `assembly/naming.rs` goes
/// red the moment that file moves, and that red means nothing while looking exactly like a real
/// one. So they became string literals. Then a bulk path rewrite (`crate::assembly::Wall` ->
/// `crate::draft::Wall`) swept every `.rs` file in the crate and **edited the literals**, because
/// a fixture that looks like code is code to a script. Names that no module has are immune to
/// both.
///
/// The `std::ops::Deref` case keeps a real name on purpose: `ops` is a module of this crate, and
/// reading that line as a reference to it put three phantom edges into the first measurement.
#[test]
fn the_parser_reads_what_it_should() {
    let mods: BTreeSet<String> = ["alpha", "beta", "gamma", "ops"]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    let cases: [(&str, &str, &[&str]); 6] = [
        // A pattern match through a full path: the head is `crate`, so the module is the second
        // segment, and the type and variant after it are not modules.
        (
            "            let crate::beta::Wall::Plane(c) = *w else {",
            "alpha",
            &["beta"],
        ),
        // A borrowed type in a signature.
        ("    cut_rims: &crate::alpha::CutRims,", "beta", &["alpha"]),
        // A brace import names its module even though the braces stop the path.
        (
            "use crate::gamma::{Profile2d, SketchPlane};",
            "alpha",
            &["gamma"],
        ),
        // ★ A foreign root that shares a name with a module of this crate. This is the case the
        // whole calibration exists for.
        ("impl std::ops::Deref for Ring {", "beta", &[]),
        // A glob re-export names its module although no second segment follows.
        ("pub use gamma::*;", "alpha", &["gamma"]),
        // A module never names itself, and one line can name two.
        (
            "    let x = alpha::f(beta::g(), gamma::h());",
            "alpha",
            &["beta", "gamma"],
        ),
    ];
    for (line, owner, want) in cases {
        let got = targets(line, owner, &mods);
        assert_eq!(got, want, "parsing {line:?} as {owner}");
    }
}

/// ★★★ **What this restructure established, asserted — and only that.**
///
/// The graph is **not** a DAG and this does not pretend otherwise: two module pairs still
/// point at each other. What is settled is the boolean pipeline's shape and the direction of
/// everything under the engine, so that is what is locked. Each of these was a real edge before
/// the work and is zero after it; an editor who reintroduces one is undoing something, not
/// adding to it.
#[test]
fn the_pipeline_runs_one_way() {
    let e = edges();
    let named = |owner: &str, target: &str| {
        e.get(owner)
            .and_then(|r| r.get(target))
            .copied()
            .unwrap_or(0)
    };
    let mut broken = Vec::new();

    // `draft` is the vocabulary both halves speak, so it sits under them. It needs names
    // (`combinatorics`), the class tables a face lives in (`planes`) and the judge wrapper
    // (`tolerant`) -- and nothing above.
    for above in [
        "arrangement",
        "assembly",
        "boolean",
        "ops",
        "bands",
        "cyl_chart",
        "nesting",
    ] {
        let n = named("draft", above);
        if n > 0 {
            broken.push(format!(
                "draft -> {above} ({n}): the floor is reaching upward"
            ));
        }
    }
    // The assembly is the last stage. It reads what the engine handed it; it does not ask the
    // engine anything.
    let n = named("assembly", "arrangement");
    if n > 0 {
        broken.push(format!(
            "assembly -> arrangement ({n}): the last stage is calling the engine"
        ));
    }
    // Four modules live under the engine: the names and the class table it reads
    // (`combinatorics`, `planes`), and the two helpers it calls (`bands`, `nesting`), which
    // speak their own vocabulary back -- a `CylinderDef`, a `CylRow` -- rather than the
    // engine's. The engine reads all four; none of them knows it exists. A name pointing up
    // here means something below has started taking the caller's words again.
    //
    // `cyl_chart` is deliberately absent: it names `arrangement` 14 times, and a pair with no
    // answer yet cannot be asserted to zero.
    for below in ["bands", "combinatorics", "nesting", "planes"] {
        let n = named(below, "arrangement");
        if n > 0 {
            broken.push(format!(
                "{below} -> arrangement ({n}): something under the engine is naming the engine"
            ));
        }
    }
    // The class table is under the names that use it. `planes` holds per-face and
    // per-plane-class rows; `combinatorics` is what names points and edges from them. A table
    // calling the vocabulary built on it runs backwards -- which is what the plane setup did,
    // by holding each solid's edge incidence for a consumer it never read it for.
    let n = named("planes", "combinatorics");
    if n > 0 {
        broken.push(format!(
            "planes -> combinatorics ({n}): the class table is naming what is built on it"
        ));
    }
    // The front door is the first stage and lives above the engine, not inside it.
    let n = named("arrangement", "boolean");
    if n > 0 {
        broken.push(format!(
            "arrangement -> boolean ({n}): the engine is calling the front door"
        ));
    }

    assert!(
        broken.is_empty(),
        "the boolean pipeline runs boolean -> arrangement -> assembly, all of them reading \
         draft. These edges run the other way:\n  {}",
        broken.join("\n  ")
    );
}

/// **Everything else, frozen.** The seven surviving cycles are not claimed to be right -- they are
/// claimed to be *known*, and this is the list.
///
/// ★ When it fires there are two answers and no third. Either the edge goes away, or it joins the
/// table **and the commit body says why that direction is correct**. A gate that people only ever
/// add to is a list, not a gate.
#[test]
fn no_module_edge_appears_that_is_not_recorded() {
    // owner -> the modules it may name. Measured on the commit that split the front door from the
    // assembly; counts deliberately left out, because a count that must be updated to add a call
    // is a gate that gets edited until it is quiet.
    const RECORDED: &[(&str, &[&str])] = &[
        (
            "arrangement",
            &[
                "assembly",
                "bands",
                "combinatorics",
                "cyl_chart",
                "draft",
                "nesting",
                "par",
                "phase",
                "planes",
                "reuse",
                "tolerant",
                "transform",
            ],
        ),
        (
            "assembly",
            &[
                "combinatorics",
                "draft",
                "nesting",
                "planes",
                "realize",
                "tolerant",
                "transform",
            ],
        ),
        ("bands", &["combinatorics", "draft", "planes", "tolerant"]),
        (
            "boolean",
            &["arrangement", "assembly", "draft", "reject_census"],
        ),
        ("combinatorics", &["planes", "tolerant"]),
        (
            "cyl_chart",
            &[
                "arrangement",
                "assembly",
                "bands",
                "combinatorics",
                "draft",
                "planes",
                "tolerant",
            ],
        ),
        ("draft", &["combinatorics", "planes", "tolerant"]),
        ("error", &["reject_census"]),
        ("exact", &["ops", "rotated_vertex"]),
        ("nesting", &["combinatorics", "planes", "tolerant"]),
        (
            "ops",
            &[
                "boolean",
                "exact",
                "planes",
                "realize",
                "rotated_vertex",
                "transform",
            ],
        ),
        ("planes", &["par", "phase", "rotated_vertex"]),
        ("realize", &["planes", "rotated_vertex"]),
        (
            "reuse",
            &["combinatorics", "draft", "planes", "rotated_vertex"],
        ),
        ("tolerant", &["planes"]),
        ("transform", &["realize"]),
    ];
    let recorded: BTreeMap<&str, BTreeSet<&str>> = RECORDED
        .iter()
        .map(|(o, ts)| (*o, ts.iter().copied().collect()))
        .collect();
    let mut fresh = Vec::new();
    for (owner, row) in edges() {
        let known = recorded.get(owner.as_str());
        for (target, n) in row {
            if !known.is_some_and(|k| k.contains(target.as_str())) {
                fresh.push(format!("{owner} -> {target} ({n})"));
            }
        }
    }
    assert!(
        fresh.is_empty(),
        "a module edge that the table does not record:\n  {}\n\nRemove it, or add it to \
         RECORDED and say in the commit body why that direction is right.",
        fresh.join("\n  ")
    );
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
