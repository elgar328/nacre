//! **Every exemption from `clippy.toml`'s disallowed calls is listed here, with its reason.**
//!
//! `clippy.toml` names calls the workspace makes only in the places that own them — the f64 cache
//! write doors (product code reaches them through `nacre-ops`' realization funnel alone), rayon's
//! parallel entry points (`par.rs` alone), and the two ring predicates a shared rule sits above
//! (`point_in_ring` under `cell_inside`'s retry, `side_of` under `ring_against_plane`'s walk) — and
//! the pre-commit hook's `clippy -D warnings` refuses one anywhere else. An exemption is one attribute
//! (`allow`, `expect`, a `cfg_attr` around either), so a caller who meets the lint could add one
//! and be through. This makes that a second, visible step: every exemption in the workspace's
//! sources is counted per file against [`ALLOWED`], and a new one fails here until it is written
//! into the table with why it is not a way around the rule. A lint *group* that would silence the
//! lint without naming it (`clippy::all`, `clippy::style`, `warnings`) is held to zero everywhere.
//! [`the_config_still_names_every_call`] holds the other end: a path dropped from `clippy.toml`
//! is as quiet a way through as an exemption. And a level set outside the sources — a `[lints]`
//! table in a `Cargo.toml`, `rustflags` in a cargo config — silences the lint with no attribute
//! at all, measured: one `disallowed_methods = "allow"` line in `nacre-ops`'s manifest let a
//! planted door call through clippy and this table both; so none is set
//! ([`no_manifest_or_cargo_config_sets_the_lints_level`]).
//!
//! A second `clippy.toml` in a crate replaces the root's for it, so there is only the root's
//! ([`the_root_clippy_toml_is_the_only_one`]).
//!
//! The lock binds the commits the pre-commit hook runs on. `--no-verify`, an edited hook (it is
//! not versioned), a user's own `~/.cargo/config.toml` or `CLIPPY_CONF_DIR` in the environment are
//! past what a test in the tree can see.
//!
//! Textual, like `comparators.rs`: it reads the attribute's spelling, which is fixed, and
//! [`the_scan_reads_what_it_should`] pins that reading on literal lines.

use std::collections::{BTreeMap, BTreeSet};

/// Every path `clippy.toml` disallows.
const DISALLOWED: [&str; 19] = [
    "nacre_topo::Model::push_vertex",
    "nacre_topo::Model::refine_vertex_cache",
    "nacre_topo::Model::push_plane",
    "nacre_topo::Model::push_plane_through",
    "nacre_topo::Model::push_cylinder",
    "nacre_topo::Model::refine_surface_cache",
    "nacre_topo::Model::push_edge",
    "nacre_topo::Model::rebuild_edge_cache",
    "nacre_topo::Model::push_plane_unregistered",
    "nacre_topo::Model::push_cylinder_unregistered",
    "nacre_ops::realize::push_vertex_prerealized",
    "rayon::iter::IntoParallelIterator::into_par_iter",
    "rayon::iter::IntoParallelRefIterator::par_iter",
    "rayon::slice::ParallelSlice::par_chunks",
    "rayon::iter::ParallelBridge::par_bridge",
    "rayon::join",
    "rayon::scope",
    "nacre_ops::combinatorics::ring_walk::point_in_ring",
    "nacre_ops::combinatorics::ring_walk::side_of",
];

/// Every exemption, per file (relative to `crates/`): how many, and why it is not a way around.
const ALLOWED: &[(&str, usize, &str)] = &[
    (
        "nacre-ops/src/realize.rs",
        1,
        "the cache write doors: this module is the realization funnel (a module `expect`)",
    ),
    (
        "nacre-topo/src/model/mod.rs",
        1,
        "the cache write doors: the empty model's three world planes, whose names derive their \
         cache (a statement `expect`, and a `debug_assert` that the cache is `Realized`)",
    ),
    (
        "nacre-ops/src/assembly/reconstruct.rs",
        1,
        "`push_vertex_prerealized`: the seam table's cache, which `realize::point_cache` gave for \
         the very definition pushed",
    ),
    (
        "nacre-ops/src/par.rs",
        1,
        "rayon: this module is where parallel iteration lives (a module `expect`, parallel \
         build only)",
    ),
    (
        "nacre-ops/src/nesting.rs",
        1,
        "`point_in_ring`: `ask`, one witness's answer; the retry over witnesses is \
         `cell_inside`'s loop",
    ),
    (
        "nacre-ops/src/combinatorics/ring_walk.rs",
        1,
        "`side_of`: `ring_against_plane`, the ring walk itself",
    ),
    (
        "nacre-ops/src/combinatorics/rays.rs",
        1,
        "`side_of`: a point question in `point_on_ring` (is the literal vertex on the edge's line)",
    ),
    (
        "nacre-ops/src/combinatorics/in_faces.rs",
        1,
        "`side_of`: a point question in `point_in_component` (does the query lie on the face's \
         plane)",
    ),
    (
        "nacre-ops/src/arrangement/aliases.rs",
        1,
        "`side_of`: a point question in `seed_from_operands` (does a pierce corner lie on a class)",
    ),
    (
        "nacre-ops/src/assembly/naming.rs",
        1,
        "`side_of`: a point question in `name_result_vertices` (is a pierce point on the other \
         class of a shared ruling)",
    ),
    (
        "nacre-ops/src/assembly/self_touch.rs",
        1,
        "`side_of`: a point question in `self_touch_reject` (does a result node lie on a face's \
         plane)",
    ),
    (
        "nacre-ops/src/lib.rs",
        1,
        "test builds only (`cfg_attr(test, …)`): fixtures build models through the doors; the \
         library is linted without it",
    ),
    (
        "nacre-topo/src/lib.rs",
        1,
        "test builds only, as in nacre-ops",
    ),
    (
        "nacre-validate/src/lib.rs",
        1,
        "test builds only, as in nacre-ops",
    ),
    (
        "nacre-ops/tests/instruments.rs",
        1,
        "an integration-test binary: fixtures build models through the doors",
    ),
    (
        "nacre-ops/tests/invariants.rs",
        1,
        "an integration-test binary, as above",
    ),
    (
        "nacre-topo/tests/cylinder_truth.rs",
        1,
        "an integration-test binary, as above",
    ),
];

/// `t` with its string literals emptied and any `//` comment cut, so a lint named in a reason,
/// a message or a comment is not an exemption. A char literal (`'"'`, `'\"'`) is kept whole, so
/// its quote does not open a string that would swallow the rest of the line.
fn code_of(t: &str) -> String {
    let chars: Vec<char> = t.chars().collect();
    let mut out = String::new();
    let (mut in_str, mut escaped) = (false, false);
    let mut k = 0;
    while k < chars.len() {
        let c = chars[k];
        k += 1;
        if in_str {
            match (escaped, c) {
                (true, _) => escaped = false,
                (false, '\\') => escaped = true,
                (false, '"') => {
                    in_str = false;
                    out.push('"');
                }
                _ => {}
            }
            continue;
        }
        match c {
            '\'' => {
                // `'x'` or `'\x'`: copy it through; anything else (a lifetime) is one char.
                let len = match chars.get(k..) {
                    Some([x, '\'', ..]) if *x != '\\' => 2,
                    Some(['\\', _, '\'', ..]) => 3,
                    _ => 0,
                };
                out.push(c);
                out.extend(&chars[k..k + len]);
                k += len;
            }
            '"' => {
                in_str = true;
                out.push('"');
            }
            '/' if chars.get(k) == Some(&'/') => break,
            _ => out.push(c),
        }
    }
    out
}

/// Whether `name` occurs in `code` as a whole path (not the tail or head of a longer one).
fn has_path(code: &str, name: &str) -> bool {
    let ident = |c: char| c.is_ascii_alphanumeric() || c == '_' || c == ':';
    code.match_indices(name).any(|(at, _)| {
        !code[..at].chars().next_back().is_some_and(ident)
            && !code[at + name.len()..].chars().next().is_some_and(ident)
    })
}

/// How many exemptions line `t` makes: the lint named (`allow`, `expect`, a `cfg_attr` around
/// either — the name is the tell, whatever wraps it). The group spellings are counted apart.
fn exemptions_in(t: &str) -> usize {
    usize::from(has_path(&code_of(t), "clippy::disallowed_methods"))
}

/// Whether line `t` silences a whole lint group that holds `disallowed_methods`.
fn silences_a_group(t: &str) -> bool {
    let code = code_of(t);
    let in_attr = code.contains("allow(") || code.contains("expect(");
    in_attr
        && (has_path(&code, "clippy::all")
            || has_path(&code, "clippy::style")
            || code.contains("allow(warnings")
            || code.contains("expect(warnings"))
}

/// Every `.rs` under `crates/*/src` and `crates/*/tests` (relative to `crates/`): tests are where
/// the doors are legitimately called, so their exemptions are counted too.
fn workspace_sources() -> Vec<(String, std::path::PathBuf)> {
    let crates = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf();
    let mut dirs: Vec<std::path::PathBuf> = Vec::new();
    for entry in std::fs::read_dir(&crates).expect("crates/") {
        let krate = entry.expect("entry").path();
        for sub in ["src", "tests"] {
            if krate.join(sub).is_dir() {
                dirs.push(krate.join(sub));
            }
        }
    }
    let mut out = Vec::new();
    while let Some(dir) = dirs.pop() {
        for entry in std::fs::read_dir(&dir).expect("dir") {
            let path = entry.expect("entry").path();
            if path.is_dir() {
                dirs.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") {
                let rel = path
                    .strip_prefix(&crates)
                    .expect("under crates/")
                    .to_string_lossy()
                    .replace('\\', "/");
                out.push((rel, path));
            }
        }
    }
    out.sort();
    out
}

#[test]
fn every_exemption_from_a_disallowed_call_is_listed() {
    let mut found: BTreeMap<String, usize> = BTreeMap::new();
    let mut groups = Vec::new();
    for (rel, path) in workspace_sources() {
        let text = std::fs::read_to_string(&path).expect("read");
        let mut n = 0;
        for (k, line) in text.lines().enumerate() {
            n += exemptions_in(line);
            if silences_a_group(line) {
                groups.push(format!("{rel}:{}: {}", k + 1, line.trim()));
            }
        }
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
        "an exemption from clippy.toml's disallowed calls outside the list -- don't add one: \
         take the road the lint's reason names (for a cache write door, nacre-ops' realization \
         funnel in realize.rs). If a new place truly owns the call, add it to ALLOWED with why \
         (and drop an entry whose exemption is gone)"
    );
    assert!(
        groups.is_empty(),
        "a lint group silenced, which turns the disallowed calls off without naming them: \
         {groups:?}"
    );
}

/// `clippy.toml` still disallows every path, and nothing else: dropping one is as quiet a way
/// through as an exemption, and one added there belongs in [`DISALLOWED`] with its exemptions here.
#[test]
fn the_config_still_names_every_call() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|c| c.parent())
        .expect("the workspace root");
    let config = std::fs::read_to_string(root.join("clippy.toml")).expect("clippy.toml");
    let named: BTreeSet<&str> = config
        .lines()
        .filter(|l| !l.trim_start().starts_with('#'))
        .filter_map(|l| l.split("path = \"").nth(1))
        .filter_map(|rest| rest.split('"').next())
        .collect();
    let want: BTreeSet<&str> = DISALLOWED.into_iter().collect();
    assert_eq!(named, want, "clippy.toml's disallowed paths moved");
}

/// Whether line `t` of a manifest or cargo config sets the lint's level (or a group's) — a
/// `[lints]` table, a `[workspace.lints]` one, or `rustflags` — which silences the disallowed
/// calls with no attribute in any source. `table` is the header the line sits under.
fn manifest_silences(table: &str, t: &str) -> bool {
    let code = t.split('#').next().unwrap_or("");
    let word = |w: &str| {
        code.match_indices(w).any(|(at, _)| {
            let edge = |c: Option<char>| {
                !c.is_some_and(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
            };
            edge(code[..at].chars().next_back()) && edge(code[at + w.len()..].chars().next())
        })
    };
    let names_the_lint = word("disallowed_methods") || word("disallowed-methods");
    let in_lints = table.starts_with("[lints") || table.starts_with("[workspace.lints");
    let groups = word("all") || word("style") || word("warnings");
    // Any `rustflags` at all: its array may span lines, and the flags that matter (`-A …`) would
    // then sit on lines this reads without the key.
    let flags = word("rustflags") || code.contains("cap-lints");
    names_the_lint || (in_lints && groups) || flags
}

/// The manifests and cargo configs a lint level can be set in: the root's and every crate's
/// `Cargo.toml`, and `.cargo/config(.toml)` if one is ever added.
fn lint_level_files() -> Vec<std::path::PathBuf> {
    let crates = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf();
    let root = crates.parent().expect("the workspace root").to_path_buf();
    let mut out = vec![
        root.join("Cargo.toml"),
        root.join(".cargo/config.toml"),
        root.join(".cargo/config"),
    ];
    for entry in std::fs::read_dir(&crates).expect("crates/") {
        out.push(entry.expect("entry").path().join("Cargo.toml"));
    }
    out.retain(|p| p.is_file());
    out.sort();
    out
}

/// A level set in a manifest or a cargo config is an exemption no source shows, so there is none.
#[test]
fn no_manifest_or_cargo_config_sets_the_lints_level() {
    let mut offenders = Vec::new();
    for path in lint_level_files() {
        let text = std::fs::read_to_string(&path).expect("read");
        let mut table = String::new();
        for (k, line) in text.lines().enumerate() {
            let t = line.trim();
            if t.starts_with('[') {
                table = t.to_string();
            } else if manifest_silences(&table, t) {
                offenders.push(format!("{}:{}: {t}", path.display(), k + 1));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "a manifest or cargo config sets the level of clippy.toml's disallowed calls (or of a \
         group holding them), which no exemption table sees -- set it nowhere: {offenders:?}"
    );
}

/// The root's `clippy.toml` is the only one. clippy reads the nearest config from a crate's own
/// folder upward, so a `clippy.toml` in a crate (or in `crates/`) — even an empty one — replaces
/// the root's for that crate and drops every disallowed path with it, measured: an empty
/// `crates/nacre-ops/clippy.toml` let a planted door call through clippy and every other test here.
#[test]
fn the_root_clippy_toml_is_the_only_one() {
    let crates = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf();
    let root = crates.parent().expect("the workspace root").to_path_buf();
    let mut dirs = vec![crates.clone()];
    for entry in std::fs::read_dir(&crates).expect("crates/") {
        let path = entry.expect("entry").path();
        if path.is_dir() {
            dirs.push(path);
        }
    }
    let mut shadows: Vec<_> = dirs
        .iter()
        .flat_map(|d| [d.join("clippy.toml"), d.join(".clippy.toml")])
        .collect();
    shadows.push(root.join(".clippy.toml"));
    shadows.retain(|p| p.exists());
    assert!(
        shadows.is_empty(),
        "a clippy config besides the root's, which replaces it (and its disallowed calls) for \
         the crates under it: {shadows:?}"
    );
}

/// ★ **The rule is exercised on strings, not on the tree** — a scan that only ever sees a tree it
/// agrees with has not been shown to say anything.
#[test]
fn the_scan_reads_what_it_should() {
    assert_eq!(exemptions_in("#[allow(clippy::disallowed_methods)]"), 1);
    assert_eq!(exemptions_in("#![allow(clippy::disallowed_methods)]"), 1);
    assert_eq!(
        exemptions_in("#![cfg_attr(test, allow(clippy::disallowed_methods))]"),
        1,
        "inside a cfg_attr"
    );
    assert_eq!(
        exemptions_in("    clippy::disallowed_methods,"),
        1,
        "a multi-line `expect`'s lint line"
    );
    assert_eq!(
        exemptions_in("#[allow(clippy::too_many_lines, clippy::disallowed_methods)]"),
        1,
        "beside another lint"
    );
    assert_eq!(
        exemptions_in("// clippy::disallowed_methods is the lint"),
        0
    );
    assert_eq!(
        exemptions_in("/// `clippy::disallowed_methods` in a doc"),
        0
    );
    assert_eq!(
        exemptions_in("assert!(x, \"no clippy::disallowed_methods here\");"),
        0,
        "inside a string"
    );
    assert_eq!(
        exemptions_in("#[allow(clippy::disallowed_methods_extra)]"),
        0,
        "a longer name"
    );
    assert!(silences_a_group("#[allow(clippy::all)]"));
    assert!(silences_a_group("#![allow(clippy::style)]"));
    assert!(silences_a_group("#![cfg_attr(test, allow(warnings))]"));
    assert!(silences_a_group("#[expect(clippy::style, reason = \"x\")]"));
    assert!(!silences_a_group("#[allow(clippy::all_the_things)]"));
    assert!(!silences_a_group("// allow(clippy::all) would hide it"));
    assert!(
        !silences_a_group("let style = clippy::style;"),
        "not in an attribute"
    );
    assert_eq!(
        exemptions_in("let q = '\"'; #[allow(clippy::disallowed_methods)]"),
        1,
        "a char literal's quote does not open a string"
    );
    assert_eq!(
        exemptions_in("let q = '\\\"'; #[allow(clippy::disallowed_methods)]"),
        1,
        "nor an escaped one"
    );
    assert_eq!(
        exemptions_in("fn f<'a>(x: &'a str) {} #[allow(clippy::disallowed_methods)]"),
        1,
        "a lifetime is not a char literal"
    );

    assert!(manifest_silences(
        "[lints.clippy]",
        "disallowed_methods = \"allow\""
    ));
    assert!(manifest_silences(
        "[workspace.lints.clippy]",
        "all = { level = \"allow\", priority = -1 }"
    ));
    assert!(manifest_silences(
        "[lints]",
        "clippy = { style = \"allow\" }"
    ));
    assert!(manifest_silences(
        "[build]",
        "rustflags = [\"-A\", \"clippy::disallowed_methods\"]"
    ));
    assert!(manifest_silences(
        "[build]",
        "rustflags = [\"--cap-lints\", \"warn\", \"-A\", \"warnings\"]"
    ));
    assert!(
        manifest_silences("[build]", "rustflags = ["),
        "an array spanning lines is caught at its key"
    );
    assert!(
        !manifest_silences("[package]", "name = \"rustflags-free\""),
        "not the key"
    );
    assert!(!manifest_silences(
        "[workspace.lints.rust]",
        "unsafe_code = \"forbid\""
    ));
    assert!(
        !manifest_silences("[features]", "all = [\"a\", \"b\"]"),
        "a feature named all"
    );
    assert!(
        !manifest_silences("[lints.clippy]", "# all = \"allow\" once"),
        "a comment"
    );
}
