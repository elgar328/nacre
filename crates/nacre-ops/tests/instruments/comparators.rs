//! **A sort handed a comparison is handed one that cannot fail — or says why it is safe.**
//!
//! A comparison that reads «could not tell» as `Equal` is not a total order, and since Rust 1.81
//! the standard sorts may **panic** on one past twenty elements (measured: 1,896 of 2,000 sorts of
//! forty, one element in ten failing every comparison it took part in; and on the arrangement's own
//! split sort, a single failing point). That idiom had been copied to seven places in `nacre-ops`,
//! one of them with a comment calling the result «any permutation». Where an exact comparison can
//! fail, `combinatorics::try_sort_by` sorts and stops at the failure.
//!
//! So every product use of the entries a comparison rides into — `sort_by`, `sort_unstable_by`,
//! `select_nth_unstable_by` and their rayon twins, and a hand-written `Ord`/`PartialOrd` (which
//! `sort`, `sort_by_key` and a `BTreeMap` all read) — is on [`ALLOWED`] with the reason it is total,
//! counted per file so a second one in an allowed file is seen. The rest of the workspace is held
//! too: the panic does not care which crate a sort is in, and every product crate ships to wasm.
//!
//! Being textual it can be talked around (a comparator passed through a generic helper is not a
//! call this reads). It is a tripwire on the ordinary way to write the defect, not a proof — so
//! [`the_scan_reads_what_it_should`] pins the rule against literal lines rather than against
//! files, which would go red whenever a file moves and mean nothing when it did.

use std::collections::BTreeMap;

/// The methods a comparison closure rides into, matched as whole identifiers before `(`.
const ENTRIES: [&str; 5] = [
    "sort_by",
    "sort_unstable_by",
    "select_nth_unstable_by",
    "par_sort_by",
    "par_sort_unstable_by",
];

/// Every product use, per file (relative to `crates/`), with why its comparison is a total order.
const ALLOWED: &[(&str, usize, &str)] = &[
    (
        "nacre-ops/src/arrangement/rulings.rs",
        1,
        "an exact `Rat` axis parameter, then the derived order of the event",
    ),
    (
        "nacre-ops/src/assembly/grouping.rs",
        1,
        "f64 coordinates under the finite-vertex invariant, which the `expect` states",
    ),
    (
        "nacre-ops/src/assembly/reconstruct.rs",
        1,
        "the same coordinate keys as `grouping`'s, under the same `expect`",
    ),
    (
        "nacre-ops/src/bands.rs",
        1,
        "a class index, then an exact `Rat` span end",
    ),
    (
        "nacre-store/src/lib.rs",
        2,
        "`Handle`'s `Ord`/`PartialOrd`: its integer index",
    ),
    ("nacre-tess/src/lib.rs", 3, "`f64::total_cmp`"),
    (
        "nacre-tess/src/polygon/monotone.rs",
        2,
        "the sweep's order: a tie it cannot decide is `Equal`, and `self_touch` runs first and \
         refuses every coincidence but one sanctioned pair of twins, so a tie is between two \
         elements and consistent; the outgoing edges: half-plane, then exact `orient2d`",
    ),
];

/// How many entries line `t` (trimmed) makes: calls of an [`ENTRIES`] method as a whole
/// identifier, plus one for a hand-written `Ord`/`PartialOrd`.
fn entries_in(t: &str) -> usize {
    if t.starts_with("//") {
        return 0;
    }
    let ident = |c: char| c.is_ascii_alphanumeric() || c == '_';
    let mut n = 0;
    for name in ENTRIES {
        let call = format!("{name}(");
        let mut from = 0;
        while let Some(at) = t[from..].find(&call) {
            let start = from + at;
            if !t[..start].chars().next_back().is_some_and(ident) {
                n += 1;
            }
            from = start + call.len();
        }
    }
    if t.starts_with("impl") && (t.contains(" Ord for ") || t.contains(" PartialOrd for ")) {
        n += 1;
    }
    n
}

/// Whether `lines[k]` opens a test module written inline — the scan stops there, as a product
/// file's inline tests are not product code. A `#[cfg(test)]` over a module declared elsewhere
/// (`mod tests;`, `#[path]`) is not one.
fn inline_tests_at(lines: &[&str], k: usize) -> bool {
    lines[k].trim() == "#[cfg(test)]"
        && lines.get(k + 1).is_some_and(|n| {
            let n = n.trim();
            n.starts_with("mod ") && n.ends_with('{')
        })
}

/// Entries in one file's text.
fn entries_in_file(text: &str) -> usize {
    let lines: Vec<&str> = text.lines().collect();
    let mut n = 0;
    for k in 0..lines.len() {
        if inline_tests_at(&lines, k) {
            break;
        }
        n += entries_in(lines[k].trim());
    }
    n
}

/// Every product source of the workspace: `crates/*/src/**/*.rs`, minus any `tests` folder.
fn product_sources() -> Vec<(String, std::path::PathBuf)> {
    let crates = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .to_path_buf();
    let mut out = Vec::new();
    for entry in std::fs::read_dir(&crates).expect("crates/") {
        let src = entry.expect("entry").path().join("src");
        if !src.is_dir() {
            continue;
        }
        for path in crate::walk::files_under(&src, |name| name != "tests") {
            if path.extension().is_some_and(|e| e == "rs") {
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
fn every_comparison_a_sort_is_handed_is_total_or_listed() {
    let mut found: BTreeMap<String, usize> = BTreeMap::new();
    for (rel, path) in product_sources() {
        let n = entries_in_file(&std::fs::read_to_string(&path).expect("read"));
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
        "a sort is handed a comparison, or an `Ord` is written by hand, outside the list -- if \
         the comparison can fail, sort with `combinatorics::try_sort_by`, which stops at the \
         failure; if it cannot, add the file and count to ALLOWED with the reason it is total \
         (and drop an entry whose sort is gone)"
    );
}

/// **The rule is exercised on strings, not on the tree** — a scan that only ever sees a tree it
/// agrees with has not been shown to say anything.
#[test]
fn the_scan_reads_what_it_should() {
    assert_eq!(entries_in("v.sort_by(|a, b| a.cmp(b));"), 1);
    assert_eq!(entries_in("v.sort_unstable_by(f);"), 1);
    assert_eq!(entries_in("v.select_nth_unstable_by(3, f);"), 1);
    assert_eq!(entries_in("v.par_sort_by(f);"), 1, "rayon's twin counts");
    assert_eq!(
        entries_in("a.sort_by(f); b.sort_by(g);"),
        2,
        "two on a line"
    );
    assert_eq!(
        entries_in("combinatorics::try_sort_by(&mut v, f)?;"),
        0,
        "the sort that stops at a failure is not one"
    );
    assert_eq!(
        entries_in("v.sort_by_key(|a| a.0);"),
        0,
        "a key sort reads `Ord`"
    );
    assert_eq!(entries_in("// v.sort_by(f) would panic"), 0);
    assert_eq!(entries_in("/// `sort_by` with a flag."), 0);
    assert_eq!(entries_in("impl<T> Ord for Handle<T> {"), 1);
    assert_eq!(entries_in("impl PartialOrd for Key {"), 1);
    assert_eq!(
        entries_in("impl Ordinal for Key {"),
        0,
        "a trait whose name starts with Ord"
    );

    let inline = [
        "fn f() {}",
        "#[cfg(test)]",
        "mod tests {",
        "    v.sort_by(f);",
        "}",
    ];
    assert_eq!(
        entries_in_file(&inline.join("\n")),
        0,
        "an inline test module is not product code"
    );
    let mounted = [
        "#[cfg(test)]",
        "#[path = \"../tests/x.rs\"]",
        "mod tests;",
        "v.sort_by(f);",
    ];
    assert_eq!(
        entries_in_file(&mounted.join("\n")),
        1,
        "a test module mounted from elsewhere does not end the product code"
    );
}
