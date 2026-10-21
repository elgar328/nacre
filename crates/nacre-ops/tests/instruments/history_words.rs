//! **No comment or document tells how the code got here — the sentence says what holds today.**
//!
//! overview's document rules 1, 6 and 10: a comment states the present, and the story of how it
//! came to be lives in the commit message. A story has a few fixed spellings, and this file holds
//! them ([`PHRASES`], [`PHRASES_KO`]) together with three shapes: a work-unit cell number, a family
//! code, and a commit hash. The examples live only in string literals, so this file does not trip
//! itself.
//!
//! What it reads: every `.rs` and `.md` file [`crate::markers::prose_files`] reads, except
//! `CHANGELOG.md`, whose content is history. In a `.rs` file only the comment part of a line counts
//! ([`crate::lex::split_line`]), so a message or a string is not prose; a `.md` file counts whole.
//! `#` comments (shell, Python, toml, the hook) are not read. In `docs/design.md` the section of
//! rejected approaches is skipped, because rule 2 puts the history and the commits of a rejected
//! approach exactly there; the section is found by its heading ([`STORY_SECTION`]), so renaming the
//! heading turns this red rather than quietly reading less.
//!
//! Being textual it is a tripwire, not a proof: a story told without these spellings is judged by
//! rules 1 and 6. [`the_scan_reads_what_it_should`] pins the rule against literal lines.

use std::collections::BTreeMap;

/// Spellings that only a story uses, matched in lower case at word boundaries.
///
/// Spellings with a present-tense use are left out on purpose: the plain verb of purpose ("a box
/// used to reject"), "no longer" ("lives no longer than one operation"), "the fix" ("the fix for an
/// exceeded bound is …"), "previously", "originally".
const PHRASES: &[&str] = &[
    "it used to",
    "this used to",
    "that used to",
    "which used to",
    "they used to",
    "first draft",
    "earlier draft",
    "earlier version",
    "previous version",
    "old version",
    "earlier spelling",
    "before the fix",
    "before this fix",
    "learned to",
    "once broke",
    "once failed",
    "formerly",
    "has since",
    "have since",
];

/// The Korean spellings, matched at a left boundary only: a particle follows the word directly and
/// Hangul is alphanumeric, so a right boundary would miss the common spelling.
const PHRASES_KO: &[&str] = &["예전에", "과거에는"];

/// The heading of `docs/design.md`'s section of rejected approaches.
const STORY_SECTION: &str = "## 가지 말 것";

/// Lines that may tell a story, per file (relative to the repository root), with why.
const ALLOWED: &[(&str, usize, &str)] = &[];

fn word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// Whether `p` occurs in `lower` with no word character before it (and, when `right`, after it).
fn has_word(lower: &str, p: &str, right: bool) -> bool {
    lower.match_indices(p).any(|(at, _)| {
        !lower[..at].chars().next_back().is_some_and(word)
            && !(right && lower[at + p.len()..].chars().next().is_some_and(word))
    })
}

/// A work unit's number: "cell", then an optional dash or space, then a digit or a circled number.
fn has_cell_number(lower: &str) -> bool {
    lower.match_indices("cell").any(|(at, _)| {
        if lower[..at].chars().next_back().is_some_and(word) {
            return false;
        }
        let mut rest = lower[at + 4..].chars();
        let next = match rest.next() {
            Some('-' | ' ') => rest.next(),
            other => other,
        };
        next.is_some_and(|c| c.is_ascii_digit() || ('\u{2460}'..='\u{2473}').contains(&c))
    })
}

/// A family code: "family #" and a digit.
fn has_family_code(lower: &str) -> bool {
    lower.match_indices("family #").any(|(at, p)| {
        lower[at + p.len()..]
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_digit())
    })
}

/// A commit hash: a token of 7–40 lower-case hex digits holding both a letter and a digit. Tokens
/// run over `[A-Za-z0-9_.-]` with a `.` or `-` at either end trimmed (a hash ends sentences), so a
/// number in exponent notation stays one token with a `.` or `-` inside it; a bare `1e10000` is a
/// number too.
fn has_hash(text: &str) -> bool {
    text.split(|c: char| !(c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-')))
        .map(|t| t.trim_matches(|c| matches!(c, '.' | '-')))
        .any(|t| {
            let hex = t.chars().all(|c| matches!(c, '0'..='9' | 'a'..='f'));
            let exponent = t.split_once('e').is_some_and(|(a, b)| {
                !a.is_empty()
                    && !b.is_empty()
                    && (a.to_owned() + b).chars().all(|c| c.is_ascii_digit())
            });
            (7..=40).contains(&t.len())
                && hex
                && t.chars().any(|c| c.is_ascii_digit())
                && t.chars().any(|c| c.is_ascii_lowercase())
                && !exponent
        })
}

fn tells_a_story(text: &str) -> bool {
    let lower = text.to_lowercase();
    PHRASES.iter().any(|p| has_word(&lower, p, true))
        || PHRASES_KO.iter().any(|p| has_word(&lower, p, false))
        || has_cell_number(&lower)
        || has_family_code(&lower)
        || has_hash(text)
}

/// How many lines of `text` tell a story, read as the file `rel` (relative to the repository root).
fn findings(rel: &str, text: &str) -> usize {
    let rs = rel.ends_with(".rs");
    let design = rel == "docs/design.md";
    let mut skipping = false;
    let mut n = 0;
    for line in text.lines() {
        if design && line.starts_with("## ") {
            skipping = line.starts_with(STORY_SECTION);
        }
        if skipping {
            continue;
        }
        let prose = if rs {
            match crate::lex::split_line(line).1 {
                Some(comment) => comment,
                None => continue,
            }
        } else {
            line.to_owned()
        };
        if tells_a_story(&prose) {
            n += 1;
        }
    }
    n
}

#[test]
fn no_comment_or_document_tells_how_the_code_got_here() {
    let mut found: BTreeMap<String, usize> = BTreeMap::new();
    for (rel, path) in crate::markers::prose_files() {
        if !(rel.ends_with(".rs") || rel.ends_with(".md")) || rel == "CHANGELOG.md" {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("{rel} does not read as UTF-8 text: {e}"));
        let n = findings(&rel, &text);
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
        "a comment or document tells how the code got here -- say what holds today (overview's \
         document rules 1, 6 and 10): the story goes in the commit message, a rejected approach in \
         design.md's table; the list holds only lines that must keep one (and drop an entry that is \
         gone)"
    );
}

/// **The rule is exercised on strings, not on the tree** — a scan that only ever sees a tree it
/// agrees with has not been shown to say anything.
#[test]
fn the_scan_reads_what_it_should() {
    let rs = |line: &str| findings("crates/x/src/lib.rs", line);
    for p in PHRASES {
        assert_eq!(
            rs(&format!("/// {p} the other way.")),
            1,
            "{p:?} at a comment's head"
        );
        assert_eq!(
            rs(&format!("let a = 1; // and {p} there")),
            1,
            "{p:?} mid-comment"
        );
    }
    assert_eq!(
        rs("//! This used to read otherwise."),
        1,
        "case does not hide a spelling"
    );
    for p in PHRASES_KO {
        assert_eq!(rs(&format!("/// {p} 그랬다.")), 1, "{p:?}");
    }
    assert_eq!(
        rs("/// 예전에는 그랬다."),
        1,
        "a particle right after the word"
    );
    for shape in [
        "/// cell 3 opens it.",
        "/// cell-13 left it.",
        "/// cell \u{2469} did.",
    ] {
        assert_eq!(rs(shape), 1, "a cell number: {shape}");
    }
    assert_eq!(rs("/// family #3 did."), 1, "a family code");
    assert_eq!(rs("/// See 0dbab39."), 1, "a commit hash");
    for fine in [
        "/// a box used to reject a candidate",
        "/// it lives no longer than one operation",
        "/// the fix for an exceeded bound is to find it",
        "/// the cell's chamber",
        "/// the cell reader",
        "/// the *wrong* cell.",
        "/// whatever it learned at 128 bits",
        "/// the previous cache moved",
        "/// an earlier stage's record",
        "/// before the change of frame",
        "/// `r = 5.000000000000001e-8`",
        "/// `(-2.0, 2.4492935982947064e-16)`",
        "/// 1e10000",
    ] {
        assert_eq!(rs(fine), 0, "present tense: {fine}");
    }
    assert_eq!(
        rs("let s = \"// it used to\";"),
        0,
        "a string is not a comment"
    );
    assert_eq!(
        rs("let s = \"x\"; // it used to"),
        1,
        "the comment after a string is"
    );
    assert_eq!(
        findings("docs/todo.md", "it used to be"),
        1,
        "a document counts whole"
    );
    let design = "## Shape\nsee 0dbab39\n## 가지 말 것 — 재 보고 버린 길\n| a | earlier version | `0dbab39` |\n";
    assert_eq!(
        findings("docs/design.md", design),
        1,
        "the rejected approaches are skipped"
    );
    assert_eq!(findings("docs/todo.md", design), 2, "only in design.md");
}
