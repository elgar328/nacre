//! **The docs a session reads stay inside a size budget** — overview's document rule 11.
//!
//! The parts are the units overview's reading rule reads, so what a session pays grows with their
//! size. Adding a sentence never turns anything red by itself, so without a number nothing says a
//! part grew. Each part here has a ceiling, counted in characters (a token follows characters more
//! closely than bytes: a Hangul syllable is three bytes).
//!
//! Three parts also have a floor at nine tenths of the ceiling, so room freed by trimming is given
//! back rather than left to refill unseen. Neither todo part has one, because todo is a queue whose
//! finished items are deleted (rule 3); nor does the section of rejected approaches, which grows a
//! row per rejected approach (rule 2).
//!
//! design is split at two headings, each of which must start exactly one line. The second is
//! [`crate::history_words::STORY_SECTION`], matched as a prefix there and here. That section runs
//! to the end of the file here, so a section added after it is counted with the rejected
//! approaches; `history_words` reads again from the next `## ` heading. todo is split at one
//! heading ([`PARKED`]), which must also start exactly one line; the parked items run from it to the
//! end of the file.

use std::path::Path;

/// One part: its name, its file, the heading it starts at and the one it ends before (`None` = the
/// file's start or end), its ceiling in characters, and whether it has a floor.
struct Part {
    name: &'static str,
    file: &'static str,
    from: Option<&'static str>,
    to: Option<&'static str>,
    ceiling: usize,
    floored: bool,
}

const LAYERS: &str = "# 층과 하위 체계";

/// todo's parked items: a session reads the headings and the items its work touches, not the whole.
const PARKED: &str = "## 보류 — 인구가 생기면";

const PARTS: &[Part] = &[
    Part {
        name: "overview",
        file: "docs/overview.md",
        from: None,
        to: None,
        ceiling: 20_100,
        floored: true,
    },
    Part {
        name: "design: the always-read rules part",
        file: "docs/design.md",
        from: None,
        to: Some(LAYERS),
        ceiling: 22_800,
        floored: true,
    },
    Part {
        name: "design: the body",
        file: "docs/design.md",
        from: Some(LAYERS),
        to: Some(crate::history_words::STORY_SECTION),
        ceiling: 171_400,
        floored: true,
    },
    Part {
        name: "design: rejected approaches",
        file: "docs/design.md",
        from: Some(crate::history_words::STORY_SECTION),
        to: None,
        ceiling: 30_600,
        floored: false,
    },
    Part {
        name: "todo: the always-read part",
        file: "docs/todo.md",
        from: None,
        to: Some(PARKED),
        ceiling: 21_200,
        floored: false,
    },
    Part {
        name: "todo: parked items",
        file: "docs/todo.md",
        from: Some(PARKED),
        to: None,
        ceiling: 13_500,
        floored: false,
    },
];

/// The byte offset of the one line of `text` that starts with `heading`, or how many lines do
/// when that is not exactly one.
fn heading_at(text: &str, heading: &str) -> Result<usize, usize> {
    let mut at = Vec::new();
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if line.starts_with(heading) {
            at.push(offset);
        }
        offset += line.len();
    }
    match at.as_slice() {
        [one] => Ok(*one),
        _ => Err(at.len()),
    }
}

/// The part of `text` from the line starting with `from` up to the line starting with `to`.
fn slice<'a>(text: &'a str, from: Option<&str>, to: Option<&str>) -> &'a str {
    let at = |h: &str| {
        heading_at(text, h)
            .unwrap_or_else(|n| panic!("{n} lines start with {h:?}; the split needs exactly one"))
    };
    let start = from.map_or(0, at);
    let end = to.map_or(text.len(), at);
    assert!(start <= end, "{from:?} comes after {to:?}");
    &text[start..end]
}

#[test]
fn each_doc_part_is_inside_its_budget() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("the repository root, two above this crate");
    let mut table = String::new();
    let mut off = false;
    for part in PARTS {
        let text = std::fs::read_to_string(root.join(part.file))
            .unwrap_or_else(|e| panic!("{} does not read: {e}", part.file));
        let size = slice(&text, part.from, part.to).chars().count();
        let floor = part.floored.then_some(part.ceiling / 10 * 9);
        let verdict = if size > part.ceiling {
            "over the ceiling"
        } else if floor.is_some_and(|f| size < f) {
            "under the floor: lower the ceiling"
        } else {
            "ok"
        };
        off |= verdict != "ok";
        let floor = floor.map_or("none".to_string(), |f| f.to_string());
        table.push_str(&format!(
            "  {:<36} {size:>7} chars, ceiling {:>7}, floor {floor:>7}  {verdict}\n",
            part.name, part.ceiling
        ));
    }
    assert!(
        !off,
        "a doc part is outside its budget -- see overview's document rule 11\n{table}"
    );
}

#[test]
fn the_split_reads_what_it_should() {
    let doc = "# title\nrules\n# 층과 하위 체계\nbody\n## 가지 말 것 — rows\n| a |\n";
    assert_eq!(slice(doc, None, Some(LAYERS)), "# title\nrules\n");
    assert_eq!(
        slice(doc, Some(LAYERS), Some(crate::history_words::STORY_SECTION)),
        "# 층과 하위 체계\nbody\n"
    );
    assert_eq!(
        slice(doc, Some(crate::history_words::STORY_SECTION), None),
        "## 가지 말 것 — rows\n| a |\n"
    );
    assert_eq!(slice(doc, None, None), doc);
    let todo = "# todo\n## 다음\nnext\n## 보류 — 인구가 생기면\n### item\n";
    assert_eq!(slice(todo, None, Some(PARKED)), "# todo\n## 다음\nnext\n");
    assert_eq!(
        slice(todo, Some(PARKED), None),
        "## 보류 — 인구가 생기면\n### item\n"
    );
    assert_eq!(heading_at("# title\n", LAYERS), Err(0), "a missing heading");
    assert_eq!(
        heading_at("# 층과 하위 체계\n# 층과 하위 체계 again\n", LAYERS),
        Err(2),
        "a heading on two lines"
    );
    assert_eq!(
        heading_at("text # 층과 하위 체계\n# 층과 하위 체계\n", LAYERS),
        Ok("text # 층과 하위 체계\n".len()),
        "a heading counts only at a line's start"
    );
}
