//! **One line lexer** for the checks that read Rust sources as text.
//!
//! A scan that looks for a name in code must not count it inside a string or a comment, and a scan
//! that reads comments must not count a string's `//` as one; both read a line through
//! [`split_line`], so the two cannot come to disagree about where a comment starts.
//!
//! It reads one line at a time: a string that runs on past a line break has its continuation read
//! as code, and a `//` inside that continuation (a URL, say) as a comment. Block comments are not
//! read at all. The tree has neither.

/// `t` split at its comment: the code with every string literal emptied (the quotes kept), and the
/// comment from the first `//` outside a string, if there is one — so a lint named in a message or
/// a comment is not code, and a `//` inside a string is not a comment. A char literal (`'"'`,
/// `'\"'`) is kept whole, so its quote does not open a string that would swallow the rest of the
/// line.
pub fn split_line(t: &str) -> (String, Option<String>) {
    let chars: Vec<char> = t.chars().collect();
    let mut code = String::new();
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
                    code.push('"');
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
                code.push(c);
                code.extend(&chars[k..k + len]);
                k += len;
            }
            '"' => {
                in_str = true;
                code.push('"');
            }
            '/' if chars.get(k) == Some(&'/') => {
                return (code, Some(chars[k - 1..].iter().collect()));
            }
            _ => code.push(c),
        }
    }
    (code, None)
}
