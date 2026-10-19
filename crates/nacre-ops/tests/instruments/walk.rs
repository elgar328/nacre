//! **One tree walk** for the checks that read the repository by path.
//!
//! Each scan names its own roots, the folders it will not enter and the files it reads; the walk
//! underneath is this one, so two scans cannot come to walk differently without one of them
//! saying so in its own arguments.

use std::path::{Path, PathBuf};

/// Every file under `root`, descending into a folder only when `enter` accepts its name.
///
/// A folder is what [`Path::is_dir`] says, which follows symbolic links: the crates' `LICENSE-*`
/// files are links, and a walk that read the entry's own file type would drop them. `enter` is
/// asked at every depth, so a rule like «not into `tests`» holds for `src/polygon/tests` as well
/// as `src/tests`. The order is the file system's; each caller sorts by its own key, because the
/// scans that print relative paths sort them as strings and a path sorts by component, and the two
/// orders part wherever `x.rs` sits beside `x/`.
pub fn files_under(root: &Path, enter: impl Fn(&str) -> bool) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut dirs = vec![root.to_path_buf()];
    while let Some(dir) = dirs.pop() {
        for entry in std::fs::read_dir(&dir).expect("a readable folder") {
            let path = entry.expect("a readable entry").path();
            if path.is_dir() {
                if path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(&enter)
                {
                    dirs.push(path);
                }
            } else {
                out.push(path);
            }
        }
    }
    out
}
