use core::cmp::Ordering;

/// **A stable sort by a comparison that can fail** — the first failure is returned at once,
/// before any order is formed from it.
///
/// **Not `sort_by` with a flag.** A comparison that reads «could not tell» as `Equal` is not a
/// total order, and the standard sort may panic on one past twenty elements — measured: 1,896 of
/// 2,000 sorts of forty, with one element in ten failing every comparison it takes part in. An
/// exact comparison here can fail (a cylinder-pinned point against a class with no world
/// statement), so the sort has to stop rather than carry on with a guess.
///
/// On a comparison that never fails the result is the standard stable sort's — a stable order
/// is unique — so ties a caller reads (one point under two names, kept in the order it gave them)
/// come out as they did.
///
/// **On `Err` the slice is a permutation of its input**: a merge is written to scratch and
/// copied back only whole, and an insertion moves by swaps. A caller may read the slice to name
/// the failure (`split_at_crossings` asks which classes it holds).
///
/// Up to twenty elements this is the standard's own insertion sort, asking `cmp(tail, prev)` in
/// the same order, so a short sort asks what it asked before; past that, a top-down merge that
/// skips a merge whose halves already join in order.
pub(crate) fn try_sort_by<T: Copy, E, F>(v: &mut [T], mut cmp: F) -> Result<(), E>
where
    F: FnMut(&T, &T) -> Result<Ordering, E>,
{
    let mut scratch = Vec::new();
    merge_sort(v, &mut cmp, &mut scratch)
}

/// The length the standard sort still sorts by insertion.
const INSERTION: usize = 20;

fn merge_sort<T: Copy, E, F>(v: &mut [T], cmp: &mut F, scratch: &mut Vec<T>) -> Result<(), E>
where
    F: FnMut(&T, &T) -> Result<Ordering, E>,
{
    let n = v.len();
    if n <= INSERTION {
        for i in 1..n {
            let mut j = i;
            while j > 0 && cmp(&v[j], &v[j - 1])? == Ordering::Less {
                v.swap(j, j - 1);
                j -= 1;
            }
        }
        return Ok(());
    }
    let mid = n / 2;
    merge_sort(&mut v[..mid], cmp, scratch)?;
    merge_sort(&mut v[mid..], cmp, scratch)?;
    if cmp(&v[mid], &v[mid - 1])? != Ordering::Less {
        return Ok(());
    }
    scratch.clear();
    let (mut i, mut j) = (0, mid);
    while i < mid && j < n {
        // Strictly less takes the right one: an equal pair keeps the left first, which is what
        // makes the sort stable.
        if cmp(&v[j], &v[i])? == Ordering::Less {
            scratch.push(v[j]);
            j += 1;
        } else {
            scratch.push(v[i]);
            i += 1;
        }
    }
    scratch.extend_from_slice(&v[i..mid]);
    scratch.extend_from_slice(&v[j..]);
    v.copy_from_slice(scratch);
    Ok(())
}
