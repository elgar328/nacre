use super::*;
/// **What a face leaves on one line, read off the stations where its boundary meets the line.**
///
/// `stations` are sorted along the line, and at one point by the order their events apply
/// ([`LineEvent`]'s `Ord`); `same[i]` says stations `i` and `i + 1` are one point. Each station
/// carries the caller's payload `P`, of which `name` reads the point's name. The face is absent
/// before the first station; the gap after each station is a `Graze` while a run is open, `mat`
/// while the face is present, and nothing otherwise; adjacent gaps of one kind come back as one
/// piece. The road that collected the stations knows its line — a class line on a plane face, a
/// ruling on a lateral — and this does not.
///
/// ★ **A station is a point, not an event.** A boundary that passes one point twice — a face
/// pinched there — leaves two events under one name at one station: a run's end and a touch, or
/// two runs meeting end to start. They are one point and are read as one station; refusing them
/// as two coinciding features refuses every pinched face. Two **different** names at one station
/// are refused (`CoincidentFeatures`) before the sweep: one point under two names is an identity
/// the alias table did not learn, and this is no place to decide it.
///
/// Refused by name, besides: a run's end with no run open, a crossing or a second run inside an
/// open run (`RunSplit` — the run is not one interval of the boundary), and a face still present
/// or a run still open past the last station (`OddParity` — a closed boundary leaves the line as
/// often as it enters).
pub(super) fn sweep_line<P: Copy>(
    stations: &[(LineEvent, P)],
    same: &[bool],
    name: impl Fn(P) -> NodeId,
    mat: SegKind,
) -> Result<Vec<([P; 2], SegKind)>, DeclineKind> {
    debug_assert_eq!(same.len(), stations.len().saturating_sub(1));
    if (1..stations.len()).any(|i| same[i - 1] && name(stations[i - 1].1) != name(stations[i].1)) {
        return Err(DeclineKind::CoincidentFeatures);
    }
    let mut pieces: Vec<([P; 2], SegKind)> = Vec::new();
    let (mut present, mut graze): (bool, Option<bool>) = (false, None);
    // Whether the gap just before the current station left a piece, which the next gap of the same
    // kind then extends.
    let mut extends = false;
    let mut i = 0;
    while i < stations.len() {
        let at = stations[i].1;
        let mut j = i;
        loop {
            match stations[j].0 {
                LineEvent::GrazeEnd { toggle } => {
                    if graze.take().is_none() {
                        return Err(DeclineKind::RunSplit);
                    }
                    present ^= toggle;
                }
                LineEvent::Toggle => {
                    if graze.is_some() {
                        return Err(DeclineKind::RunSplit);
                    }
                    present = !present;
                }
                LineEvent::GrazeStart { body_above } => {
                    if graze.replace(body_above).is_some() {
                        return Err(DeclineKind::RunSplit);
                    }
                }
            }
            j += 1;
            if j == stations.len() || !same[j - 1] {
                break;
            }
        }
        if j < stations.len() {
            let next = stations[j].1;
            let kind = match graze {
                Some(body_above) => Some(SegKind::Graze { body_above }),
                None if present => Some(mat),
                None => None,
            };
            if let Some(kind) = kind {
                match pieces.last_mut() {
                    Some((end, k)) if extends && *k == kind => end[1] = next,
                    _ => pieces.push(([at, next], kind)),
                }
            }
            extends = kind.is_some();
        }
        i = j;
    }
    if present || graze.is_some() {
        return Err(DeclineKind::OddParity);
    }
    Ok(pieces)
}
