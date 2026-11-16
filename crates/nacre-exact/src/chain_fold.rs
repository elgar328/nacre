//! **What a motion chain folds into, and why a statement it carries can have no rational image.**

use super::*;

/// **Why a statement carried along a [`ChainFold`] has no image in `Rat`.**
///
/// The two are asked in this order — whether the image is rational at all, then whether `Rat`
/// holds it — so a statement the chain puts off the rationals is never named for its width.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Carry {
    /// The chain carries this statement to an irrational image: no width would state it.
    Irrational,
    /// The image is rational and `Rat` cannot hold it, or a step on the way overflowed.
    Width,
}

/// **A motion chain, folded** — the map a statement written against the chain's leaf is carried
/// to the world by, asked one statement at a time.
///
/// Every carry answers for the statement it is given: its exact image, or why there is none
/// ([`Carry`]). A direction is asked two ways because two questions are asked of it: as a vector
/// ([`ChainFold::dir_rat`] — its length kept, for a caller that measures along it) and as a line
/// ([`ChainFold::dir_along`] — some positive multiple of its image, for a caller that only needs
/// where it points).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainFold {
    axis: AxisAffine,
}

impl ChainFold {
    /// The fold of a chain of exact nodes — translations, quarter turns and axis mirrors.
    pub fn axis(axis: AxisAffine) -> ChainFold {
        ChainFold { axis }
    }

    /// This fold with the exact node `outer` applied after it.
    pub fn then_axis(&self, outer: &AxisAffine) -> ChainFold {
        ChainFold {
            axis: self.axis.then(outer),
        }
    }

    /// `+1` for a proper map, `−1` for one carrying an odd number of reflections.
    pub fn det(&self) -> i8 {
        self.axis.det()
    }

    /// The image of the point `p`.
    pub fn point_rat(&self, p: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        self.axis.point_rat(p).ok_or(Carry::Width)
    }

    /// The image of the vector `d`, its length kept.
    pub fn dir_rat(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        self.axis.dir_rat(d).ok_or(Carry::Width)
    }

    /// A positive multiple of the image of the direction `d` — where it points, not how long it
    /// is.
    pub fn dir_along(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        self.dir_rat(d)
    }

    /// The image of the plane `n·x + d = 0`, canonicalized ([`canonical_plane_coeffs`]).
    pub fn plane_coeffs(&self, c: [Rat; 4]) -> Result<[Rat; 4], Carry> {
        self.axis.plane_coeffs(c).ok_or(Carry::Width)
    }

    /// The point whose image is `p`.
    pub fn point_rat_inverse(&self, p: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        self.axis.inverse().point_rat(p).ok_or(Carry::Width)
    }

    /// A positive multiple of the direction whose image is `d`.
    pub fn dir_along_inverse(&self, d: [Rat; 3]) -> Result<[Rat; 3], Carry> {
        self.axis.inverse().dir_rat(d).ok_or(Carry::Width)
    }

    /// The image of an `f64` vector — exact, since the map only permutes and negates.
    pub fn dir_f64(&self, v: [f64; 3]) -> Option<[f64; 3]> {
        Some(self.axis.dir_f64(v))
    }
}
