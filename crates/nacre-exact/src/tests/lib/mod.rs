use super::*;

use proptest::prelude::*;

// Small ranges keep checked arithmetic inside i128, so these exercise the
// algebraic laws, not the overflow path (which its own test above pins).
prop_compose! {
    fn small_rat()(num in -1000i128..=1000, den in 1i128..=1000) -> Rat {
        Rat::new(num, den).unwrap()
    }
}

mod angles_and_trig;
mod cylinder_strip;
mod decimals_and_wide;
mod planes;
