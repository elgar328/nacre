# An internal layer of nacre

This crate is one layer of [nacre](https://crates.io/crates/nacre), an exact b-rep CAD
kernel in pure Rust. It is published only so that the `nacre` facade can depend on it.

**Depend on [`nacre`](https://crates.io/crates/nacre) instead.** The facade re-exports
everything a user needs, and its API is the only one nacre promises: the layers beneath
it may be renamed, merged or split in any release.

Source and documentation: <https://github.com/elgar328/nacre>
