# Nacre

**Nacre is an exact b-rep CAD kernel.** Every topological decision is exact — it never fails silently.

- **No global tolerance.** Dimensions and angles are kept as exact rationals, and positions are defined as intersections of surfaces, not stored as coordinates. A coordinate, when needed, is a cache — an f64 value with its own local tolerance, recomputable at any higher precision on demand. Tolerance is thus local to each site and measured — there is no global constant to tune.

- **Rounding cannot corrupt topology.** Each point is defined by the surfaces meeting at it, and indirect predicates decide on that definition. The definition yields a tolerance that says how much precision the sign needs — f64 where that settles it, higher where it doesn't. The decision is always exact; it just pays for precision only where a site demands it.

- **Never silently guesses.** When a decision stays ambiguous even at full precision, it is never resolved by picking — interactively the kernel asks, in batch a policy decides. And anything outside the supported coverage is a named error, not a plausible-looking wrong solid.

- **The log is the truth.** Objects are immutable once written and referenced by handle, so identity is an integer comparison, never a coordinate one. A model is the bit-exact replay of its operation log — coordinates are only a cache, regenerable at any precision.

> Named for *nacre*, mother-of-pearl: it grows one layer at a time and never rewrites a layer beneath. The kernel treats geometry the same way, and unusually, its topology too.

## Relation to Fornjot

Nacre began after [Fornjot](https://github.com/hannobraun/fornjot) — Hanno Braun's Rust b-rep CAD kernel, developed from 2020 until it was shut down in 2026 — and shares much of its outlook: code-first mechanical CAD, immutable objects referenced by handle, and clear errors in place of quietly wrong results. No code is shared.

The two differ on what counts as the truth. Fornjot's [final experiment](https://github.com/hannobraun/fornjot/tree/main/experiments/2025-12-03) made approximated geometry the uniform representation, with topology recorded alongside as the geometry is built — an explicit trade of exactness for simplicity. Nacre makes the opposite trade: exact geometry is the truth, and approximations are a cache derived from it.
