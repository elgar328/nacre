# Releasing

The steps for publishing nacre to crates.io. The development rules live in `overview.md`.

## What is released

- **The facade `nacre` and its normal-dependency closure**, all at the same version. The
  layers are mapped in `design.md` 「크레이트 구조」 (crate structure). `nacre-oracle` is
  marked `publish = false`, so `cargo publish --workspace` skips it; any new dev-only crate
  must be marked the same way.
- **Only `nacre`'s API is promised.** The `nacre-*` crates are internal layers and may be
  merged or split in any release. Users import only the facade's modules and its `prelude`.
- **The playground does not use releases.** It depends on nacre-kit and nacre by path and
  builds against their `main` branches, so `main` must always build.

## Versions

- Development happens directly on `main`; there is no separate `dev` branch to squash-merge.
  Commit bodies are the record of each change, and the playground builds from `main`.
- **While the version is `0.0.z`, every release bumps `z`.** Cargo treats `^0.0.z` as an exact
  version, so any release may break anything — which reflects where the project is today.
- **Between releases, `main` carries the next version with a `-dev` suffix** (after `0.0.1`,
  `0.0.2-dev`), so a build from `main` never claims to be a published version. Only release
  commits have a plain version.

### Where the version is written

**Only in the root `Cargo.toml`.**

- `[workspace.package] version` sets the version of every crate (each declares
  `version.workspace = true`).
- `[workspace.dependencies]` declares the dependencies between layers, each with a version
  requirement, and the crates use them via `{ workspace = true }`. Cargo needs this second
  copy: a published crate cannot depend on a path alone, and a version requirement cannot
  refer to the workspace version. As a result the version string appears once per layer, all
  within that one file.
- **Dev-dependencies between siblings have no version**: `{ path = "..." }` only, never
  `{ workspace = true }`. The comment above `[workspace.dependencies]` explains why: the
  tests form a dependency cycle, and versioned dev-dependencies would leave no valid publish
  order. `nacre-oracle` is never published and depends on its siblings by path only.

To change the version:

```sh
OLD=0.0.1-dev NEW=0.0.1
sed -i '' "s/\"${OLD//./\\.}\"/\"$NEW\"/g" Cargo.toml
grep -n "\"$OLD\"" Cargo.toml crates/*/Cargo.toml   # prints nothing
cargo check --workspace                                # updates Cargo.lock
```

## CHANGELOG

`CHANGELOG.md` at the root follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
It contains release notes for people using the `nacre` facade; it is not a development log.
The reasoning behind each change is in the commit bodies.

- Write the entries at release time from `git log vA.B.C..`. Include anything a `nacre` user
  can observe — a new operation, a refusal that now succeeds, a changed result — regardless of
  which layer the change is in. Leave out tests, instruments, docs and refactors with no
  visible effect.
- Write from the user's point of view: one concise, complete sentence per entry, **on a single
  line**, because GitHub release notes render every line break. Start an entry with
  **Breaking:** when existing code has to change, and explain how.
- Use these categories in this order, omitting empty ones: **Added** · **Changed** ·
  **Deprecated** · **Removed** · **Fixed** · **Security**.
- List the most recent release first, and always keep `## [Unreleased]` at the top.

## Before you start

- A Cargo version that supports `cargo publish --workspace` (check with `cargo publish --help`).
- `cargo login` with a crates.io token that has the `publish-new` and `publish-update` scopes
  (plus `yank`, for handling a broken release).
- The GitHub CLI, logged in with push rights (`gh auth status`).

## Release checklist

`X.Y.Z` is the new version and `vA.B.C` the previous release tag. Steps 1–5 are local and can
be undone. **From step 6 on, everything is public and permanent**: a pushed tag is never
moved, and a version number published to crates.io can never be reused.

1. **Start from an up-to-date `main`.**
   ```sh
   git switch main && git pull --ff-only
   git status --short                 # prints nothing
   git log --oneline origin/main..    # prints nothing
   ```

2. **Run every gate** listed in `overview.md` 「관문」 (gates), including the `cargo doc`
   baseline. If nacre-kit is being released as well, run its gates too.

3. **Finalize the version and the CHANGELOG.**
   - Change `X.Y.Z-dev` to `X.Y.Z` as described in
     [Where the version is written](#where-the-version-is-written).
   - Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD` (today's date) and add a new, empty
     `## [Unreleased]` above it. Update the links at the bottom:
     ```markdown
     [Unreleased]: https://github.com/elgar328/nacre/compare/vX.Y.Z...HEAD
     [X.Y.Z]: https://github.com/elgar328/nacre/compare/vA.B.C...vX.Y.Z
     ```

4. **Do a dry run of the publish.**
   ```sh
   cargo publish --workspace --dry-run --allow-dirty
   ```
   Each crate is packaged and then built from its packaged sources, and the output for each
   ends with `aborting upload due to dry run`. If anything fails, fix it before going on.

5. **Commit and tag.**
   ```sh
   git commit -am "Release X.Y.Z"
   git tag -a vX.Y.Z -m "nacre X.Y.Z"
   git describe --exact-match         # vX.Y.Z
   git show --stat HEAD               # only Cargo.toml, Cargo.lock and CHANGELOG.md
   ```
   To start over: `git tag -d vX.Y.Z; git reset --hard origin/main`.

6. **Push** the commit and the tag together, so that either both land or neither does:
   ```sh
   git push --atomic origin main vX.Y.Z
   ```

7. **Publish to crates.io** from the clean release commit (without `--allow-dirty`):
   ```sh
   cargo publish --workspace
   ```
   Cargo uploads the crates in dependency order. **If it stops partway**, the crates already
   uploaded remain published; publish the remaining ones with
   `cargo publish -p <crate> -p <crate> …`, which also orders them by dependency. A 429 error
   means crates.io's rate limit was hit; the message says when you can retry.

8. **Check the published crate.**
   ```sh
   T=$(mktemp -d) && cd $T && cargo new --lib probe && cd probe
   cargo add nacre@=X.Y.Z && cargo build
   cd / && rm -rf $T
   ```
   It can take about a minute for a new version to become available. Also check the docs.rs
   build (<https://docs.rs/crate/nacre/X.Y.Z>) and the README on the crates.io page.

9. **Create the GitHub Release** from the CHANGELOG section:
   ```sh
   gh release create vX.Y.Z --title "vX.Y.Z" --verify-tag --notes-file <(
     sed -n '/^## \[X\.Y\.Z\]/,/^## \[/{ /^## \[X\.Y\.Z\]/d; /^## \[/d; p; }' CHANGELOG.md
     echo "**Full Changelog**: https://github.com/elgar328/nacre/compare/vA.B.C...vX.Y.Z"
   )
   ```

10. **Release nacre-kit** if it is part of this release, following [nacre-kit](#nacre-kit).
    **Kit's next-cycle commit must land before step 11.**

11. **Start the next development cycle.** Change `X.Y.Z` to `X.Y.(Z+1)-dev` the same way, then:
    ```sh
    git commit -am "Start next development cycle"
    git push origin main
    ```

## A broken release

A version on crates.io cannot be replaced, and tags are never moved. Release the fix as the
next patch version and **yank** the broken version from every published crate. Cargo never
selects a yanked version for a new resolution, but projects whose `Cargo.lock` already lists
it still build.

```sh
for c in $(cargo metadata --no-deps --format-version 1 | python3 -c \
    'import sys,json; print(" ".join(p["name"] for p in json.load(sys.stdin)["packages"] if p["publish"] is None))'); do
  cargo yank --version X.Y.Z "$c"
done
```

The crate list comes from the manifests (every crate not marked `publish = false`), so it
stays correct as crates are added.

## nacre-kit

nacre-kit is a separate repository outside the workspace. It is released at **the same
version** as nacre: kit `X.Y.Z` depends on exactly nacre `X.Y.Z`. Between releases, kit's
version also carries `-dev`.

- **On kit's `main`, the nacre dependency is a path dependency only**, with no `version`. When
  a path dependency also specifies a version, Cargo checks that requirement against the path
  crate's actual version. As soon as nacre's `main` moves to `-dev`, the check fails, and the
  playground, which builds both `main` branches together, breaks.
- Publishing requires the version, so it appears **only in kit's release commit**. After
  nacre's step 7:
  1. Set kit's version to `X.Y.Z` and add `version = "X.Y.Z"` to its nacre dependency.
  2. Run kit's gates and `cargo publish --dry-run`. Commit `Release X.Y.Z`, tag `vX.Y.Z`,
     push with `git push --atomic origin main vX.Y.Z`, then run `cargo publish`.
  3. **Start kit's next cycle:** set its version to `X.Y.(Z+1)-dev`, remove `version` from the
     nacre dependency, then commit and push.
  4. Only then go on to nacre's step 11.
