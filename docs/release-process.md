# Release process

<!--
This template is the lightweight release process companion to the Claude Project
Discipline (architecture reference, decisions, spikes, plans, handoffs). Adjust
the <PROJECT> placeholders, then prune the sections that don't apply.
-->

Lightweight by design. <PROJECT> is a <single-binary | library | service> project; releases are <git tags + a CHANGELOG entry | crates.io / npm / PyPI publications | container images>. <Optionally: no pre-built binaries until/unless that demand surfaces.>

## When to release

Cut a release when:

- A user-visible feature or fix lands that you want consumers to be able to pin against.
- Schema or Command-Line Interface (CLI) changes happened and consumers need a stable reference point.
- A meaningful chunk of work in `docs/plans/` is closed (e.g. a phase exits all criteria).

There is no fixed cadence.

## Versioning

[Semantic Versioning (SemVer)](https://semver.org/), starting once the first numbered release ships. Pre-1.0 versions (`0.x.0`) carry the standard SemVer caveat: minor bumps may include breaking changes; specific stable contracts are documented in their own decision records and follow stricter rules.

Stable contracts (additive-only; breaking changes require a major bump and a deprecation window):

- `<contract 1>` — `<file:line or doc reference>`, governed by Decision <NNNN>.
- `<contract 2>` — `<file:line or doc reference>`, governed by Decision <MMMM>.

CLI flags, log message formats, and human-readable text outputs are **not** stable parseable contracts; consumers that need to script against them should use the structured contracts above.

## Release notes

Write the notes in the `[Unreleased]` section of [`CHANGELOG.md`](../CHANGELOG.md) as work lands. Open them with a short **"What this means for you"** block: a few benefit-framed bullets saying what a user can now *do*, before the technical `### Added` / `### Changed` detail. That block becomes the release intro.

## Steps

For maintainers cutting a release:

1. **Verify CI is green** on `main` for all jobs. If a runner is offline, hold the release until it's restored.

2. **Roll the `[Unreleased]` section in [`CHANGELOG.md`](../CHANGELOG.md) into a numbered version block.** Format follows [Keep a Changelog](https://keepachangelog.com/). Update the link references at the bottom. Leave a fresh empty `[Unreleased]` section.

3. **Bump `version`** in your project manifest (`Cargo.toml` / `package.json` / `pyproject.toml` / etc.) and re-build to update the lockfile.

4. **Commit:** `chore: release v<X.Y.Z>` with whatever attribution conventions your project uses.

5. **Tag:** `git tag -a v<X.Y.Z> -m "v<X.Y.Z>"` then `git push --tags`.

6. **Create a release** on the platform of your choice from the tag. Body = the CHANGELOG section for that version.

## Release candidates (optional)

<Delete this section if <PROJECT> doesn't cut RCs.>

During the RC window, `main` **is** the candidate. It stays at `X.Y.Z-rc.N` with the notes still under `[Unreleased]`; each new RC is a `main` commit bumped to the next `-rc.N` and tagged `vX.Y.Z-rc.N`. Publish RCs as prereleases, never marked "Latest", and validate the documented install commands against them.

**Don't roll `main` to the final version before the RC validates.** Only promotion rolls `[Unreleased]` into a dated `[X.Y.Z]` block (step 2 above) and bumps to `X.Y.Z`. Rolling first leaves `main` claiming a version that has no release.

## What does NOT need to change at release time

- Submodule pins (unless deliberately bumping a vendored dep).
- Area docs in `docs/architecture/` (unless the release changes the design they describe).
- `LICENSE` (unless re-licensing).

## Pre-release checklist (one-time, before the first numbered release)

These items block the first release. As they land, they're checked off in the same Pull Request (PR) that landed them.

- [ ] <First load-bearing release blocker>
- [ ] <Second load-bearing release blocker>
- [ ] All Continuous Integration (CI) jobs green on `main`.
- [ ] License posture confirmed; `NOTICE` file in place if vendoring third-party code.
- [ ] End-to-end test exercising the headline feature.
- [ ] User-facing docs cover the install path.

## Cross-references

- [`CHANGELOG.md`](../CHANGELOG.md) — release log.
- [`docs/architecture/decisions/<NNNN>-…md`](architecture/decisions/<NNNN>-…md) — stability contract for `<contract 1>`.
