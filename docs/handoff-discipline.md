# Handoff discipline

Handoffs in this project are **ephemeral working memory**, not historical record. They exist to bridge a single session boundary -- when you stop working and someone else (Claude or human) picks up -- and they are deleted once the work they describe is resolved.

This document defines what a handoff is, what it isn't, when to write one, and exactly what to do when one is resolved.

For the broader four-document discipline (Architecture Decision Records, plans, spikes, handoffs), see <https://robtaylor.github.io/claude-project-discipline/>.

## Why this discipline exists

Decision rationale, technical context, and project state all have natural homes:

- **Architecture Decision Records** (ADRs) in `docs/adr/` capture architectural decisions and their *why*.
- **Design docs** (`docs/<topic>.md`) capture *how* things work.
- **Plan docs** (`docs/plans/<topic>.md`) capture *what's left* and the next workstream slices.

When that content lives in a handoff instead, two things go wrong:

1. **It's not where contributors look.** A new contributor reading the README -> ADR chain shouldn't have to dig through a stack of resolved handoff docs to find load-bearing decisions or the current state of a workstream.
2. **It rots out of sync with reality.** Handoffs are point-in-time snapshots. A "STATUS: RESOLVED" banner doesn't help when the thing referenced has moved or changed; the canonical doc is what should hold the current truth.

The discipline closes this gap by **forcing migration before deletion**. Every load-bearing piece of a handoff lands in its proper home (ADR / plan / design doc) before the handoff file is removed.

## What a handoff IS

A handoff *lives* in its own dedicated directory, separate from the persistent plan docs whose content it eventually feeds: a single markdown file at `docs/handoffs/<topic>-handoff.md` containing exactly what the next session needs to pick up where you left off:

- **Goal & next-up** -- what this session was trying to do, and what the *very next* concrete action is.
- **Done this session** -- commits landed, with one-line summaries.
- **Open follow-ups** -- the work that wasn't done, with enough scope detail to start cold.
- **Critical context** -- gotchas, surprising findings, environment specifics that aren't obvious from the code or docs *yet*.
- **Verification** -- the command(s) the next session runs to confirm the work is in the state you say it is.

Exactly one handoff exists at a time. There's no chain of resolved predecessors to wade through.

## What a handoff IS NOT

- **Not a decision log.** Decisions go in ADRs. If you find yourself writing "we chose X over Y because Z" in a handoff, that paragraph belongs in an ADR (or an existing ADR's "Consequences" / "Walk-back" section).
- **Not a design doc.** "How clock arrival flows from <upstream tool> through the IR into the GPU constraint buffer" is a design topic; it lives in `docs/<topic>.md`, not in a handoff's "Critical context" section.
- **Not a status dashboard for the project.** Workstream status lives in plan docs. A handoff cites those, doesn't reproduce them.
- **Not a historical record.** `git log` is the historical record. Handoffs that survive past their resolution turn into noise that misleads new contributors.

## When to write one

Write a handoff at the end of any session that:

1. **Leaves work in a partial state** that someone else might pick up cold.
2. **Captures non-obvious context** the next session needs.
3. **Documents the next concrete step** with enough scope to start without re-discovering it.

If the session ended at a clean stopping point (everything merged, all decisions documented in ADRs/plans, nothing surprising), don't write a handoff. The plan doc already says what's next.

## Resolution: fold, then delete

The two-location split is deliberate: handoffs *live* at `docs/handoffs/<topic>-handoff.md` while in flight; their *content* migrates into the persistent docs (`docs/adr/`, `docs/plans/`, design docs under `docs/`) at resolution. The handoff file then gets removed; nothing about the work is lost because everything load-bearing has a permanent home elsewhere.

When a handoff's work is done -- whether in the next session or several sessions later -- every load-bearing piece of it must be migrated to its proper home **before the handoff file is deleted**:

| If the handoff says... | It belongs in... |
|---|---|
| "We chose approach X over Y because Z" | The relevant ADR's Decision/Consequences section, or a new ADR if no fit exists |
| "Future scope for WS-N: do A then B then C" | The plan doc's WS-N section |
| "Gotcha: tool Z's API X behaves Y" | A code comment near the call site, or a design doc if the gotcha cuts across files |
| "Build dep Z is required on Linux" | The build script's apt-suggestion / Brewfile / README install section |
| "Subsystem A doesn't yet do B" | Plan doc as a new open item, or an ADR-tracked walk-back if it's a deferred design choice |
| "Run `<command>` to verify" | The verification block in the relevant plan doc, or a test-running section in `CLAUDE.md` |

After migration, the handoff file is removed in the same commit as the migration:

```sh
git rm docs/handoffs/<topic>-handoff.md
git add <files-receiving-the-migrated-content>
git commit -m "docs: resolve <topic> handoff -- fold into <where-it-went>"
```

The commit message records what migrated where -- that's the audit trail. `git log -- docs/handoffs/` then shows the project's handoff history (one add, one delete per session) without needing the files themselves to live forever.

## Template

See `docs/handoffs/_template.md` for the skeleton, or copy from <https://github.com/robtaylor/claude-project-discipline/blob/main/starter/docs/handoffs/_template.md>.
