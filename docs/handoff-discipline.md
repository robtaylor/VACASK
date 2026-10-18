# Handoff discipline

Handoffs in this project are **ephemeral working memory**, not historical record. They exist to bridge a session boundary (when you stop working and someone else, Claude or human, picks up) and they are deleted once the work they describe is resolved.

This document defines what a handoff is, what it isn't, when to write one, and exactly what to do when one is resolved.

For the broader discipline (architecture reference, decisions, spikes, plans, handoffs), see <https://robtaylor.github.io/claude-project-discipline/>.

## Why this discipline exists

Current state, rationale, and project state all have natural homes:

- **Area docs** (`docs/architecture/<area>.md`) capture *what* the code does now.
- **Decision records** (`docs/architecture/decisions/`) capture *why* it got that way.
- **Plan docs** (`docs/plans/<topic>.md`) capture *what's left* and the next workstream slices.

When that content lives in a handoff instead, two things go wrong:

1. **It's not where contributors look.** A new contributor reading the README → architecture map → decisions chain shouldn't have to dig through a stack of resolved handoff docs to find load-bearing decisions or the current state of a workstream.
2. **It rots out of sync with reality.** Handoffs are point-in-time snapshots. A "STATUS: RESOLVED" banner doesn't help when the thing referenced has moved or changed; the canonical doc is what should hold the current truth.

The discipline closes this gap by **forcing migration before deletion**. Every load-bearing piece of a handoff lands in its proper home (area doc / decision / plan) before the handoff file is removed.

## What a handoff IS

A handoff *lives* in its own dedicated directory, separate from the persistent plan docs whose content it eventually feeds: a markdown file at `docs/handoffs/<topic>-handoff.md` containing exactly what the next session needs to pick up where you left off:

- **Goal & next-up** — what this session was trying to do, and what the *very next* concrete action is.
- **Done this session** — commits landed, with one-line summaries.
- **Open follow-ups** — the work that wasn't done, with enough scope detail to start cold.
- **Critical context** — gotchas, surprising findings, environment specifics that aren't obvious from the code or docs *yet*.
- **Verification** — the command(s) the next session runs to confirm the work is in the state you say it is.

**One handoff per active thread of work**, not one globally. A *thread* is a distinct arc someone could pick up cold (a migration, a feature, a triage sweep). Concurrent threads get concurrent handoffs, each named by topic. There's still no *chain*: a thread holds at most one live handoff, and a resolved one is folded-and-deleted (see below), not archived behind its successor.

**Treat the handoff count as a WIP signal.** If the number of live handoffs climbs past roughly **2× the number of people actively working**, the team has more parallel threads than it can hold context on. Read it as a prompt to *resolve, consolidate, or drop* threads, not to open more. `uv run docs/scripts/doctool.py lint --people N` checks this.

## What a handoff IS NOT

- **Not a decision log.** If you find yourself writing "we chose X over Y because Z" in a handoff, that paragraph belongs in a decision record.
- **Not a design doc.** "How subsystem A works" belongs in subsystem A's area doc, not in a handoff's "Critical context" section.
- **Not a status dashboard for the project.** Workstream status lives in plan docs. A handoff cites those, doesn't reproduce them.
- **Not a historical record.** `git log` is the historical record. Handoffs that survive past their resolution turn into noise that misleads new contributors.

## When to write one

Write a handoff at the end of any session that:

1. **Leaves work in a partial state** that someone else might pick up cold.
2. **Captures non-obvious context** the next session needs.
3. **Documents the next concrete step** with enough scope to start without re-discovering it.

If the session ended at a clean stopping point (everything merged, all decisions recorded, area docs current, nothing surprising), don't write a handoff. The plan doc already says what's next.

## Before you hand off: sweep for undocumented decisions

A handoff is written at the moment your session context is richest and about to evaporate. That makes it the right checkpoint for one deliberate question:

> **Did I make any architectural decision this session that isn't yet captured in a decision record, or change the code in a way the area docs don't yet describe?**

Run the sweep *before* writing the handoff, not after. Walk the session's diffs and discussion for choices that pass the decision test: *"if someone asked 'why is it like this?' six months from now, would the answer be in the code, the commit, or a decision record?"* For each one whose answer is "a decision record":

- **No decision covers it** → write a new one now, while the rationale and the alternatives you weighed are still in working memory.
- **It changes an existing decision** → write a new record that supersedes the old one, and archive the old one with `uv run docs/scripts/doctool.py archive <old-record>`. Don't edit the old record's body.

Then check the area docs. If the session changed what the code does, edit the area doc now, so its present-tense claims stay true.

Only *then* write the handoff. This keeps the decision out of the handoff's "Critical context", where it would be lost on fold-and-delete, and puts it in its permanent home at the moment you can write it best.

Two failure modes this closes:

1. **Decisions captured nowhere.** Without the sweep, a choice made mid-session lives only in the chat transcript and the diff. The transcript is gone next session; the diff shows *what* changed, never *why this and not that*.
2. **Decisions smuggled into the handoff.** The shortcut is to jot "we went with X because Y" into the handoff and move on. But the handoff is deleted at resolution, and if fold-and-delete misses that paragraph the rationale dies with the file. Writing the decision up front means there's nothing to migrate later.

## Resolution: fold, then delete

The two-location split is deliberate: handoffs *live* at `docs/handoffs/<topic>-handoff.md` while in flight; their *content* migrates into the persistent docs (`docs/architecture/`, `docs/plans/`) at resolution. The handoff file then gets removed; nothing about the work is lost because everything load-bearing has a permanent home elsewhere.

When a handoff's work is done (whether in the next session or several sessions later), every load-bearing piece of it must be migrated to its proper home **before the handoff file is deleted**:

| If the handoff says... | It belongs in... |
|---|---|
| "We chose approach X over Y because Z" | A new decision record (superseding an old one if it changes it) |
| "Subsystem A now does B" | Subsystem A's area doc |
| "Subsystem A doesn't yet do B" | The area doc's Implementation status, plus a plan item or issue |
| "Future scope for WS-N: do A then B then C" | The plan doc's WS-N section |
| "Gotcha: tool Z's API X behaves Y" | A code comment near the call site, or the area doc if the gotcha cuts across files |
| "Build dep Z is required on Linux" | The build script's apt-suggestion / Brewfile / README install section |
| "Run `<command>` to verify" | The verification block in the relevant plan doc, or a test-running section in `CLAUDE.md` |

After migration, the handoff file is removed in the same commit as the migration:

```sh
git rm docs/handoffs/<topic>-handoff.md
git add <files-receiving-the-migrated-content>
git commit -m "docs: resolve <topic> handoff — fold into <where-it-went>"
```

The commit message records what migrated where; that's the audit trail. `git log -- docs/handoffs/` then shows the project's handoff history (one add, one delete per thread) without needing the files themselves to live forever.

## Template

See `docs/handoffs/_template.md` for the skeleton, or copy from <https://github.com/robtaylor/claude-project-discipline/blob/main/starter/docs/handoffs/_template.md>.

## Tooling

The `create_handoff` and `resume_handoff` skills (from various Claude Code orchestration toolkits) generate and consume handoffs. They're optional; the discipline above is the load-bearing artefact. A handoff written by hand following the template is just as valid.

If you use one of those skills, expect them to default to YAML format under `thoughts/shared/handoffs/` with database indexing. **That doesn't apply to this project.** Override it: produce markdown at `docs/handoffs/<topic>-handoff.md` and skip the database step. The skill activation is informational; the project's convention takes precedence. (The CLAUDE.md snippet from the upstream discipline includes this override.)
