# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Repository layout

Look in copilot-instructions.md in the .github directory. 

# Writing documentation

Look in copilot-instructions.md in the .github directory for instructions. 

# Binary 

It is in workspace_directory/../build.VACASK/Release/simulator/vacask

Make sure it is freshly built before running. 

If you run VACASK, make /tmp your working directory before running.

## Project docs discipline

This project uses the docs discipline from
<https://robtaylor.github.io/claude-project-discipline/>. Five kinds of
document, each with a clear home and a clear lifecycle:

| Kind | Lives in | Lifetime | Answers |
|---|---|---|---|
| **Architecture reference** | `docs/architecture/<area>.md` | Living; edited in place as the code moves | *What* does the code do now? |
| **Decision** | `docs/architecture/decisions/NNNN-*.md` | Forever; immutable; superseded ones move to `archive/` | *Why* did we choose this? |
| **Spike** | `docs/architecture/decisions/spikes/<topic>.md` | Forever, marked Resolved | *Did this idea work?* |
| **Plan** | `docs/plans/<topic>.md` | Updated as work lands; finished ones move to `plans/archive/` | *What's next, in what order?* |
| **Handoff** | `docs/handoffs/<topic>-handoff.md` | Ephemeral; folded into the others, then deleted | *What's in flight right now?* |

The load-bearing rule: **information has exactly one home, and the
handoff is the only doc that gets deleted.** Everything that survives a
session migrates from the handoff into an area doc / decision / spike /
plan / `CHANGELOG.md` *before* the handoff file is removed.

Start at `docs/architecture/README.md`: it maps each area to its
reference doc and the decisions behind it. VACASK's user manual is the
rest of `docs/` (the `cir-*`, `cmd-*`, `dev-*`, `expr-*`, `input-*`,
`startup-*` pages, indexed from `docs/index.md`); the discipline docs
above live alongside it.

### Where things go

| If you're about to write… | It belongs in… |
|---|---|
| "X does Y" (present tense, true today) | The area doc in `docs/architecture/` |
| "We chose X over Y because Z" | A new decision record (past tense) |
| "X is decided but not built yet" | The area doc's `## Implementation status`, plus a plan or issue |
| "What needs to happen next, in what order, with exit criteria" | The relevant plan doc |
| "I want to validate <assumption> before accepting Decision NNNN" | A new spike, time-boxed |
| "Here's what's in flight on thread T and what the next session should pick up" | `docs/handoffs/<t>-handoff.md` |
| "Here's what shipped in this version" | `CHANGELOG.md` `[Unreleased]`, led by "What this means for you" |

### Decisions are immutable

A decision record is written in past tense, as the state at decision
time, and its body is never edited after acceptance. Only the
`**Status:**` line and the `**Current architecture:**` forward-link
change. When the code moves, edit the area doc. When a decision changes,
write a new record that supersedes it and archive the old one. Don't add
dated `## Amendment` blocks.

Status is a claim about the code: verify it against the implementation
before setting it.

### Tooling

`docs/scripts/doctool.py` (stdlib only, run with `uv run`):

- `uv run docs/scripts/doctool.py lint` checks links, the decision index,
  archive placement, forward-links, spike/plan statuses, and handoff
  count. Run it before committing doc changes.
- `uv run docs/scripts/doctool.py archive <file.md>` moves a finished plan
  or superseded decision to `archive/` and repoints every reference to it
  (and the moved file's own links). Never move these files by hand.
- `uv run docs/scripts/doctool.py new decision "<title>"` creates the next
  numbered record and its index row. Also `new spike|plan|handoff|area`.

### Smell tests

If you notice these, the discipline has been broken. Fix it before
proceeding:

- **More live handoffs than ~2× the people actively working** → too many
  threads in flight. Resolve, consolidate, or drop some; don't open more.
- **A decision record describes what the code does today** → that's the
  area doc's job. Move the present-tense claim there.
- **A decision record contains "first do X, then do Y"** → that's a plan.
- **A handoff contains "we chose X over Y because Z"** → that's a
  decision. Write the record.
- **A plan says Implemented but sits in `docs/plans/`** → archive it.
- **A plan hasn't been updated in months but the work has evolved** →
  update it in place; past states live in `git log`.

### Handoff write/read protocol

Keep **one handoff per active thread of work**, not one globally.

**At session end (work in progress, partial state):**

1. **Decision sweep first.** Before writing the handoff, scan this
   session's diffs and discussion for any architectural decision not yet
   recorded. Write the record *now*, while the rationale is fresh. If it
   changes an existing decision, write a superseding record. Then check
   that the area docs describe what the code now does.
2. **Then write the handoff** at `docs/handoffs/<topic>-handoff.md`
   using `docs/handoffs/_template.md`: goal, next-up,
   done-this-session, open follow-ups, critical context, verification
   command. Skip it if the session ended at a clean stopping point.

**At session start (resuming work):** read the handoff for the thread
you're resuming. Run its verification command first to confirm the
claimed state matches reality, then proceed with its next-up item.

**At resolution:** migrate every load-bearing piece of the handoff into
its proper home, then `git rm docs/handoffs/<topic>-handoff.md` in the
same commit as the migration. Commit message: `docs: resolve <topic>
handoff — fold into <where-it-went>`.

The migration table is in `docs/handoff-discipline.md`.

### Override for skill toolkit defaults

Various Claude Code skill toolkits (`create_handoff`, `resume_handoff`,
etc.) default to YAML format under `thoughts/shared/handoffs/` with
database indexing. **That doesn't apply here.** Override:

- Format: markdown.
- Location: `docs/handoffs/<topic>-handoff.md`.
- No database indexing.

The skill activation is informational; this project's convention takes
precedence.

### When NOT to write a handoff

- Session ended at a clean stopping point (everything merged, decisions
  recorded, area docs current, nothing surprising left over).
- Work was a single small commit with no follow-ups.
- The plan doc already captures what's next at sufficient detail.

### When NOT to write a decision

- Mechanical changes (renames, typo fixes, dependency bumps).
- Choices the code itself makes self-evident.
- Things adequately captured by a commit message.

Test: *if someone asked "why is it like this?" six months from now,
would the answer be in the code, the commit, or a decision record?* Only
write the record if the answer is "a decision record".

### Writing docs

Prose in the docs follows `docs/prose-style.md`: direct voice, no
AI-marker words, no mid-sentence interjections, and present-tense claims
that track the code. Read it before writing or editing a doc.
