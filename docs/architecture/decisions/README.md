# Decisions

A decision record captures a choice worth understanding later: the context, the options considered, and the rationale for the one taken. It's the *why* behind a piece of the architecture. The *what the code does now* lives in the [architecture reference](../README.md), one doc per area. Each decision links forward to the area doc it shaped, and each area doc links back to the decisions behind it.

For the full convention, see [Claude Project Discipline — Decision records](https://robtaylor.github.io/claude-project-discipline/decisions.html).

## The rule: immutable and past-tense

A decision record is a record of a choice at the moment it was made, not a description of the system now. Its body is **never rewritten** when the code moves; edit the area doc instead. Only two lines change after acceptance: the `**Status:**` line and the `**Current architecture:**` forward-link.

When a decision itself changes, **write a new record that supersedes the old one**. Set the old one's status to `Superseded by Decision NNNN (YYYY-MM-DD)` and move it to [`archive/`](archive/):

```sh
uv run docs/scripts/doctool.py archive docs/architecture/decisions/NNNN-<name>.md
```

No dated `## Amendment` blocks. A change to *what the code does* is an area-doc edit; a change to *why* is a new decision.

## Three kinds of claim, aged differently

1. **Why the decision was made.** Past tense, as the state at decision time. Lives here, never rewritten.
2. **What the current architecture is.** Present tense, kept in step with the code. Lives in the area doc.
3. **What's decided but not built yet.** Lives in the area doc's `## Implementation status`, with a plan or issue, and drains into (2) as each piece lands.

Write all three in the same present tense and a reader can't tell *why we did it* from *what the code does* from *what we still intend*.

## Status legend

- **Proposed** — drafted, not yet ratified. Open the body with what the code does today, so the proposal can't be mistaken for reality.
- **Accepted (YYYY-MM-DD)** — the decision is in effect in the code.
- **Accepted (partial) (YYYY-MM-DD)** — ratified and partly built; the area doc's Implementation status says which parts.
- **Superseded by Decision NNNN (YYYY-MM-DD)** — replaced; lives in `archive/`.
- **Withdrawn (YYYY-MM-DD)** — abandoned without replacement; lives in `archive/`.

**Status must match the code.** A Status line is a claim about the codebase. Verify it against the implementation before setting it; don't trust the previous status or that a feature "shipped". Don't bump Proposed → Accepted just because a design merged.

## Index

| # | Title | Status |
|---|---|---|
| [0001](0001-gpu-accelerated-device-evaluation.md) | GPU-Accelerated Device Evaluation | Accepted (partial) (2026-05-28) |

(`new decision` appends rows for you. `uv run docs/scripts/doctool.py lint` checks this table against the files.)

## Adding a new decision

```sh
uv run docs/scripts/doctool.py new decision "Short title"
```

That copies `0000-template.md` to the next free number and adds the index row. The four-digit number is a stable ID, cited as "Decision NNNN" from docs, code, and commits. Never reuse or renumber it, including after a record moves to `archive/`. Then add the `**Current architecture:**` link to the area doc this decision shapes, and a link back from that area doc's "Decisions behind this" section.
