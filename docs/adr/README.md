# Architecture Decision Records

Numbered records of non-obvious design choices, kept forever, amended in place when scope expands, never deleted.

For the full convention, see [Claude Project Discipline -- ADRs](https://robtaylor.github.io/claude-project-discipline/adr.html).

## Format

```markdown
# ADR NNNN -- Short title

**Status:** Accepted (YYYY-MM-DD).

## Context
What was the situation? What constraints forced a decision?

## Decision
What we chose, with enough specificity to act on.

## Consequences
What this buys us, what it costs us, what would trigger a walk-back.

## Walk-back options
(Optional, recommended.) Conditions to revisit.

## Links
Cross-references.
```

When scope expands, amend in place rather than rewriting:

```markdown
**Status:** Accepted (2026-05-28). Scope expanded YYYY-MM-DD -- see Decision §N.
```

When superseded, update Status on the old ADR and add a new ADR that references it. Never delete the old one.

## Index

| ADR  | Title | Status |
|------|-------|--------|
| 0001 | [GPU-accelerated device evaluation](0001-gpu-accelerated-device-evaluation.md) | Accepted 2026-05-28 |

## Template

Copy `0000-template.md` to `NNNN-<short-name>.md` for the next available number.
