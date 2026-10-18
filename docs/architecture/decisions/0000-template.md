# Decision NNNN — <Short title that names the decision>

**Status:** Proposed.

**Current architecture:** [<Area>](../<area>.md).

<!--
Status words and their meaning: see README.md § Status legend.

This record is immutable once accepted. Only the Status line and the
Current architecture link change. When the code moves, edit the area doc.
When the decision changes, write a new record and archive this one:
  uv run docs/scripts/doctool.py archive docs/architecture/decisions/NNNN-<name>.md

Write the body in past tense, as the state at decision time. While this is
Proposed, open the Context with one line on what the code does today.
-->

**TL;DR.** In the context of <use case>, facing <concern>, we chose <option> to achieve <quality>, accepting <downside>.

## Context

What was the situation? What constraints forced a decision? Be concrete: name the systems, the data formats, the friction. A reader six months from now should be able to reconstruct the world this decision was made in.

## Decision

What we chose, with enough specificity to act on. Numbered sub-points if multiple commitments.

1. <Commitment 1>
2. <Commitment 2>
3. <Commitment 3>

## Alternatives considered

The options we looked at and rejected, with the reason for each. Future readers should not have to re-litigate options that were already weighed and dismissed.

- **<Option A>** — rejected because <specific reason tied to this project's constraints, not generic preference>.
- **<Option B>** — rejected because <…>.
- **<Do nothing / defer>** — rejected because <…>.

## Consequences

What this bought us, what it cost us. Be honest about the costs.

- <Positive consequence>
- <Negative consequence or cost>
- <Operational impact (CI, build, runtime, …)>

## Walk-back options

Conditions under which this decision should be revisited, and how to revisit it cheaply.

- **If <X happens>** — <how to walk back>.
- **If <Y measurement comes back unexpected>** — <how to walk back>.

## Links

- `spikes/<spike>.md` — <if this decision rested on a spike>
- Decision <MMMM> — <related or superseded decision>
