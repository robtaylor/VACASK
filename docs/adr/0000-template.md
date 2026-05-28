# ADR NNNN — <Short title that names the decision>

**Status:** Proposed.

<!--
Status lifecycle:
  Proposed → Accepted (YYYY-MM-DD) → [Superseded by ADR MMMM (YYYY-MM-DD) | Withdrawn (YYYY-MM-DD)]

When scope expands, amend the Status line in place:
  **Status:** Accepted (YYYY-MM-DD). Scope expanded YYYY-MM-DD — see Decision §N.

Add an inline Amendment section near the top — see https://robtaylor.github.io/claude-project-discipline/adr.html
-->

**TL;DR.** In the context of <use case>, facing <concern>, we chose <option> to achieve <quality>, accepting <downside>.

## Context

What was the situation? What constraints forced a decision? Be concrete — name the systems, the data formats, the friction. A reader six months from now should be able to reconstruct the world this decision was made in.

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

What this buys us, what it costs us. Be honest about the costs.

- <Positive consequence>
- <Negative consequence or cost>
- <Operational impact (CI, build, runtime, …)>

## Walk-back options

Conditions under which this decision should be revisited, and how to revisit it cheaply.

- **If <X happens>** — <how to walk back>.
- **If <Y measurement comes back unexpected>** — <how to walk back>.

## Links

- `<related-design-doc.md>` — <why it's relevant>
- ADR <MMMM> — <related decision>
- `<spike.md>` — <if this decision was informed by a spike>
