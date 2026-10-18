# Plan — <Phase or topic>: <one-line goal>

**Status:** Proposed.

<!--
Status lifecycle:
  Proposed → Active → Implemented (YYYY-MM-DD)
  (also: Deferred; Withdrawn (YYYY-MM-DD) if dropped)
The first word must be one of these; doctool.py lint checks it.
Update in place; don't stack past states.
When Implemented or Withdrawn, move this file to docs/plans/archive/ with
  uv run docs/scripts/doctool.py archive docs/plans/<this-file>.md
so links to and from it are repointed.
-->

## Goal

What this phase / topic delivers, in 1–3 sentences. Reference the decision(s) this plan implements.

## Prerequisites

- Decision NNNN <accepted | proposed>
- Fixture / data / dependency in place
- Previous phase's exit criteria met

## Where things stand (YYYY-MM-DD)

The single most-recent status snapshot. Updated in place. Past states are in `git log`, not stacked here.

- Workstream A: <one-line current state>
- Workstream B: <one-line current state>

## Workstreams

### WS1 — <Name>

**Status:** <In flight | Shipped commit `abc1234` | Blocked on X>

<Concrete scope. Reference decisions and area docs; don't duplicate them.>

**Deliverables:**

- <Specific artefact: a binary, a schema, a Continuous Integration (CI) gate, a test fixture>
- <Another artefact>

**Exit criteria:**

- <Verifiable condition>
- <Another verifiable condition>

### WS2 — <Name>

**Status:** ...

(...)

## Phase exit criteria

When all of these are true, this plan closes:

- [ ] <Condition>
- [ ] <Condition>
- [ ] All workstreams shipped or explicitly deferred to a successor plan.

## References

- `architecture/decisions/NNNN-….md` — relevant decision
- `architecture/<area>.md` — current design of the area this changes
- `<predecessor-plan>.md` — what this plan picks up from
