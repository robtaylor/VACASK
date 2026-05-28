# Plan — <Phase or topic>: <one-line goal>

**Status:** Proposed.

<!--
Status lifecycle:
  Proposed → Active → Closed (YYYY-MM-DD)
Update in place; don't stack past states.
-->

## Goal

What this phase / topic delivers, in 1–3 sentences. Reference the Architecture Decision Record(s) (ADRs) this plan implements.

## Prerequisites

- ADR NNNN <accepted | proposed>
- Fixture / data / dependency in place
- Previous phase's exit criteria met

## Where things stand (YYYY-MM-DD)

The single most-recent status snapshot. Updated in place. Past states are in `git log`, not stacked here.

- Workstream A: <one-line current state>
- Workstream B: <one-line current state>

## Workstreams

### WS1 — <Name>

**Status:** <In flight | Shipped commit `abc1234` | Blocked on X>

<Concrete scope. Reference ADRs and design docs; don't duplicate them.>

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

- `<adr/NNNN-…>.md` — relevant decision
- `<design-doc>.md` — relevant design
- `<predecessor-plan>.md` — what this plan picks up from
