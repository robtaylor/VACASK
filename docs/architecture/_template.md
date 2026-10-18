# <Area name>

*Reference: the current design as it stands in the code. The why lives in the decision records linked from each section; the not-yet-built work is under [Implementation status](#implementation-status). Present tense means "true today". If a claim here is stale, the code moved and this doc didn't.*

<!--
Rules for this doc:
- Every present-tense sentence is a claim about the code. Check it before writing it.
- No rationale ("we chose X because Y"): link to the decision instead.
- No dead history ("previously this used Z"): that's git log and decisions/archive/.
- Keep each section self-contained. Restate the subject rather than "as above".
-->

One paragraph: what this area does and where it sits in the system. Rationale for the overall shape: [Decision NNNN](decisions/NNNN-<name>.md).

## <Component or concept>

What it does today. File references where they help (`src/<path>`). A decision link wherever a reader will ask "why?".

## <Component or concept>

(...)

## Constraints

Hard limits a contributor needs to know: sizes, supported platforms, unsupported inputs.

## Implementation status

**Implemented:**

| Part | Where |
|---|---|
| <part> | `src/<path>` |

**Deferred:**

| Part | Gap | Tracking |
|---|---|---|
| <part> | <what's missing, specifically> | `docs/plans/<plan>.md`, issue #N |

## Decisions behind this

- [Decision NNNN — <title>](decisions/NNNN-<name>.md)
