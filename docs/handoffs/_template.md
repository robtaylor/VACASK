# Handoff — <Topic> (one-line summary of what this session left open)

**Created:** YYYY-MM-DD
**Working tree:** clean | <state if not clean>
**Branch:** main | <branch>

<!--
Reminder: a handoff is ephemeral. At resolution, every load-bearing piece
below migrates into a docs/adr/, docs/plans/, docs/spikes/, or design-doc
home, and this file is then `git rm`'d in the same commit as the migration.

See docs/handoff-discipline.md for the migration table.
-->

## Goal & next-up

**Goal of this session:** <what you were trying to do, in 1-3 sentences>

**Next session should pick up:** <the very next concrete action, by name. Reference the plan doc section if applicable.>

**Verification command:**

```sh
<commands the next session runs to confirm this handoff's claimed state>
# Expect: <what success looks like>
```

## Done this session

| Commit | Subject | Notes |
|---|---|---|
| `<sha>` | <subject> | <one-line note> |

## Open follow-ups (priority-ordered)

### 1. <Item name> (<rough size>)

<Concrete scope. Enough detail to start cold. Link to existing plan/ADR/design-doc sections rather than reproducing them.>

### 2. <Item name> (<rough size>)

(...)

## Critical context

<Things the next session needs to know that aren't yet in the code/docs.
Be honest about what's truly load-bearing — anything obvious from
`git log` or a quick `grep` doesn't belong here.>

## References

- [`<predecessor-handoff if any>`](<path>) -- predecessor (rare; usually a smell)
- [`<plan doc>`](<path>) -- current workstream state
- [`<ADR>`](<path>) -- relevant decision

## Migration note

<When this resolves, what migrates where. Helps the next session do the
fold-and-delete cleanly. Example:

- Open follow-up 1 -> docs/plans/<plan>.md as new WS-X.Y
- Critical context bullet 2 -> code comment near <file:line>
- Verification command -> docs/plans/<plan>.md verification block
- Then `git rm docs/handoffs/<this-file>.md` in the migration commit>
