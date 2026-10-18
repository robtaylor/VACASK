# Architecture

The current design of VACASK, by area. Each area is a **reference** doc: what the code does today, present tense, no rationale and no dead history. The *why* behind each choice lives in the [decision records](decisions/README.md), linked from the section it explains. The *not-yet-built* work lives in each doc's Implementation status section and in [plans](../plans/).

If you're an agent or a new reader, start here: this page maps every area to its current-state doc and the decisions behind it.

| Area | Reference | Decisions behind it |
|---|---|---|
| GPU device evaluation | [gpu-device-eval](gpu-device-eval.md) | [0001](decisions/0001-gpu-accelerated-device-evaluation.md) |

VACASK's wider architecture (parser, elaboration, the NR core, analyses, OSDI device layer, output) is not yet captured here. Add an area doc per part as work touches it; this GPU row is the first.
