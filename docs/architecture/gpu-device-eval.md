# GPU device evaluation

*Reference: the current design as it stands in the code. The why lives in the decision records linked from each section; the not-yet-built work is under [Implementation status](#implementation-status). Present tense means "true today". If a claim here is stale, the code moved and this doc didn't.*

VACASK evaluates every device model inside the Newton-Raphson loop: each instance turns its terminal voltages into currents, charges, and Jacobian entries, which load into the sparse MNA system the solver factorises. For large circuits built from PSP103 MOSFETs this evaluation is the dominant cost. This area covers the work to offload batched device evaluation to the GPU while the NR loop, timestep control, and convergence heuristics stay in C++. The rationale for that shape is [Decision 0001](decisions/0001-gpu-accelerated-device-evaluation.md).

The GPU backend is not built yet. What exists in the tree today is the instrumentation and the f64 reference dump that measure the target and ground the accuracy work. The code generation and Metal kernel live in a separate prototype tree, not in VACASK. [Implementation status](#implementation-status) draws the line.

## Device evaluation and loading

`Circuit::evalAndLoad` (`lib/circuit.cpp:1373`) walks the device list, and for each device type whose instance count is non-zero calls that device's `evalAndLoad`. For an OSDI device this evaluates the compiled model and loads the resistive residual and Jacobian straight into the matrix. Evaluation and loading are one fused operation per device type, and the accounting times them together as `tevalload` ("Eval and load time", `lib/acct.cpp:25`).

On c6288 (10,112 PSP103 instances, 25,380 unknowns) this fused eval+load costs 10.3 ms per NR iteration and is about 90% of NR time. The reactive (charge) terms and the resistive terms are both computed in this path.

## Instrumentation

Two timing facilities measure where transient time goes.

- **Per-step breakdown** prints `nr_iter`, `eval`, `factor`, and `solve` for each accepted step. It fires when `tran_debug >= 2` through `dumpStepTiming` in `TranCore::coroutine` (`lib/coretran.cpp:1211`).
- **Per-device-type eval+load** accumulates a time and call count for each device type. It is compiled in behind the `devacct` flag (`include/acct.h`), summed in `Circuit::evalAndLoad`, and printed by `Accounting::dumpDevTimes` (`lib/acct.cpp:66`). `devacct` is currently `true` in the source, so this accounting is always on. The plan tracks reverting it to `false` before merge.

Aggregate totals come from `print stats`, which reports parse, elaboration, eval+load, NR, factor, and solve times separately (`Accounting::dumpTotal`, `lib/acct.cpp`).

## f64 device-eval reference dump

A flag-gated sink writes, per PSP103 instance at accepted transient timepoints, the OSDI input voltages, the absolute internal-node voltages, and the f64 resistive residual and Jacobian. The transient core creates the sink when `VACASK_EVAL_DUMP` names an output path, and `VACASK_EVAL_DUMP_STRIDE=N` keeps every Nth timepoint. The pieces are `EvalDumpSink` (`include/evaldump.h`), `OsdiInstance::dumpEvalIO` (`lib/osdiinstance.cpp:1381`), and a post-acceptance pass in `TranCore` (`lib/coretran.cpp:1108`).

This dump is the ground truth the generated f32 kernel is validated against. It is spike scaffolding rather than a shipping feature, and it is removed or promoted when the Metal backend lands.

## Measured baseline

These numbers are from a fresh instrumented c6288 run (`tran stop=2n step=2p`, Release build on an Apple M4 Pro, 2026-07-15). They describe VACASK's own behaviour today.

| Quantity | Value |
|---|---|
| Wall | 144.3 s |
| Elaboration (one-time) | 103.1 s |
| NR solve | 39.9 s |
| PSP103 eval+load | 36.1 s total, 10.3 ms per NR iteration |
| Refactor + solve | 3.5 s |

Elaboration is 71% of this benchmark's wall because the run is only 2 ns. It is a one-time cost, so a longer transient amortises it and NR returns to the bulk of wall. The GPU-eval work pays off in wall-clock only on runs long enough to amortise elaboration. GPU-side prototype measurements (the Metal sync floor and the resistive-eval kernel throughput) live in the [spike](decisions/spikes/openvaf-gpu-codegen.md), since they were taken outside the VACASK tree.

## Constraints

- GPU dispatch targets large PSP103-heavy circuits. Small circuits stay on the existing CPU path, below a threshold set by kernel-launch overhead.
- Apple GPUs provide no f64. The device-eval kernel runs f32, and model initialisation stays in f64. The accuracy evidence for this split is in the spike.
- Elaboration of c6288 costs 103 s and is not accelerated. On short and medium runs it is the single largest cost.

## Implementation status

**Implemented (in VACASK):**

| Part | Where |
|---|---|
| Per-step transient timing (`tran_debug>=2`) | `lib/coretran.cpp` |
| Per-device eval+load accounting (`devacct`) | `include/acct.h`, `lib/circuit.cpp`, `lib/acct.cpp` |
| f64 device-eval reference dump | `include/evaldump.h`, `lib/evaldump.cpp`, `lib/osdiinstance.cpp`, `lib/coretran.cpp` |

**Deferred (not in the tree):**

| Part | Gap | Tracking |
|---|---|---|
| Batched eval interface (gather/eval/scatter) | not started; CPU reference first | `docs/plans/gpu-device-eval.md` WS1 |
| Metal backend | MIR→MSL emitter and a validated psp103 `.metallib` exist as a prototype in a separate tree, not wired into VACASK; the reactive path is stubbed | plan WS3 |
| CUDA backend | not started | plan WS2 |
| GPU sparse solver | deferred; the CPU KLU solve is not the bottleneck | plan WS4 |

## Decisions behind this

- [Decision 0001 — GPU-Accelerated Device Evaluation](decisions/0001-gpu-accelerated-device-evaluation.md)
