# Plan — GPU device evaluation: accelerate evalAndLoad() via batched GPU dispatch

**Status:** Active.

<!--
Status lifecycle:
  Proposed -> Active -> Closed (YYYY-MM-DD)
Update in place; don't stack past states.
-->

## Goal

Accelerate VACASK's device evaluation (89-90% of NR time for large circuits) by dispatching batched device eval to GPU, keeping VACASK's proven NR loop, timestep control, and convergence algorithms unchanged. Implements [Decision 0001](../architecture/decisions/0001-gpu-accelerated-device-evaluation.md).

## Prerequisites

- Decision 0001 accepted
- VACASK builds and passes benchmarks on gpu-acceleration branch
- Profiling data collected (done -- see findings in Decision 0001)

## Where things stand (2026-07-15)

- Phase 0 (Instrumentation): **done.** WS0.1 + WS0.2 shipped (per-step timing under `tran_debug>=2`, per-device times via `devacct`); WS0.3 baseline collected across all 5 benchmarks (results below).
- **Spike `openvaf-gpu-codegen`: RESOLVED (2026-07-15)** — codegen works (OpenVAF MIR→MSL), f32 eval is accurate, and the **surgical eval-only architecture is vindicated**. The decider is the **sync floor: 0.083 ms**, not 80 ms — the GPU-resident-NR-loop's premise (~80 ms/step host sync) was a JAX/XLA artifact. Decision 0001's original architecture stands. Resolves open questions 1 and 2 below.
- **What is NOT yet known (the measurement audit, 2026-07-15):** the headline "~10× eval" is a *compute-only, resistive-only* microbenchmark vs the CPU's *full* eval+load (resistive+reactive, scatter included) — an optimistic upper bound, not an eval-replacement number. And the one wall-clock datapoint (WS0.3) is **71% one-time elaboration**, so it tops out at 1.33× wall; the GPU-eval win only converts to wall-clock on runs long enough to amortise c6288's 103 s elaboration. Grounding these is the **[real-world measurement ladder](#real-world-measurement-ladder-2026-07-15)** below.
- Phase 1 (Abstract eval interface): **not started — now the critical path.** WS1's gather/scatter is the last unquantified cost between here and a real speedup.
- Phase 2 / Phase 3 (GPU backends): not started. **Metal now leads** (see WS2/WS3).

### Consequences of the spike for this plan

1. **WS3 (Metal) and WS2 (CUDA) swap priority.** The spike built the Metal path — a working MIR→MSL emitter (`emit_msl3.py`) and a validated psp103 `.metallib`. Metal is also the dev machine and needs no explicit transfers. CUDA becomes the follow-on, near-free from the same MIR frontend via LLVM retarget.
2. **WS4 (GPU sparse solver) stays deferred** — on compute grounds (factor+solve ≈ 3.5 s absolute on c6288 vs eval's 36 s). The one thing that could pull it forward is *assembly*: if WS1 shows scatter into KLU is costly, an on-device solver would let the Jacobian skip the scatter. Contingent on the WS1 number, not assumed.
3. **WS1 gains an exit criterion**: measure the Jacobian scatter cost (566k f32/iteration on c6288 into KLU compressed-column). That number, not the kernel, now decides the achievable speedup.

### Real-world measurement ladder (2026-07-15)

> **Noted before the docs-discipline migration so it is not lost.** The spike proved the *architecture* (surgical) via the 0.083 ms sync floor. It did **not** measure real-world performance: the "~10×" is a resistive-only, compute-only, best-of-N kernel microbenchmark compared against the CPU's *complete* eval+load, on a benchmark that is 71% one-time elaboration. Before any speedup claim is quoted as real, climb this ladder in order — each rung makes the comparison more honest, and later rungs depend on earlier ones.

1. **Representative run length.** The WS0.3 benchmark (`tran stop=2n step=2p`, ~1000 steps) is elaboration-bound (103 s elab / 40 s NR), so it caps at 1.33× wall regardless of eval. Pick a `stop=` long enough that NR ≥ ~90% of wall (elaboration amortised), so wall-clock speedup is even a meaningful quantity. Report the elaboration fraction alongside every wall number.
2. **Apples-to-apples kernel.** Emit the **reactive (ddt)** residual+Jacobian so the GPU kernel computes what the CPU `evalAndLoad` computes. Until then, compare only against a resistive-only CPU reference (not the full `tevalload`).
3. **End-to-end eval replacement, not compute.** Measure **gather + kernel + scatter-into-KLU** against the CPU's full `tevalload` (10.3 ms/iter on c6288). The scatter of 566k f32 into compressed-column is the likely bottleneck (this is WS1's deliverable). This rung produces the first defensible eval-replacement speedup.
4. **Realistic dispatch statistics.** Report mean / p95 per-dispatch (one commit+wait per NR iteration), not best-of-N. Modify `harness.mm` to emit the distribution.
5. **Whole-run wall clock.** With 1–4 in place, run the full transient and compare wall against the CPU baseline at the same run length — the only number a user actually experiences. Separate one-time (elaborate) from per-step (NR) in the report.

**Baselines to beat, from the fresh 2026-07-15 c6288 run** (`tran stop=2n step=2p`): wall 144.3 s; elaborate 103.1 s; NR 39.9 s; PSP103 eval+load 36.1 s = **10.3 ms/NR-iter** (resistive+reactive, scatter included); refactor+solve 3.5 s. Verify via `vacask --skip-embed --skip-postprocess --no-output runme.sim` in `benchmark/c6288/vacask` + `print stats`.

**Separately: elaboration (103 s, ~per-instance PSP103 init) is now the largest cost on short/medium runs and is unaddressed.** Not part of this plan's GPU-eval scope; flagged for a future decision record. Confirming it is init-dominated (vs KLU symbolic) needs a sub-timer in elaborate.

### Build note (macOS / Apple Silicon)

The `~/.local/bin/openvaf-r` was built against old LLVM and crashes with
"Unsupported stack probing method" on `.va` compilation. A clean `openvaf-r`
built from the OpenVAF `llvm21` branch (which carries the ARM64 correctness
fixes) against Homebrew `llvm@21` works. Configure VACASK with:

```sh
cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release \
  -DFLEX_INCLUDE_DIR=/opt/homebrew/opt/flex/include \
  -DOPENVAF_DIR=<openvaf>/target/release -UOPENVAF_COMPILER
```

(`-UOPENVAF_COMPILER` clears the cached `find_program` result; `FLEX_INCLUDE_DIR`
satisfies FindFLEX's singular var.)

## Workstreams

### WS0 -- Instrumentation & baseline (1-2 weeks)

**Status:** Done -- WS0.1 + WS0.2 shipped commit `6d6378c`; WS0.3 baseline collected 2026-05-28 (results below).

Add per-timestep timing breakdown and per-device-type eval profiling. Establishes the baseline that all subsequent work is measured against.

0.1. Add per-timestep timing to `TranCore::coroutine()` (gated behind `tran_debug >= 2`). Use existing `acctPrevPoint` slot to compute deltas. Output per step: eval_time, factor_time, solve_time, nr_iters, accepted/rejected.

0.2. Enable `devacct = true` in `acct.h` and measure per-device-type eval times. Need to understand: how much time is PSP103 vs resistor vs vsource?

0.3. Run all benchmarks (rc, graetz, ring, c6288, mul64 if available) with instrumentation. Record per-step distributions, not just totals.

**Deliverables:**

- Per-step timing output in simulator
- Baseline timing data for all benchmarks

**Exit criteria:**

- Can produce per-step eval/factor/solve breakdown for any benchmark
- Per-device-type eval times available

#### WS0.3 baseline results (2026-05-28)

Methodology: clean run per benchmark (`vacask --skip-embed --skip-postprocess
--no-output runme.sim`, `print stats`), Release build on Apple M-series, single
run. Per-step distribution from a second run with `options tran_debug=2`; the
1e6-step cases (rc, graetz, mul) truncated to a ~2000-step representative sample
to avoid 100 MB logs (noted below). Times in seconds. devacct compiled in.

**Aggregate (clean runs -- authoritative):**

| Benchmark | Unknowns | Wall (s) | Steps | NR iters | Eval (s) | Refactor (s) | Solve (s) | NR total (s) | Eval/NR | Dominant device, t/call |
|-----------|----------|----------|-------|----------|----------|--------------|-----------|--------------|---------|-------------------------|
| rc        | 3        | 1.05     | 1.005M| 2.01M    | 0.444    | 0.131        | 0.045     | 0.817        | 54%     | R/C, ~4.4e-8            |
| mul       | --       | 1.06     | 500k  | 1.00M    | 0.566    | 0.176        | 0.089     | 0.932        | 61%     | diode, 2.4e-7          |
| graetz    | 9        | 2.42     | 1.0M  | 2.0M     | 1.100    | 0.425        | 0.124     | 1.857        | 59%     | diode, 2.5e-7          |
| ring      | 47       | 1.43     | 26070 | 81878    | 1.001    | 0.070        | 0.024     | 1.115        | **90%** | psp103, 1.2e-5         |
| c6288     | 25380    | 145.1    | 1024  | 3506     | 37.08    | 3.093        | 0.508     | 40.95        | **90.6%** | psp103, 1.06e-2      |

**Per-step distribution (`tran_debug=2`; % of eval+factor+solve):**

| Benchmark | Sample steps | Eval % | Factor % | Solve % | Eval/step mean | Eval/step max |
|-----------|--------------|--------|----------|---------|----------------|---------------|
| rc        | 2016 (trunc) | 71.6   | 21.0     | 7.5     | 4.4e-7         | 1.1e-5        |
| graetz    | 2003 (trunc) | 67.5   | 25.4     | 7.1     | 1.2e-6         | 2.8e-6        |
| mul       | 2063 (trunc) | 65.9   | 23.6     | 10.5    | 1.2e-6         | 7.5e-6        |
| ring      | 26071 (full) | 91.3   | 6.5      | 2.2     | 3.8e-5         | 1.5e-4        |
| c6288     | 1033 (full)  | 90.4   | 8.2      | 1.4     | 5.2e-2         | 1.3e-1        |

(Per-step % uses eval+factor+solve as denominator and so runs higher than the
aggregate Eval/NR, which divides by full NR time including loop overhead.)

**Findings:**

1. **Eval dominance tracks device-model complexity, not unknown count.** ring
   has only 47 unknowns yet is already 90% eval, because a PSP103 MOSFET eval
   costs ~1.2e-5 s/call vs ~4e-8 for R/C and ~2.5e-7 for a diode -- 300-1000x
   more. c6288 confirms the same 90% at 25k unknowns / 10k transistors.
2. **For the GPU-target circuits (PSP-MOSFET-heavy: ring, c6288), device eval
   is ~90% of NR and KLU factor+solve is only 8-10%.** Offloading eval addresses
   the dominant cost; this confirms WS4 (GPU sparse solver) is correctly
   deferred -- even a 10x eval speedup leaves the solver at a few seconds
   absolute for c6288.
3. **Simple-device circuits (rc/graetz/mul) sit at 54-65% eval** and are *not*
   the target -- their factor/solve fraction is too large for eval offload to
   pay off, and they are tiny in absolute terms (~1-2 s total).
4. **Per-step eval cost on c6288 is 17.6-132 ms**, scaling with NR iterations
   per step (2-6). This is the granularity a GPU dispatch must beat after
   launch/transfer overhead -- comfortably above typical kernel-launch latency.

These confirm Decision 0001's decision-time table (the load-bearing ring/c6288 rows
match within rounding). The small-circuit rows differ from the decision's original
profiling (different options/run), which is expected and immaterial -- those
circuits are not GPU targets.

### WS1 -- Abstract the device eval interface (2-3 weeks)

**Status:** Not started.

Refactor `circuit.evalAndLoad()` from sequential iteration into gather->eval->scatter pattern, then implement a batched eval interface with CPU reference implementation.

1.1. Extract device eval from evalAndLoad into gather (collect node voltages into contiguous arrays per model type), eval (call OSDI on batch), scatter (write currents/charges/Jacobian back to MNA).

1.2. Define `BatchedDeviceEval` interface:
```cpp
struct BatchedDeviceEval {
    // Input: node voltages (n_instances x n_terminals)
    // Input: model params (shared), instance params (n_instances x n_instance_params)
    // Output: currents, charges, Jacobian entries (per instance)
    virtual void eval(/* ... */) = 0;
};
```

1.3. CPU reference implementation using existing OSDI calls. Must produce bit-identical results to current sequential path.

**Deliverables:**

- Refactored evalAndLoad with gather/scatter
- BatchedDeviceEval interface
- CPU reference implementation

**Exit criteria:**

- All benchmarks produce identical results via batched CPU path
- No performance regression on CPU path
- **Gather/scatter cost measured on c6288** (added 2026-07-15). The GPU kernel is no longer the unknown — this is. Gather should be cheap by construction: only ~36 of each instance's 2090 kernel inputs change per NR iteration (the other 2054 are bias-independent cached init values, uploaded once), so it is a ~1.5 MB unified-memory write. **Scatter is the risk**: 566k f32 Jacobian entries into KLU compressed-column, every iteration. Report it against the 10.3 ms/iteration CPU eval+load it must fit inside — if scatter approaches a few ms, the compute-only ~10× win collapses and WS3 should not start until the assembly path is redesigned (see Decision 0001 Risks). This is rung 3 of the measurement ladder.

### WS2 -- GPU backend: CUDA (3-4 weeks)

**Status:** Not started. **Deprioritised behind WS3 (2026-07-15)** — the spike built and validated the Metal path, and CUDA is near-free from the same MIR frontend once WS1+WS3 prove the architecture end-to-end. Option A (OpenVAF CUDA codegen) and Option C (CubeCL) below are superseded by the spike's resolution: the answer is a MIR-walking emitter with an NVPTX backend.

Wire GPU device eval into `OpNRSolver::buildSystem()` with CUDA backend, starting with PSP103.

Three options for kernel generation (decision deferred to spike):
- **Option A**: OpenVAF CUDA codegen (long-term right answer, significant compiler work)
- **Option B**: Translate OSDI to CUDA (faster prototype)
- **Option C**: CubeCL kernels (rejected per Decision 0001)

Memory: allocate at elaboration, pin host, double-buffer voltages, use CUDA graphs.

Jacobian assembly: scatter on GPU (atomic adds) or download and scatter on CPU.

**Deliverables:**

- c6288 running with GPU device eval on CUDA
- CPU/GPU path selection based on circuit size threshold

**Exit criteria:**

- c6288 GPU results match CPU within tolerance
- Measurable speedup on device eval for c6288

### WS3 -- GPU backend: Metal/Apple Silicon (2-3 weeks)

**Status:** Not started. **Now the lead GPU backend** (2026-07-15, was behind WS2).

Direct Metal compute shaders (not WGPU, not CubeCL) fed by the MIR→MSL emitter the spike built. Leverage unified memory: no explicit transfers, and cached init values uploaded once rather than per iteration.

Carry-in from the spike (all in `~/Code/ChipFlow/vajax/spikes/msl-codegen/`, untracked scratch — graduating it into VACASK proper is part of this WS, and warrants a `docs/gpu-codegen.md` design doc):

- `emit_msl3.py` — MIR→MSL emitter (SCCP + backward-liveness DCE + dominator phi resolution); psp103 eval → 18.9k lines MSL, compiles clean.
- `psp103_v3.metal` / validated `.metallib`; `harness.mm` Metal runner; `bench_psp103_scale.py` (the compute-only ~10× microbenchmark) + `trivial.metal` (the 0.083 ms sync floor).
- **Known gaps to close here:** the reactive (ddt) path is stubbed (v3 emits resistive only — required for a transient-correct kernel); f64 host-side init must feed the f32 kernel's cache; deterministic (bucketed) scatter-add, since f32 addition is non-associative under duplicate-index parallel scatter.

**Deliverables:**

- Same benchmarks running on Apple Silicon GPU
- Codegen graduated from spike scratch into the VACASK build

**Exit criteria:**

- Benchmarks pass on Apple Silicon with GPU path
- c6288 GPU results match CPU within tolerance
- **Speedup measured end-to-end** per the [measurement ladder](#real-world-measurement-ladder-2026-07-15) — not the compute microbenchmark. The spike's ~10× eval-compute is an optimistic upper bound (resistive-only, scatter excluded); the real target is gather+full-eval+scatter vs CPU `tevalload` (10.3 ms/iter), on a run long enough to amortise elaboration.

### WS4 -- GPU sparse solver (future, optional)

**Status:** Deferred.

Only pursue if profiling shows solver becoming bottleneck after eval is on GPU. For c6288, solver is currently 10% of NR -- even with 10x eval speedup, solver becomes ~50% but is still only 4.5s absolute.

### WS5 -- Advanced optimizations (future)

**Status:** Deferred.

Multi-simulation batching (Monte Carlo/corners on GPU), CUDA graphs for full NR loop, device eval bypass on GPU (predicate kernel to skip unchanged instances).

## Pre-merge cleanup (spike scaffolding on `gpu-acceleration`)

Carried over from the resolved spike's handoff. None of this blocks WS1; all of it blocks **merging `gpu-acceleration` to `main`**:

- **`devacct = true` is compiled into `include/acct.h:11`** (commit `6d6378c`) — always-on per-device timing overhead. Revert to `false` or runtime-gate it. This is also a suspect in the wall-clock-vs-NR gap flagged in Decision 0001's Risks (c6288 wall 145.1 s vs NR 40.95 s); triage that before assuming the overhead is benign.
- **The f64 eval-dump is removable spike scaffolding**: `EvalDumpSink` (`include/evaldump.h`, `lib/evaldump.cpp`), `OsdiInstance::dumpEvalIO` (`lib/osdiinstance.cpp`), and the TranCore post-acceptance pass (`lib/coretran.cpp`). Env-gated (`VACASK_EVAL_DUMP=<path>`, `VACASK_EVAL_DUMP_STRIDE=N`), hot path costs one branch. It produced the f64 ground truth that validated the f32 kernel, and it will be needed again to validate c6288 in WS3 — so keep it until WS3 closes, then remove it (cf. `devacct`) or promote it to a supported debug facility.

## Open questions

1. ~~**OpenVAF codegen strategy**~~ — **RESOLVED (spike, 2026-07-15):** a MIR-walking emitter with pluggable backends. Not an OpenVAF compiler extension, not OSDI post-processing. See Decision 0001 design choice #2.
2. ~~**Which GPU first?**~~ — **RESOLVED (spike, 2026-07-15): Metal.** The emitter and a validated psp103 `.metallib` already exist; unified memory removes explicit transfers; it is the dev platform. CUDA follows from the same frontend.
3. **Build system**: CMake + Rust FFI via corrosion, or separate Rust static lib? *(Still open — but note the spike's emitter is Python, so the shipping form of the codegen is itself undecided: build-time codegen vs a runtime FFI layer.)*
4. **Minimum viable circuit size for GPU** — *partially answered*: the dispatch floor is **0.083 ms** and psp103 eval at N=126 is 0.357 ms, so break-even against CPU eval (~1.05 µs/instance) is on the order of a few hundred PSP103 instances. vajax's 500-node threshold is the right ballpark; confirm against the real gather/scatter in WS1.
5. **Jacobian format**: VACASK uses KLU compressed column. Most efficient GPU assembly path? *(Now the highest-value open question — see WS1.)*

## Phase exit criteria

When all of these are true, this plan closes:

- [x] WS0 baseline data collected (2026-05-28)
- [x] Open questions 1-2 resolved (spike `openvaf-gpu-codegen`, 2026-07-15: MIR→MSL emitter; Metal first)
- [ ] WS1 batched eval interface working with CPU reference, **with gather/scatter cost measured on c6288**
- [ ] WS2 or WS3 GPU backend demonstrating speedup on c6288
- [ ] All benchmarks pass with GPU path enabled
- [ ] Open question 3 resolved (build system / shipping form of the codegen)
- [x] Wall-clock-vs-NR gap triaged (2026-07-15): it is one-time elaboration (103 s), amortises over run length; see the measurement ladder
- [ ] Real-world measurement ladder climbed (rungs 1–5) — an end-to-end, apples-to-apples, representative-length speedup exists

## References

- [Decision 0001](../architecture/decisions/0001-gpu-accelerated-device-evaluation.md) -- architectural decision
- Key source files: `lib/coretran.cpp` (timestep loop), `lib/nrsolver.cpp` (NR loop), `lib/coreopnr.cpp` (buildSystem/evalAndLoad), `lib/circuit.cpp` (evalAndLoad iterates devices), `include/acct.h` (timing)
