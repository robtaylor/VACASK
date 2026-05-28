# Plan — GPU device evaluation: accelerate evalAndLoad() via batched GPU dispatch

**Status:** Active.

<!--
Status lifecycle:
  Proposed -> Active -> Closed (YYYY-MM-DD)
Update in place; don't stack past states.
-->

## Goal

Accelerate VACASK's device evaluation (89-90% of NR time for large circuits) by dispatching batched device eval to GPU, keeping VACASK's proven NR loop, timestep control, and convergence algorithms unchanged. Implements [ADR 0001](../adr/0001-gpu-accelerated-device-evaluation.md).

## Prerequisites

- ADR 0001 accepted
- VACASK builds and passes benchmarks on gpu-acceleration branch
- Profiling data collected (done -- see findings in ADR 0001)

## Where things stand (2026-05-28)

- Phase 0 (Instrumentation): **done.** WS0.1 + WS0.2 shipped (per-step timing under `tran_debug>=2`, per-device times via `devacct`); WS0.3 baseline collected across all 5 benchmarks (results below).
- Phase 1 (Abstract eval interface): not started
- Phase 2 (CUDA backend): not started
- Phase 3 (Metal backend): not started

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

These confirm ADR 0001's decision-time table (the load-bearing ring/c6288 rows
match within rounding). The small-circuit rows differ from the ADR's original
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

### WS2 -- GPU backend: CUDA (3-4 weeks)

**Status:** Not started.

Wire GPU device eval into `OpNRSolver::buildSystem()` with CUDA backend, starting with PSP103.

Three options for kernel generation (decision deferred to spike):
- **Option A**: OpenVAF CUDA codegen (long-term right answer, significant compiler work)
- **Option B**: Translate OSDI to CUDA (faster prototype)
- **Option C**: CubeCL kernels (rejected per ADR 0001)

Memory: allocate at elaboration, pin host, double-buffer voltages, use CUDA graphs.

Jacobian assembly: scatter on GPU (atomic adds) or download and scatter on CPU.

**Deliverables:**

- c6288 running with GPU device eval on CUDA
- CPU/GPU path selection based on circuit size threshold

**Exit criteria:**

- c6288 GPU results match CPU within tolerance
- Measurable speedup on device eval for c6288

### WS3 -- GPU backend: Metal/Apple Silicon (2-3 weeks)

**Status:** Not started.

Leverage unified memory (no explicit transfers). Options: direct Metal compute shaders, WGPU, or CubeCL.

**Deliverables:**

- Same benchmarks running on Apple Silicon GPU

**Exit criteria:**

- Benchmarks pass on Apple Silicon with GPU path
- Measurable speedup

### WS4 -- GPU sparse solver (future, optional)

**Status:** Deferred.

Only pursue if profiling shows solver becoming bottleneck after eval is on GPU. For c6288, solver is currently 10% of NR -- even with 10x eval speedup, solver becomes ~50% but is still only 4.5s absolute.

### WS5 -- Advanced optimizations (future)

**Status:** Deferred.

Multi-simulation batching (Monte Carlo/corners on GPU), CUDA graphs for full NR loop, device eval bypass on GPU (predicate kernel to skip unchanged instances).

## Open questions

1. **OpenVAF codegen strategy**: Extend OpenVAF compiler to emit GPU kernels, or separate VA-to-GPU tool? Needs a spike.
2. **Which GPU first?** CUDA (better tooling, HPC standard) or Metal (unified memory, dev platform)?
3. **Build system**: CMake + Rust FFI via corrosion, or separate Rust static lib?
4. **Minimum viable circuit size for GPU**: VAJAX used 500 nodes. Need kernel launch overhead benchmark.
5. **Jacobian format**: VACASK uses KLU compressed column. Most efficient GPU assembly path?

## Phase exit criteria

When all of these are true, this plan closes:

- [x] WS0 baseline data collected (2026-05-28)
- [ ] WS1 batched eval interface working with CPU reference
- [ ] WS2 or WS3 GPU backend demonstrating speedup on c6288
- [ ] All benchmarks pass with GPU path enabled
- [ ] Open questions 1-3 resolved (via spikes or ADR amendments)

## References

- [`../adr/0001-gpu-accelerated-device-evaluation.md`](../adr/0001-gpu-accelerated-device-evaluation.md) -- architectural decision
- Key source files: `lib/coretran.cpp` (timestep loop), `lib/nrsolver.cpp` (NR loop), `lib/coreopnr.cpp` (buildSystem/evalAndLoad), `lib/circuit.cpp` (evalAndLoad iterates devices), `include/acct.h` (timing)
