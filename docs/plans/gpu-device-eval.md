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

- Phase 0 (Instrumentation): not started -- aggregate profiling done via acct.h, per-step timing not yet added
- Phase 1 (Abstract eval interface): not started
- Phase 2 (CUDA backend): not started
- Phase 3 (Metal backend): not started

## Workstreams

### WS0 -- Instrumentation & baseline (1-2 weeks)

**Status:** Not started.

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

- [ ] WS0 baseline data collected
- [ ] WS1 batched eval interface working with CPU reference
- [ ] WS2 or WS3 GPU backend demonstrating speedup on c6288
- [ ] All benchmarks pass with GPU path enabled
- [ ] Open questions 1-3 resolved (via spikes or ADR amendments)

## References

- [`../adr/0001-gpu-accelerated-device-evaluation.md`](../adr/0001-gpu-accelerated-device-evaluation.md) -- architectural decision
- Key source files: `lib/coretran.cpp` (timestep loop), `lib/nrsolver.cpp` (NR loop), `lib/coreopnr.cpp` (buildSystem/evalAndLoad), `lib/circuit.cpp` (evalAndLoad iterates devices), `include/acct.h` (timing)
