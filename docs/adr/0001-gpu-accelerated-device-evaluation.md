# ADR-0001: GPU-Accelerated Device Evaluation

## Status

Proposed

## Date

2026-05-28

## Context

VACASK is a SPICE circuit simulator written in C++. Profiling shows that for
circuits at scale, the Newton-Raphson (NR) solver dominates total simulation
time (78-99.5%), and within NR, **device evaluation is the primary bottleneck**:

| Benchmark | Unknowns | NR % of total | Eval % of NR | Solver % of NR |
|-----------|----------|---------------|--------------|----------------|
| rc        | 3        | 78%           | 31%          | 34%            |
| graetz    | 9        | 86%           | 46%          | 40%            |
| ring      | 47       | 99%           | 89%          | 9%             |
| c6288     | 25,380   | 99.5%         | 90%          | 10%            |

For c6288 (10,112 PSP103 transistors, 25k unknowns), device eval takes 41.3s
vs 4.5s for KLU factorization+solve. Device eval is embarrassingly parallel
across device instances — each transistor's eval is independent.

### Prior art: VAJAX

The VAJAX project (JAX-based SPICE simulator) validated the GPU acceleration
approach but revealed framework overhead as the dominant cost:

- **JAX host-sync**: `lax.while_loop` compiles to XLA's `kWhile` thunk, which
  forces a CPU round-trip every NR iteration (1.3s/step on Tesla T4)
- **Sprux Metal solver**: Raw solve was fast (7.6ms for c6288) but JAX callback
  overhead dominated (122ms/step total)
- **MLX backend**: 200ms/step from Python loop overhead vs JAX 0.2ms/step

The lesson: the algorithms work, but Python/JAX/XLA framework overhead defeats
the GPU benefit for anything smaller than massive circuits.

### Alternatives considered

**CubeCL (Rust GPU kernel language)**: Attractive for cross-platform GPU
(CUDA/Metal/Vulkan via WGPU). Supports control flow in kernels (while loops,
conditionals). However, no sparse linear algebra, no scientific computing
ecosystem, and Burn's sparse tensor PR was closed as stale. Too risky as a
foundation.

**Rewrite in Rust**: OpenVAF is already Rust, but rewriting the entire simulator
loses VACASK's battle-tested NR loop, convergence heuristics (pnjlim/fetlim),
and full feature set (noise, AC, parameter sweeps).

## Decision

Add GPU-accelerated device evaluation to VACASK via a C/Rust FFI layer, keeping
VACASK's existing C++ simulation loop intact.

### Architecture

```
VACASK C++ simulation loop (unchanged)
  |
  |-- evalAndLoad() ---> GPU device eval dispatcher (new)
  |                        |-- CUDA backend (cuBLAS/custom kernels)
  |                        |-- Metal backend (via WGPU or direct)
  |                        |-- CPU fallback (existing OSDI path)
  |
  |-- jac.factor()/solve() ---> KLU (unchanged initially)
  |                              |-- future: GPU sparse solver
```

### Key design choices

1. **Surgical replacement**: Only replace `evalAndLoad()` with GPU dispatch.
   NR loop, timestep control, convergence — all stay in C++.

2. **OpenVAF as GPU codegen source**: OpenVAF already compiles Verilog-A to
   native code. Extend it (or post-process its output) to generate GPU kernels
   for batched device evaluation.

3. **Batch all instances of same model**: Group all PSP103 instances (or
   resistors, diodes, etc.) and evaluate them in one GPU kernel launch.
   This is the `vmap` pattern from VAJAX, implemented natively.

4. **Data stays on device**: Upload circuit topology and model parameters once.
   Node voltages and Jacobian values stay in GPU memory across NR iterations.
   Only download waveform samples periodically.

5. **CPU/GPU threshold**: Small circuits stay on CPU (existing path). GPU
   dispatch only when instance count justifies kernel launch overhead.
   VAJAX used 500 nodes; actual threshold TBD via benchmarking.

6. **Solver stays on CPU initially**: KLU at 10% of NR time is not the
   bottleneck. GPU sparse solver is a future optimization, not a prerequisite.

### Memory architecture

**Unified memory (Apple Silicon, AMD APU)**:
- Shared address space, no explicit transfers
- GPU kernels read/write node voltages directly
- Simplest path; Sprux already validated this model

**Discrete GPU (NVIDIA CUDA)**:
- Upload once at sim start: netlist topology, model parameters, sparsity pattern
- Keep on-device across NR iterations: voltage vector, Jacobian values
- Download waveform samples in async batches
- Use CUDA graphs to record eval→assemble sequence, replay per NR step

## Consequences

### Benefits
- 10-50x speedup potential on device eval for large circuits (c6288, mul64)
- Reuses VACASK's proven simulation infrastructure
- Incremental: can ship CPU-only and add GPU backends over time
- OpenVAF GPU codegen benefits the broader OpenVAF ecosystem

### Risks
- OpenVAF GPU codegen is novel engineering (no prior art for VA→GPU)
- CUDA graphs with dynamic control flow (variable NR iterations) need care
- Maintaining C++/Rust FFI boundary adds build complexity
- GPU memory limits may constrain maximum circuit size

### What this does NOT change
- NR convergence algorithms (pnjlim, fetlim)
- Adaptive timestepping (LTE, predictor)
- Output format (rawfiles)
- Netlist parsing
- Any analysis type other than transient (initially)
