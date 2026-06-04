# ADR-0001: GPU-Accelerated Device Evaluation

## Status

Accepted (2026-05-28). **Partially under revision (2026-06-01)** — see
[Amendment 2026-06-01](#amendment-2026-06-01-revisions-under-spike-validation).

## Amendment 2026-06-01 (revisions under spike validation)

Work in [spike `openvaf-gpu-codegen`](../spikes/openvaf-gpu-codegen.md) has put
four of this ADR's choices under active revision. They are **not yet re-decided**
(the spike is Open); recorded here so the Accepted decision is not silently
contradicted. At spike resolution these fold back into the sections below.

1. **Design choice #1 (surgical eval-only replacement) → likely GPU-resident NR
   loop.** vajax evidence shows the per-iteration host↔device sync, not compute,
   dominates (`while_loop` ~80 ms/step vs fused ~0.2 ms/step on M-series).
   Offloading *only* `evalAndLoad()` re-incurs that round-trip every NR iteration.
   A resident loop (eval+assembly+solve+convergence on-device, sync at accepted
   timepoints) avoids it.
2. **"Does NOT change: NR convergence (pnjlim/fetlim), timestepping" → would move
   to GPU** under a resident loop (the cost: porting those heuristics; vajax did
   this in f32 with convergence-management hacks).
3. **Design choice #6 (solver stays on CPU; GPU solver not a prerequisite) →
   becomes a prerequisite** for a resident loop — an on-device sparse solver
   (Sprux on Metal, cuDSS on CUDA) is needed to avoid the per-iter Jacobian
   round-trip, not because solve is the compute bottleneck.
4. **Design choice #2 (OpenVAF codegen: "extend or post-process") → specifically
   MIR→MSL on Metal.** No open LLVM→AIR backend; the LLVM→SPIR-V→MSL bridge
   hard-fails at f64; we own f32 lowering regardless, so a MIR-walking emitter
   with pluggable backends (MSL now; NVPTX/CUDA near-free via LLVM retarget) is
   simpler. **Metal runs pure f32** (verified in vajax) — f64 emulation is not a
   prerequisite. Codegen feasibility confirmed: psp103 eval generates + compiles
   to a `.metallib`. **f32 *accuracy* now measured — see the Q2 addendum below.**

### Q2 f32-numerics measurement (2026-06-03, literal-kernel-confirmed 2026-06-04) — the accuracy input to choices #2/#4

The spike measured f32 device-eval accuracy for PSP103 against VACASK's own f64
OSDI eval (gilbert, 126 operating points; validation gate passed — JAX-f64 ≈
VACASK-f64 to ~5e-6 residual / ~9e-5 Jacobian rel err). The literal
`psp103_v3.metal` kernel was then run on the M4 Pro and compared to the same f64
dump (2026-06-04). Result:

- **Currents (residual) are f32-safe** — rel err max ~1e-5 (kernel 1.3e-5, proxy
  1.1e-5), at the f64 cross-impl floor.
- **The resistive Jacobian is f32-safe** — the *literal kernel* gives rel err max
  **7.8e-5** (mean 4.0e-6) vs VACASK-f64. The JAX proxy's earlier "f32-lossy ~1%
  mean / 14% max" was a **JAX-CPU-f32 artifact** (full eval graph in f32:
  overflow-in-cast on intermediates, no DCE → systematic gm/gds conductance bias).
  On the exact entries the proxy flagged at ~14%, kernel-f32 ≈ VACASK-f64 to ~1e-6.
- **init must NOT run in f32** — full-f32 (init included) is catastrophic
  (residual rel err ~5e5); PSP103 init has out-of-f32-range intermediates (the
  15 `1e±100` guard constants). The kernel feeds an f64-computed cache.

**Refines amendment item #4's "Metal runs pure f32".** A generated f32 device-eval
kernel is accurate for both currents AND the resistive Jacobian (~1e-5 / ~8e-5 vs
f64) — so a resident-NR f32 loop does **not** need f64/compensated Jacobian
assembly for the resistive part; only **init** needs f64/compensated. (Convergence
≠ accuracy still holds as a general caution; here the measured accuracy is good.)
Scope of the measurement: gilbert, resistive Jacobian only — the reactive (ddt)
Jacobian (still stubbed) and c6288 are not yet covered. Full fold into the
Decision / Consequences sections waits for spike resolution (Q3 resident-loop
speed + reactive + c6288 still open). Source: spike Findings 2026-06-03 / 2026-06-04
+ Outcome (PARTIAL).

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

## Walk-back options

- **If GPU kernel launch overhead exceeds eval savings for all practical circuits** -- fall back to CPU-only batched eval (WS1 still has value for cache locality). The gather/scatter refactor is useful regardless.
- **If OpenVAF GPU codegen proves infeasible** -- use OSDI-to-CUDA translation (Option B in plan WS2) as a pragmatic alternative.
- **If Rust FFI complexity becomes unmanageable** -- write GPU dispatch layer in pure C++ with CUDA/Metal APIs directly, avoiding the FFI boundary.

### What this does NOT change
- NR convergence algorithms (pnjlim, fetlim)
- Adaptive timestepping (LTE, predictor)
- Output format (rawfiles)
- Netlist parsing
- Any analysis type other than transient (initially)
