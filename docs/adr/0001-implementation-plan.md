# GPU Device Evaluation — Implementation Plan

## Phase 0: Instrumentation & Baseline (1-2 weeks)

### 0.1 Per-step timing output
Add per-timestep timing breakdown to `TranCore::coroutine()` (gated behind
`tran_debug >= 2`). Use the existing `acctPrevPoint` slot to compute deltas.

Output per step: eval_time, factor_time, solve_time, nr_iters, accepted/rejected.

### 0.2 Device eval profiling
Enable `devacct = true` in `acct.h` and measure per-device-type eval times.
Need to understand: how much time is PSP103 vs resistor vs vsource?

### 0.3 Establish baseline numbers
Run all benchmarks (rc, graetz, ring, c6288, mul64 if available) with
instrumentation. Record per-step distributions, not just totals.

**Deliverable**: Baseline timing data showing where every microsecond goes.

---

## Phase 1: Abstract the Device Eval Interface (2-3 weeks)

### 1.1 Extract device eval from evalAndLoad
Currently `circuit.evalAndLoad()` iterates over all device instances
sequentially. Refactor to:

1. **Gather phase**: Collect input data (node voltages) for all instances of
   each model type into contiguous arrays
2. **Eval phase**: Call device eval function (currently OSDI) on the batch
3. **Scatter phase**: Write results (currents, charges, Jacobian entries)
   back to the MNA matrix

This is the same gather→eval→scatter pattern VAJAX uses.

### 1.2 Define the batched eval interface
```cpp
struct BatchedDeviceEval {
    // Input: node voltages for all instances (n_instances x n_terminals)
    // Input: model parameters (shared across instances)
    // Input: instance parameters (n_instances x n_instance_params)
    // Output: currents (n_instances x n_branches)
    // Output: charges (n_instances x n_branches)
    // Output: Jacobian entries (n_instances x n_jacobian_entries)
    
    virtual void eval(/* ... */) = 0;
};
```

### 1.3 CPU reference implementation
Implement the batched interface using existing OSDI calls. This must produce
bit-identical results to the current sequential path. Use this as the
correctness reference for GPU backends.

**Deliverable**: Refactored eval with gather/scatter, validated against current
output.

---

## Phase 2: GPU Backend — CUDA (3-4 weeks)

### 2.1 CUDA kernel for PSP103
Start with the highest-value target: PSP103 device eval. Options:

**Option A — OpenVAF CUDA codegen**: Extend OpenVAF compiler to emit CUDA
kernels from Verilog-A. Each kernel evaluates one instance; launch N threads
for N instances. This is the long-term right answer but significant compiler
work.

**Option B — Translate OSDI to CUDA**: The OSDI interface already has a
structured eval function. Write a tool that takes OSDI metadata (parameter
lists, node mappings) and generates a CUDA wrapper that calls into
device-specific eval logic. Less elegant but faster to prototype.

**Option C — CubeCL kernels**: Write the eval kernel in CubeCL (Rust),
getting CUDA+Metal+Vulkan from one source. Risk: CubeCL maturity for
scientific computing.

### 2.2 Memory management
- Allocate device buffers at elaboration time (circuit topology is static)
- Pin host memory for async transfers
- Double-buffer voltage vectors for overlap

### 2.3 Jacobian assembly on GPU
After device eval, Jacobian entries need to be scattered into the sparse
matrix. Options:
- Scatter on GPU (atomic adds to COO/CSR) — fast but needs careful handling
  of duplicate indices (VAJAX's bucketed scatter-add pattern)
- Download Jacobian entries, scatter on CPU — simpler, may be fast enough
  if eval dominates

### 2.4 Integration with NR loop
Wire the GPU eval path into `OpNRSolver::buildSystem()`:
- If GPU available and circuit exceeds threshold → GPU path
- Otherwise → existing CPU path
- Validate: GPU and CPU paths produce same NR convergence

**Deliverable**: c6288 running with GPU device eval on CUDA, matching CPU
results within tolerance.

---

## Phase 3: GPU Backend — Metal/Apple Silicon (2-3 weeks)

### 3.1 Metal compute shaders or WGPU
Leverage Sprux experience. Options:
- Direct Metal compute shaders (fastest, Apple-only)
- WGPU (portable, covers Metal+Vulkan+DX12)
- CubeCL (write once, compile to Metal/CUDA/Vulkan)

### 3.2 Unified memory advantage
On Apple Silicon, skip all the transfer management. GPU kernels read/write
the same memory as CPU. This simplifies Phase 2.3 significantly.

**Deliverable**: Same benchmarks running on Apple Silicon GPU.

---

## Phase 4: GPU Sparse Solver (future, optional)

Only pursue if profiling shows solver becoming the bottleneck after
device eval is on GPU. For c6288, solver is currently 10% of NR —
even with 10x eval speedup, solver becomes ~50% but is still only 4.5s.

Options:
- cuDSS (NVIDIA) — proven in VAJAX via spineax
- Sprux (Metal) — proven but JAX overhead killed it; native C++ would fix that
- CUDA graphs — record factor→solve sequence, replay without host sync

---

## Phase 5: Advanced Optimizations (future)

### 5.1 Multi-simulation batching
Run N Monte Carlo / corner simulations simultaneously on GPU. Each sim
has independent NR state. This multiplies GPU utilization without needing
larger circuits.

### 5.2 CUDA graphs for NR loop
Record the entire eval→assemble→solve→converge sequence as a CUDA graph.
Replay without host interaction. Eliminates kernel launch overhead.
Requires handling variable NR iteration count (conditional graph replay).

### 5.3 Device eval bypass on GPU
VACASK's bypass mechanism skips re-evaluating devices whose terminal
voltages haven't changed significantly. The GPU version could use a
predicate kernel to identify which instances need re-eval, then compact
the workload.

---

## Open Questions

1. **OpenVAF codegen strategy**: Should OpenVAF learn to emit GPU kernels
   directly, or should we build a separate VA→GPU compiler? The former is
   architecturally cleaner but requires deep changes to OpenVAF's backend.

2. **Which GPU first?** CUDA has better tooling and is the HPC standard.
   Metal/Apple Silicon has unified memory (simpler) and is our dev platform.
   Prototype on Metal for simplicity, optimize on CUDA for production?

3. **Build system**: CMake + Rust FFI adds complexity. Use corrosion
   (CMake-Rust integration)? Separate Rust workspace linked as static lib?

4. **Minimum viable circuit size for GPU**: VAJAX used 500 nodes. Need to
   benchmark kernel launch overhead to find the real crossover.

5. **Jacobian format**: VACASK uses KLU's compressed column format. GPU eval
   produces per-instance entries. What's the most efficient way to assemble?
   COO accumulate → CSC convert? Or write directly to CSC with known positions?
