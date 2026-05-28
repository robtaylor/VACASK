# Spike — Can we generate a GPU-resident batched device-eval (and NR loop) for c6288 from OpenVAF, and does it beat the CPU NR?

**Status:** Open. Time-box: ≤ 1 week. Abort if Q1 (codegen path) has no
credible route by mid-week, or if the existing vajax-MLX c6288 run can't be
reproduced on this machine for a baseline.

<!--
Resolution states:
  Resolved (YYYY-MM-DD) — YES; <outcome>
  Resolved (YYYY-MM-DD) — NO; <outcome>
  Resolved (YYYY-MM-DD) — PARTIAL; <outcome>
-->

## Question

For c6288 (10k PSP103 transistors, 25k unknowns), can we generate GPU device-eval
kernels from OpenVAF and run a **GPU-resident Newton-Raphson loop** that beats the
CPU NR — on this hardware (Apple M4 Pro, Metal 4, **no CUDA**)?

## Why this is in question

This resolves **Open Question #1** in [the plan](../plans/gpu-device-eval.md) and
**revises a load-bearing decision in [ADR 0001](../adr/0001-gpu-accelerated-device-evaluation.md)**:

- ADR 0001 key design choice #1 is *"surgical replacement: only replace
  `evalAndLoad()`, keep the NR loop, timestep control, and convergence in C++."*
- The vajax experience concluded the opposite is optimal: **full GPU-resident
  kernels for the whole NR loop**, because the per-iteration host↔device sync —
  not compute — is the dominant cost. JAX/XLA's adaptive-timestep `lax.while_loop`
  forced a CPU round-trip every step (~80 ms/step) vs a fully-fused `fori_loop`
  (~0.2 ms/step). Offloading *only* eval re-incurs that round-trip every NR
  iteration (transfer sparse Jacobian out + solution back, ×3506 for c6288).

So the architectural fork is: **eval-only offload (CPU keeps NR + heuristics, pays
per-iter sync)** vs **resident NR loop (no per-iter sync, but NR control +
convergence heuristics must be ported to GPU)**. This spike picks the fork with
evidence.

## Established facts (2026-05-28)

### Codegen feasibility (from OpenVAF source `/tmp/claude/openvaf-llvm21` + IR dumps)

- **MIR is the clean tap point**, not OSDI. `CompiledModule` (`sim_back/src/lib.rs:245`)
  exposes `eval: Function` (with autodiff-expanded Jacobian values embedded),
  `intern: HirInterner` (Param → physical meaning), and `dae_system: DaeSystem`
  (sparse Jacobian structure). All before any LLVM/ABI decisions.
- **MIR is fully target-independent**: pointer-free SSA CFG over f64/i32/bool
  scalars; `mir/src/lib.rs:1-18` states codegen to hardware "is not a goal." 69
  opcodes, all scalar arithmetic + transcendentals. Pointer width enters only in
  `mir_llvm`.
- **Autodiff is forward-mode and operates on MIR** (`mir_autodiff/src/lib.rs:19`,
  `auto_diff(Function, ...) -> map of derivative Values`). Derivatives are MIR
  values → retargetable. **It is analytical/symbolic, not numeric AD.** Reuse it;
  do NOT re-autodiff (vajax: JAX AD gave `gds≈1e-16` in cutoff vs enforced `1e-9`
  floor → broke convergence on series-NMOS stacks).
- **Eval ABI** (`osdi/src/eval.rs:87`): `int eval(void* handle, void* instance,
  void* model, void* sim_info)`. The pure math body is GPU-friendly; the
  GPU-hostile parts are (a) struct-pointer GEP load/store for I/O — replace with
  flat per-instance device-buffer indexing; (b) flags-conditional stores — hoist
  / mask; (c) `$limit`/`$display`/`analysis` host callbacks — already conditional
  in MIR and no-op-able (`builder.rs:661`).
- **Scale**: psp103 → 63,603 lines MIR (510 KB OSDI). Large single kernel;
  register pressure is a real risk to measure. resistor → 844 lines (trivial).

### GPU-codegen strategy ranking (hardware-conditioned)

| Strategy | For CUDA | For Metal (this machine) | Effort |
|---|---|---|---|
| **A. Retarget LLVM module to GPU triple** | HIGH (nvptx64 + libdevice) | **N/A** — no upstream LLVM Metal/AIR target | Moderate |
| **B. Walk MIR → MSL / SPIR-V** | works | **the Metal path** | Higher (SSA→MSL has no phi; need SSA destruction) |
| C. Post-process compiled OSDI .so | LOW | LOW | — (ABI already baked) |

On the M4 Pro the route is **Strategy B → MSL** (or SPIR-V via Vulkan/MoltenVK).
**f64 is the hard wall**: MSL `double` is constrained on Apple GPUs → vajax used
**f32 factor + f64 iterative refinement** (Sprux pattern). Budget for f32 numerics.

#### LLVM→Metal trawl (2026-05-28) — why MIR→MSL, not the SPIR-V bridge

A GitHub/web trawl tested whether we could reuse OpenVAF's whole LLVM pipeline to
reach Metal (Strategy A for Metal) instead of hand-writing MIR→MSL:

- **No open LLVM→AIR backend exists.** Apple's Metal compiler is a closed LLVM
  fork emitting AIR; `gzorin/LLAIR` and `philipturner/llvm-metal` only shell out
  to Xcode's `metal` on MSL *source* or repackage `.metallib` — neither lowers
  arbitrary IR. Refuted.
- **The real bridge is LLVM IR → SPIR-V → SPIRV-Cross/naga → MSL.** All pieces
  mature for int/f32 OpenCL compute. The in-tree LLVM SPIR-V backend is official
  as of **LLVM 20**; our openvaf-r is on the **LLVM 21** branch, so version is
  fine.
- **It collapses at f64.** SPIRV-Cross `spirv_msl.cpp` hard-throws *"double types
  are not supported in buffers in MSL"*; Metal/MSL has no `double`; MoltenVK
  `shaderFloat64 = false`; fp64 emulation on Apple GPUs is documented-hard
  (`philipturner/metal-float64`, archived — compiler optimizes away double-single).
- **Decisive:** f64→f32 + iterative refinement must be done *by us on every path*
  (only IREE has a turnkey demote pass, = adopt MLIR wholesale). Once we own type
  lowering, the SPIR-V chain's only advantage (free codegen) is neutralized, so
  **hand-writing MIR→MSL is competitive and simpler** — we control f32/compensated
  arithmetic directly.

**Design implication:** build a **MIR-walking codegen frontend with pluggable
backends**. MSL backend for Metal (this machine); the same frontend gives
**NVPTX/CUDA nearly free via LLVM retarget (Strategy A), where f64 is native**.
**New open risk surfaced:** f32 accuracy for **PSP103 device physics itself**
(not just the linear solve) — vajax needed f64 residuals + refinement; Q2/Q3 must
measure eval accuracy in f32 against the CPU OSDI reference.

### Prior-art results (vajax — all GPU wins are Tesla T4 / CUDA; Metal has no published NR win)

| Circuit | Metric | GPU | CPU (JAX) | VACASK CPU | Note |
|---|---|---|---|---|---|
| c6288 | ms/step | **19.8** (T4, cuDSS) | 88–164 | 76–289 | ~3–5× over VACASK CPU, T4 only |
| mul64 | ms/step | **648** (T4, f32 cuDSS+refine) | 8325 | timeout | 12.8× over JAX-CPU |
| ring | ms/step | 1.49 | 0.51 | 0.109 | GPU **loses** small |
| rc | ms/step | 0.24 | 0.01 | 0.002 | GPU loses badly small |

- **GPU only pays above ~500–5000 nodes** (`gpu_threshold=500`). Below that,
  launch + vmap + COO-assembly overhead dominates.
- **The Metal killer is per-step host sync**: while_loop ~80 ms/step vs fused
  fori_loop ~0.2 ms/step (`vajax/docs/transient_step_dependencies.md:581`). This
  is the empirical basis for "go resident."
- **Sprux** (`vajax-mlx-perf/vajax/sprux/`): Metal sparse LU, f32 factor + 10-step
  f64 refine, double-buffered `beginSolve/endSolve` C++ API. Built but **not yet
  integrated into a JIT NR loop** — integration needs breaking the while_loop into
  explicit GPU-resident iterations. (For CUDA the equivalent is **cuDSS**, not
  cuDNN.)
- Metal needs **deterministic scatter-add** (f32 add non-associative under
  duplicate-index parallel scatter → bucketed scatter, `commit 5467202`).

## Approach (cheapest first)

- **Q1 — Empirical CPU-vs-GPU NR baseline on THIS machine.** Reproduce the
  existing vajax-MLX c6288 run on the M4 Pro. Gives a real Metal GPU-vs-CPU
  ms/step number today (even though it's the non-resident JAX/MLX version) and
  validates whether the resident-loop motivation holds locally.
  - Step 1: stand up `vajax-mlx-perf` env (uv, MLX, OpenVAF build).
  - Step 2: run c6288 (and ring as the small-circuit control), record ms/step,
    eval %, host-sync %.
- **Q2 — Proof-of-codegen: generate a batched MSL eval kernel for one device
  from MIR.** Start with `resistor` (844-line MIR) to prove the MIR→MSL walk +
  SSA destruction, then psp103 to expose register-pressure/f64 reality.
  - Step 1: emit MSL from `CompiledModule.eval` for resistor; validate one
    instance vs OSDI CPU eval (bit-ish within f32).
  - Step 2: batch over N instances; compare eval throughput vs CPU OSDI.
  - Step 3: attempt psp103; measure kernel size / register pressure / f64 impact.
- **Q3 — Resident-loop slice.** Wire Sprux solve + the eval kernel + assembly
  into a minimal GPU-resident NR iteration for c6288; sync only at accepted
  timepoints; compare ms/step vs VACASK CPU (40.95 s NR baseline).

## Decision matrix

| Outcome | Means | Action |
|---|---|---|
| Q1: GPU loses to CPU even at c6288 on Metal | f32/host-sync overhead defeats Metal here | Revisit: CUDA-only target, or shelve GPU on Apple; record in ADR 0001 |
| Q1 ok, Q2 fails (psp103 won't codegen/f64 unusable) | MSL path infeasible for big models | Pivot to SPIR-V/Vulkan or NVPTX-on-Linux; amend ADR scope |
| Q1 ok, Q2 ok, Q3 wins | Resident Metal NR viable | **Amend ADR 0001**: resident NR loop replaces "surgical eval-only"; promote WS4 (GPU sparse solver) from deferred to prerequisite |
| Q3 marginal | Resident helps but not enough | Hybrid: resident eval+assembly, CPU solve with batched sync; keep CPU heuristics |

## Findings

- 2026-05-28: OpenVAF MIR confirmed as clean, target-independent codegen source;
  forward-mode analytical autodiff on MIR. psp103 = 63.6k MIR lines.
- 2026-05-28: Strategy A (LLVM→GPU triple) is CUDA-only; this machine (M4 Pro) is
  Metal → Strategy B (MIR→MSL/SPIR-V). f64 is the hard wall.
- 2026-05-28: vajax CUDA c6288 ≈ 3–5× over VACASK CPU on T4; **Metal has no
  published NR win** — only bottleneck data (host-sync ~80 ms/step is the killer,
  motivating the resident-loop architecture).
- 2026-05-28: LLVM→Metal trawl — no open LLVM→AIR backend; the LLVM→SPIR-V→MSL
  bridge is mature but **hard-fails at f64** (SPIRV-Cross throws on double buffers;
  Metal has no double). We must own f64→f32 + refinement on every path, so
  **MIR→MSL (Strategy B) is confirmed**, with a pluggable backend giving NVPTX/CUDA
  near-free. New risk: **f32 accuracy of PSP103 physics** must be measured in Q2/Q3.

## Outcome

(Filled at resolution.)
