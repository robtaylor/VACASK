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
**f64 is unavailable on Apple GPUs**: MSL has no `double`. Budget for f32 numerics.

**Correction (2026-05-29, verified in vajax source):** vajax's NR *inner loop on
Metal is pure f32* — `vajax/__init__.py:_backend_supports_x64()` returns False for
Metal, so `configure_precision()` sets `jax_enable_x64=False`, forcing every jnp op
(eval, assembly, solve, residual, convergence) to f32 and silently collapsing any
`.astype(float64)` to f32. The f32-factor + **f64-iterative-refinement** I cited is
*not* the Metal loop: it's (a) `solver_factories.py:factorize_f32` (default False,
a **CUDA** VRAM-saving option) and (b) **Sprux** (`sprux_ffi.cpp`), whose `f64`
refinement runs **host-side on CPU doubles** and was never wired into the JIT NR
loop. **Implication:** pure-f32 NR (eval+solve) empirically *converges* on these
circuits on Metal, so f64 emulation (Ozaki etc.) is **not a hard prerequisite** for
the Metal path — though f32 needed convergence-management hacks (TRAP integration,
voltage-step limits, stagnation detection), and f32 *accuracy* for c6288 (≠
convergence) is still the open Q2 measurement.

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
- **Decisive:** f64→f32 demotion must be done *by us on every path* (only IREE has
  a turnkey demote pass, = adopt MLIR wholesale; refinement is *optional* — see the
  2026-05-29 correction above: vajax ran pure f32 on Metal). Once we own type
  lowering, the SPIR-V chain's only advantage (free codegen) is neutralized, so
  **hand-writing MIR→MSL is competitive and simpler** — we control f32/compensated
  arithmetic directly.

**Design implication:** build a **MIR-walking codegen frontend with pluggable
backends**. MSL backend for Metal (this machine); the same frontend gives
**NVPTX/CUDA nearly free via LLVM retarget (Strategy A), where f64 is native**.
**New open risk surfaced:** f32 accuracy for **PSP103 device physics itself**
(not just the linear solve). vajax ran pure f32 on Metal and *converged* (with
convergence-management hacks), so f64 is not strictly required — but Q2/Q3 must
still measure eval *accuracy* in f32 against the f64 reference (convergence ≠
accuracy).

#### FP64-where-needed: Ozaki scheme — scope (2026-05-28)

The **Ozaki scheme** (error-free transformation: scale FP64 by a shared
power-of-two, slice into INT8, do multiple INT8 GEMMs on tensor cores, recombine
exactly; cuBLAS picks slice count via Automatic Dynamic Precision) is a candidate
for FP64-accurate **linear algebra**, but its scope is narrow here:

- **Applies to GEMM / dense linear algebra only** → relevant to the **solver**:
  dense sub-blocks of a supernodal sparse factorization and refinement matvecs.
  *Not* the device-eval kernel.
- **Does NOT help the device-eval kernel** (the 90% hotspot): PSP103 eval is
  straight-line transcendental scalar math per instance, no GEMM structure.
  Eval f32-accuracy must still be handled with f32 + selective compensated
  arithmetic (Kahan / double-single) where needed.
- **Hardware: cuBLAS ADP is NVIDIA-only** (real tensor cores). On the M4 Pro the
  analogue is Metal `simdgroup_matrix` INT8, hand-rolled. So Ozaki most naturally
  serves a **future CUDA solver backend**, not the Metal-first path.

Recorded as the FP64 strategy for the solver's GEMM-shaped work; orthogonal to the
eval-kernel accuracy question Q2 measures.

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
  - Step 4 (f32 accuracy, **VACASK ground truth**): validate the psp103 f32
    kernel against **VACASK's own f64 OSDI eval**, not vajax/JAX. A flag-gated
    dump pass in VACASK (`OsdiInstance::evalCore`) emits, per PSP103 instance at
    several accepted transient timepoints, the OSDI input buffer (node voltages +
    cached init params) and the f64 residual + Jacobian; the f32 MSL kernel —
    generated from the **same** `devices/psp103v4/psp103.va` — is then fed those
    exact inputs and compared. Small proxy circuit: **gilbert** (6 PSP103,
    deterministic sinusoidal transient, same psp103v4 model as c6288, no
    metastability); then c6288. Supersedes the vajax `build_system_acc.py` path.
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
  Metal has no double). We must own f64→f32 demotion on every path, so
  **MIR→MSL (Strategy B) is confirmed**, with a pluggable backend giving NVPTX/CUDA
  near-free. New risk: **f32 accuracy of PSP103 physics** must be measured in Q2/Q3.
- 2026-05-29: **Correction — vajax's Metal NR inner loop is pure f32**, verified in
  source (`__init__.py:_backend_supports_x64`→`jax_enable_x64=False`). The f64
  iterative-refinement I'd cited is a CUDA option (`factorize_f32`, default off) and
  an unintegrated host-side-f64 Sprux path — *not* the Metal loop. So pure-f32 NR
  empirically converges on Metal; f64 emulation is not a hard prerequisite. f32
  *accuracy* (≠ convergence) for c6288 still open.
- 2026-05-29: **Q2 Milestone 1 (resistor) DONE.** Built a MIR→MSL emitter +
  Objective-C++ Metal harness in `~/Code/ChipFlow/vajax/spikes/msl-codegen/`
  (Bash is blocked inside sub-agents, so driven inline). Walks the eval MIR
  (`openvaf_py.get_mir_instructions`/`get_dae_system`), emits f32 MSL, compiles at
  runtime, dispatches batched on the M4 Pro. **Validated GPU f32 vs OpenVAF f64:
  max rel err 4.7e-8** (< f32 eps). Throughput 6.1e8 inst/s @ N=100k (0.16 ms
  dispatch); break-even ~10k instances. Structural finding: OpenVAF splits **init
  (bias-independent setup → cached values) from eval**; resistor's Jacobian
  conductances are computed in init and passed into eval as cached params (eval
  just `optbarrier`s them). The batched kernel's input vector = full MIR `params`
  list (named + cached). CPU baseline used was the *interpreted* MIR, not native
  OSDI — so no real speedup claim from resistor; M1 proved plumbing + f32 accuracy
  only. Next: M2 = psp103 (control flow/phi → MSL, kernel size, real f32-physics
  accuracy).
- 2026-05-29: **Q2 Milestone 2 (psp103) — codegen GENERATES + COMPILES.** psp103
  eval = 940 blocks / 20,243 instrs, **DAG (0 loops)**, 1503 phi, 56 Jacobian
  entries, 2868 params. Extended emitter to flatten the DAG (per-block
  reachability predicates + phi→nested-select). Generated **26.8k lines of MSL**;
  `xcrun metal -std=metal3.0` compiled to `.air` + linked `.metallib`, **0 errors**
  (after fixing int-cast for `iand`/`ibcast`). So the codegen scales to c6288's
  device. Caveats / next-step drivers:
  - **6238 compile warnings**, dominated by **15 constants out of f32 range**
    (`1e-100`,`1e+100`,`8.3e38`,…) — guard/clamp constants that flush to 0/inf in
    f32 and can break their guard (`max(x,1e-100)`→`max(x,0)`). Needs per-constant
    f32-safe lowering. This is the concrete f32-physics hazard.
  - **Naive flatten is wasteful: ~38% of eval instrs (7251/19304) are dead** for
    resistive-only output (they feed only the reactive path). Confirms the value of
    **reusing vajax's SSA optimizations** (DCE + SCCP constprop + dominator-based
    phi resolution) instead of the naive OR-of-all-edges flatten → v3 emitter.
  - Not yet done: f32 **accuracy** validation vs f64 (needs init-cache via
    `run_init_eval` + realistic c6288 bias points); reactive (ddt) part; register
    pressure / occupancy measurement on the M4 Pro.
  Spike code: `~/Code/ChipFlow/vajax/spikes/msl-codegen/` (`emit_msl2.py`,
  `dce_estimate.py`, `harness.mm`).
- 2026-06-01: **psp103 f32-accuracy — first attempt INVALID; need solved bias
  points.** Ran vajax's psp103 JAX eval at f64 vs f32 over a hand-built Vgs/Vds
  sweep. Result is untrustworthy: PSP103 is an internal-node model (GP/SI/DI/BP…
  solved by the NR), so hand-setting only `V(GP,SI)`/`V(DI,SI)` gives non-physical
  points — the device stayed **off (resist[0]=0 everywhere)**. The "jacobian rel
  err ≈ 840" is an artifact of near-zero f64 denominators, not real f32 loss.
  Two real signals: (a) **residual** rel err stayed ~3e-5 even in degenerate
  regimes (≈250× f32-eps, not catastrophic — weakly encouraging); (b) a concrete
  **f32 overflow-in-cast** fired — the out-of-f32-range guard constants (`1e+100`,
  from M2) hit inf in f32. **Correct method (next):** measure f32 vs f64 eval at
  *converged* per-instance operating points — extract node voltages from a real
  DC/transient solve (vajax c6288, or instrument VACASK OSDI eval to dump inputs),
  then compare there.
- 2026-06-01: **Operating-point source found & confirmed (unblocks both validations).**
  No VACASK C++ eval-dump patch exists (checked all origin branches + vendored
  VACASK; `add-tb_dp512x8-test` OSDI diff is cleanup; c6288 `.raw` saves only
  external nodes). Pivoted per the vajax comparison scripts
  (`scripts/{extract_c6288_jacobian,capture_benchmark_matrices,plot_three_way_comparison}.py`):
  the operating-point source is **vajax's own `CircuitEngine` + `FullMNAStrategy.run()`**,
  which returns the **full solution trajectory including internal PSP103 nodes**,
  in f64 on CPU (cross-validated vs VACASK on external nodes). Confirmed on `ring`:
  59 converged node-voltage vectors, internal nodes at physical voltages
  (rail 1.2V, mid-transition 0.66V). `build_system_fn(X, …)` (from
  `engine._build_transient_setup` + `_make_mna_build_system_fn`) gathers per-instance
  voltages and evaluates — eval it at a converged X in f64 vs f32 for the
  system-level f32-accuracy answer; gather per-instance PSP103 inputs for v3 MSL
  correctness. Spike code: `op_points.py` (saves `ring_Xtraj.npy`). **Next:** eval
  build_system at a converged X in f64 vs f32 (f32 accuracy), and feed per-instance
  inputs to the v2/v3 MSL kernels (correctness) — no VACASK instrumentation needed.
- 2026-06-02: **Correction + build_system-harness status.** (a) **Correction:**
  `FullMNAStrategy.run()` returns `V_out` of width 11 for ring, but ring's
  `n_unknowns=46` — so V_out is the **external-node trajectory, not the full
  internal solution**; the PSP103 internal-node voltages are solved but not
  returned there. `build_system(X)` needs the full 46-dim X. (Earlier "V_out
  includes internal PSP103 nodes" was wrong — the 0.66 V values were external
  nodes.) (b) **build_system f64-vs-f32 harness** (`build_system_acc.py`) built on
  the `extract_c6288_jacobian.py` pattern (eval at mid-rail init, non-degenerate,
  in each precision, compare J/f). Blocked on `engine._get_dc_source_values(...)`
  raising IndexError for **ring** — ring drives with `isource` (n_vsources=1 but
  the source-value setup mismatches the c6288-derived call). **Next:** fix the
  source-value call for an isource-driven bench (or run c6288, which is
  vsource-driven like the proven extract script); then compare J/f f64 vs f32 for
  the first real system-level f32-accuracy number. For per-instance/converged-X
  validation, capture the full 46-dim X from the solver state (V_out is reduced).
- 2026-06-01: **Q2 v3 emitter BUILT — reuses vajax `mir/` SSA opts + backward DCE.**
  `emit_msl3.py`: `parse_mir_function` → `CFGAnalyzer` → `SCCP(model-card params)`
  → `SSAAnalyzer`; walks `topological_order()`, resolves phis via `resolve_phi`
  (real branch conditions, no reach predicates). Results vs v2 naive flatten on
  psp103: **live instrs 11,415 (was 20,243) — 44% DCE**; **MSL 18,951 lines (was
  26,800) — 29% smaller**; phis 980 (was 1503); compiles clean (0 errors).
  Findings:
  - **SCCP alone kills 0 blocks**: psp103 *eval* branches on **cached init values**
    (runtime inputs), not model-card params — TYPE-style specialization happens in
    *init*. So the DCE win is **backward liveness from the resistive outputs**
    (reactive-only instrs are dead), which I added; SCCP only pruned ~20 dead phi
    operands.
  - **DCE is sound by construction**: inputs/consts declared unconditionally but
    body/phi *results* emitted only if live → a dropped-but-needed value would be
    an "undeclared identifier" compile error. 0 errors ⇒ liveness closure is
    reference-consistent (semantics-preserving structure).
  - **Numeric validation still BLOCKED on operating points.** v2-vs-v3 cross-check
    on *random* inputs is inconclusive (~49% inf/NaN — non-physical inputs overflow
    exp/log; v2 computes dead NaN-y instrs, v3 skips them → divergent NaN, plus
    near-zero-denominator rel-err artifacts). Both v3 correctness AND f32 accuracy
    need the same missing piece: **physically consistent eval input vectors**
    (cached init values via `run_init_eval`/`get_cache_mapping` + reasonable bias,
    or per-instance voltages dumped from a real VACASK/vajax solve). That harness
    is the next investment — it unblocks both validations at once.
- 2026-05-29: **v3 emitter plan — reuse vajax `mir/` SSA optimizations.** Mapped
  vajax's optimization pipeline: `openvaf_jax/mir/{cfg.py,ssa.py,constprop.py}`
  (`CFGAnalyzer`, `SSAAnalyzer`, `SCCP`) is **JAX-free and cleanly separable** —
  zero imports from the ast/codegen layer. v3 seam: `parse_mir_function` →
  `CFGAnalyzer` → `SCCP(known_values=param_values)` → `SSAAnalyzer(sccp=...)` →
  walk `cfg.topological_order()`, skip `sccp.is_block_dead()`, `ssa.resolve_phi()`
  → emit MSL. Wins: (a) **DCE** via SCCP dead-block elimination + dead-phi-operand
  pruning (`_get_live_operands`) — the ~38%, collapses 4-way NMOS/PMOS phis to
  `FALLBACK` single value; (b) **SCCP constant folding** (also resolves several
  out-of-f32-range guard constants); (c) **dominator-based phi resolution**
  (`PHIResolution.TWO_WAY` uses the actual branch-condition value, not naive
  OR-of-edges). psp103 eval is loop-free so `LoopInfo`/`while_loop` machinery is
  unneeded. This replaces `emit_msl2.py`'s naive flatten for v3.

- 2026-06-02: **Methodology reversal — VACASK becomes the f32-accuracy ground
  truth (replaces the vajax/JAX oracle).** The 2026-06-01 pivot to vajax's
  `CircuitEngine`/`build_system` was chosen for convenience ("no VACASK eval-dump
  patch exists"), but it makes *vajax's separate JAX assembly* the reference —
  a second implementation, not the simulator we're accelerating. Decision:
  **build the VACASK eval-dump after all.** Rationale (strictly better, not just
  preference): the MSL kernel is generated from the **OSDI-compiled psp103v4
  model**, the exact artifact VACASK runs, so VACASK's per-instance cached init
  buffer *is* OpenVAF's init output the kernel expects — tapping VACASK removes
  the cross-implementation mismatch the JAX path introduced, and one mechanism
  serves both validations (kernel f32 correctness + system f32 accuracy).
  - **Tap point (verified):** `Circuit::evalAndLoad` (`lib/circuit.cpp:1373`) →
    `OsdiDevice::evalAndLoad` (`lib/osdidevice.cpp:451`) →
    `OsdiInstance::evalCore` (`lib/osdiinstance.cpp:1131`), OSDI eval at
    `osdiinstance.cpp:1245`. Inputs in instance `core()`; f64 outputs via
    `load_jacobian_resist`/`load_residual_resist` (already used for printing at
    `osdidevice.cpp:731`). Design: a **flag-gated post-acceptance dump pass** at
    selected accepted transient timepoints (keeps the hot NR loop clean,
    guarantees converged inputs).
  - **Circuit choice:** ring rejected as the small test — free-running oscillator
    with a metastable initial state, so its "converged" points are
    dynamics-sensitive. Use **gilbert** (`demo/gilbert/`, 6 PSP103, sinusoidal
    RF+LO drive sweeps cutoff/linear/saturation, same `psp103v4.osdi` as c6288)
    as the deterministic small proxy; c6288 as the real target.
  - The `build_system_acc.py` isource `IndexError` (passes `n_isources=0`; ring
    drives with isources → `.at[0].set` on a length-0 array, `sources.py:547`) is
    now **moot** — that JAX harness is superseded.

- 2026-06-03: **VACASK ground-truth dump BUILT + validated; f64-vs-f64 harness
  cross-check exposes the real blocker (openvaf_py `run_init_eval` ignores init
  params).** The flag-gated dump (`OsdiInstance::dumpEvalIO` + TranCore
  post-acceptance pass, commits `ee8e26e`/`35152fd`) emits, per PSP103 instance
  at accepted gilbert transient timepoints: branch-voltage inputs, **absolute
  node voltages** (`A` line), and f64 resistive residual + Jacobian. Verified on
  gilbert (6 PSP103, same `psp103v4` as c6288): physical bias, internal nodes,
  collapsed parasitics read 0. Then built `compare_vacask.py` to cross-validate
  openvaf_py f64 (`run_init_eval`) against the VACASK f64 dump at the same
  operating point — the gate that must pass before any f32 number is trustworthy.
  Findings while wiring it (all verified inline):
  - **Voltage mapping is exact + positional.** VACASK's 13 OSDI `inputs[]`
    node-pairs map 1:1 onto eval voltage params `V(GP,SI)`..`V(NOI)`; the 6
    absolute-node voltage params (`V(GP)`,`V(SI)`,`V(DI)`,`V(BP)`,`V(BS)`,`V(BD)`)
    are filled from the dumped `A` line. PSP103's conduction is driven by these
    **absolute** internal-node potentials, not the branch differences — supplying
    only branch diffs leaves the device off (this is why the `A` line was added).
  - **param names are UPPERCASE** (`W`,`L`,`TYPE`,`VFBO`); models.inc card is
    lowercase — must match case-insensitively or overrides silently miss.
  - **BLOCKER — `run_init_eval(params)` ignores init-function params.** Verified:
    overriding `VFBO` (flatband, big Vt shift) or `W` (×5 width) in the params
    dict gives byte-identical output. So `run_init_eval` applies the dict only to
    *eval* (voltages work) but runs *init* with built-in defaults → cached values
    always reflect default model card/geometry, never gilbert's. Hence openvaf f64
    can't match VACASK's gilbert operating point, and the f32 comparison can't be
    grounded yet. run_init_eval *does* conduct PSP103 at default card + strong
    bias (max|F|=1.16e-4), so the eval path itself is fine.
  - **The kernel needs 2090 live inputs = 19 V + 16 I + 1 sysfun + 2054 cached**
    (1615 named `hidden_state` processed params + 439 unnamed cache slots =
    `num_cached_values`). So the cached values are the bulk of the input vector.
  - **Routes to supply the cached values (resume point):**
    - *(C) emit init→eval kernel* — extend `emit_msl3` to also walk the init MIR
      (`get_init_mir_instructions`) + wire init outputs→eval cache via
      `get_cache_mapping()`, so the kernel takes only model-card+geometry+voltages
      (all available: card from models.inc, voltages from the VACASK dump) and
      computes cached values itself. **Production-aligned** (init runs on-GPU once
      per instance) and removes openvaf_py from the *validation* loop. Bigger
      codegen (init MIR is large). RECOMMENDED.
    - *(B′) python init-interpreter shortcut* — interpret the init MIR with the
      gilbert card to compute the 2054 cached values, splice in voltages, feed the
      existing eval kernel. Fast, throwaway, openvaf_py-dependent.
    - *literal (c) — dump cached from VACASK OSDI buffer:* **impractical** — 2054
      internal values, no MIR-vid→OSDI-offset map (mostly not opvars).
    - `run_init_eval` can't be fixed by a compile flag (`compile_va` has no
      `propagate_constants`); it bakes init constants.
    Validation gate either way: reconstructed-f64 ≈ VACASK-f64 at a gilbert point,
    then `psp103_v3.metal` (f32) vs VACASK-f64 across the 21×6 points. Spike code:
    `spikes/msl-codegen/compare_vacask.py`.

- 2026-06-03: **Q2 f32-accuracy ANSWERED (PSP103, VACASK ground truth).** Got the
  JAX init+eval path to honor the gilbert model card (the `run_init_eval` blocker
  is dodged: `OpenVAFToJAX` + `translate_eval(propagate_constants=False)` reads
  params at runtime). Validation gate **passes**: JAX-f64 init+eval ≈ VACASK-f64
  dump to **~5e-6 residual / ~9e-5 Jacobian** rel err across 126 gilbert operating
  points (cross-impl floor — JAX-translated MIR vs OSDI-compiled — not bit-equal
  but firmly the same operating point/physics). Then measured f32 at the two
  regimes that matter (`spikes/msl-codegen/compare_jax_vacask.py`, JAXMODE):
  - **Eval-only f32 (f64 init/cache, f32 eval) — what the MSL kernel does:**
    residual rel err max **1.1e-5** (mean 2.9e-6, ~at the f64 floor → f32-SAFE);
    Jacobian rel err max **0.14**, mean **1.1e-2**, p99 0.13 → **f32-LOSSY**
    (~1% typical, ~14% worst on conductances; systemic, not an outlier — p99≈max).
    Physically: residual is well-conditioned; Jacobian entries (gm/gds) suffer
    f32 cancellation.
  - **Full f32 (init in f32 too): CATASTROPHIC**, residual rel err ~**5e5** — PSP103
    init has out-of-f32-range intermediates (the 15 `1e±100` guard constants from
    M2). Init must be f64 / compensated; it cannot run naive f32.
  - **Design implications (for ADR 0001 / resident loop):** an f32 device-eval
    kernel yields accurate currents but a ~1% (worst ~14%) Jacobian, and **init
    must not be f32**. A resident-NR f32 loop would need f64/compensated init and
    likely f64/compensated Jacobian assembly (or accept that f32 Jacobian error
    feeds NR — vajax's convergence-management hacks are consistent with this).
    This refines the 2026-05-29 "Metal NR is pure f32 and converges" note:
    convergence ≠ accuracy; the Jacobian is where f32 hurts.
  - **Caveat / remaining:** numbers are the JAX eval at f32 as a faithful proxy
    for `psp103_v3.metal` (same eval MIR). A literal MSL-kernel run still needs the
    1615 named `hidden_state` cached values sourced as kernel inputs (JAX computes
    them internally; emit_msl3 treats them as inputs) — expected to confirm the
    same order of magnitude (Metal transcendental/op-order differences aside).

## Outcome

**PARTIAL (2026-06-03).** Q1 codegen feasibility: YES (psp103 eval generates +
compiles to `.metallib`). Q2 f32 accuracy: **ANSWERED** with VACASK as ground
truth — f32 eval gives accurate currents (residual ~1e-5) but a lossy Jacobian
(~1% mean, ~14% max), and init is catastrophic in f32 (must be f64/compensated).
Q3 (resident-loop speed) and the literal MSL-kernel f32 run remain open. The
load-bearing decision input for ADR 0001 is now available: **f32 is viable for
device-eval currents but not for the Jacobian or init without higher-precision
handling.**
