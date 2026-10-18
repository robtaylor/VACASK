# Decision 0001 — GPU-Accelerated Device Evaluation

**Status:** Accepted (partial) (2026-05-28). Ratified and partly built — WS0
instrumentation shipped; the GPU backend (WS1+) is not built. Which parts are live
is in the area doc's [Implementation status](../gpu-device-eval.md#implementation-status).

**Current architecture:** [GPU device evaluation](../gpu-device-eval.md).

> **Format note (2026-10-19):** this record predates the immutable-decision model
> and still carries dated `## Amendment` blocks below. They are kept (they are not
> misleading), and the present-tense current state now also lives in the
> [area doc](../gpu-device-eval.md). New changes follow the immutable rule — edit the
> area doc, or write a superseding decision; don't add further amendments here.
>
> The 2026-07-15 amendment resolved the four choices put under revision on
> 2026-06-01: **three withdrawn** (the surgical architecture is vindicated by
> measurement), **one folded in** (codegen is specifically MIR→MSL in f32).

## Amendment 2026-07-15 (spike resolved)

[Spike `openvaf-gpu-codegen`](spikes/openvaf-gpu-codegen.md) is **Resolved —
YES, hypothesis inverted**. It was built to justify replacing this decision's surgical
eval-only architecture with a GPU-resident NR loop. **It disproved its own premise.**

**The measurement that settles the fork (M4 Pro, 2026-07-15 spike Finding).** The
resident-loop case rested entirely on vajax's *"per-iteration host↔device sync, not
compute, dominates — `while_loop` ~80 ms/step"*. Measured directly on a C++/Metal
path with unified-memory buffers, a full `commit` + `waitUntilCompleted` round-trip
at c6288 scale costs **0.083 ms — roughly 1000× less**. The ~80 ms was **XLA's
`kWhile` thunk**: JAX framework overhead, not Metal, and VACASK's C++ NR loop has no
XLA in it to pay that cost. **This single number decides the fork** and nothing below
weakens it: it is a pure round-trip, independent of kernel size, scatter or run
length. Against the ~10 ms of CPU PSP103 work it would sit beside each iteration, a
per-iteration sync of 0.083 ms is **under 1%** — not an architecture driver.

**The speedup is a separate, much softer claim — read it with care.** The raw
resistive-eval kernel runs the 10,112-instance c6288 batch in **1.018 ms** (sync
included) vs the CPU's **10.30 ms** per NR iteration (fresh c6288 run 2026-07-15;
`tevalload`, 36.12 s / 3506 iters). That ratio is **~10×, but it is not an
eval-replacement number** and must not be quoted as one:

| | What it includes |
|---|---|
| CPU 10.30 ms/iter (`tevalload`, `lib/acct.cpp:25`) | gather + eval **resistive AND reactive** + **scatter into the KLU matrix** — the complete `evalAndLoad` |
| GPU 1.018 ms (best-of-20 dispatch) | kernel compute only, **resistive only**, inputs pre-staged, outputs to a flat buffer — **no gather, no scatter, no reactive** |

So the ~10× compares a *complete* CPU operation against a *partial* GPU one. The
honest eval-replacement figure (GPU gather + full resistive+reactive kernel + scatter
into KLU, vs CPU 10.30 ms) is **unmeasured and certainly well below 10×**. Treat
~10× as an optimistic upper bound on raw compute, pending the measurement ladder in
[the plan](../../plans/gpu-device-eval.md).

**Resolution of the four items** (all driven by the 0.083 ms sync fact, not the
speedup):

1. **Design choice #1 (surgical eval-only replacement): CONFIRMED — revision
   withdrawn.** The per-iteration sync that was supposed to make surgical offload
   untenable is noise (<1%). Going resident would buy back ≤0.083 ms/iter in
   exchange for porting the NR loop, convergence heuristics and timestepping to the
   GPU. The trade is not close. *Caveat:* this retires the **sync** motivation for
   going resident, not every possible one. If the Jacobian **scatter into KLU**
   proves expensive (unmeasured — see item #3 and the plan's WS1), keeping the
   Jacobian on-device with a GPU solver could be revisited — but that is an
   *assembly-cost* argument, not the sync argument this amendment was built on.
2. **"Does NOT change: NR convergence (pnjlim/fetlim), timestepping": CONFIRMED —
   revision withdrawn.** These stay in C++, as originally decided.
3. **Design choice #6 (solver stays on CPU; GPU solver not a prerequisite):
   CONFIRMED — revision withdrawn.** WS4 stays deferred: the solver is not a
   compute bottleneck (c6288 factor+solve ≈ 3.5 s absolute vs eval's 36 s). The one
   thing that could pull a GPU solver back in is *assembly*, not solve speed — if
   scatter into KLU is costly, an on-device solver would let the Jacobian skip the
   scatter entirely. That is contingent on the WS1 scatter measurement, not assumed.
4. **Design choice #2 (OpenVAF codegen: "extend or post-process") → FOLDED IN as
   MIR→MSL in f32.** This is the one revision that survives, and it is now settled
   rather than pending — see *Confirmed by the spike* and the *Q2 f32-numerics
   measurement* below. Codegen feasibility is proven (psp103 eval → `.metallib`) and
   f32 accuracy is measured on the literal kernel.

**Method lesson, recorded because it nearly cost the architecture.** The *same*
error — trusting a JAX/XLA proxy for a quantity on a substrate JAX was not running
on — produced **both** the phantom "14% lossy Jacobian" (corrected 2026-06-04)
**and** the phantom "80 ms sync" that put this decision's central design choice under
revision for six weeks. Both evaporated the moment the literal artifact was run on
the literal hardware, and the sync measurement that overturned item #1 took minutes
using tooling that already existed. Do not let a proxy move this decision again.
Compounding it: the spike's own *cheapest-first* step — an empirical Metal baseline
on this machine — was skipped in favour of expensive codegen. It was the step
designed to catch exactly this, and it would have done so on day one.

### Confirmed by the spike — supporting facts behind the Decision

These findings are settled and constrain any GPU path we build. They were true
independently of the surgical-vs-resident fork and remain so now it has resolved;
they are recorded here because each is a decision input that the Decision section
below now depends on.

- **OpenVAF MIR — not OSDI — is the codegen tap point.** `CompiledModule`
  (`sim_back/src/lib.rs:245`) exposes `eval` (autodiff-expanded), `intern`
  (Param → physical meaning) and `dae_system` (sparse structure) before any
  LLVM/ABI decision; MIR is pointer-free, target-independent SSA. This is the
  concrete form of design choice #2's "extend or post-process". Post-processing
  compiled OSDI (`.so`) is rejected — the ABI is already baked in.
- **Reuse OpenVAF's analytical autodiff; never re-differentiate.** It is
  forward-mode *symbolic* AD over MIR, so derivatives retarget for free. vajax
  re-autodiffed in JAX and got `gds≈1e-16` in cutoff against an enforced `1e-9`
  floor, breaking convergence on series-NMOS stacks. A GPU path that
  re-differentiates inherits that bug.
- **Apple GPUs have no f64 and no route restores it.** MSL has no `double`;
  SPIRV-Cross hard-throws on double buffers; MoltenVK reports
  `shaderFloat64=false`; fp64 emulation is documented-hard. Since f64→f32 demotion
  must be owned by us on *every* route, the LLVM→SPIR-V→MSL bridge's only advantage
  (free codegen) is neutralised — which is what makes hand-written MIR→MSL
  (Strategy B) the simpler choice rather than merely a fallback, with NVPTX/CUDA
  near-free from the same frontend.
- **init and eval are separate MIR functions, and the split is load-bearing.** init
  is bias-independent setup producing cached values; eval consumes them. Of psp103's
  2090 live kernel inputs, **2054 are cached init outputs**. Because init cannot run
  in f32 (see Q2 below) while eval can, this split is also **the precision
  boundary**: f64/compensated init, f32 eval.
- **Metal Jacobian assembly needs deterministic scatter-add.** f32 addition is
  non-associative, so a parallel scatter over duplicate indices is run-to-run
  nondeterministic; vajax fixed this with a bucketed scatter.
- **Ozaki/tensor-core FP64 emulation is solver-only and NVIDIA-only.** It is
  GEMM-shaped, so it cannot serve device eval (straight-line transcendental scalar
  math per instance). It is a candidate for a future CUDA *solver* backend only —
  not the Metal-first path.

### Metal viability: eval win demonstrated; end-to-end still unproven (2026-07-15)

Prior to 2026-07-15 **no Metal NR win had ever been demonstrated**, by us or in
vajax — every vajax GPU win is Tesla T4 / CUDA, and on small circuits Metal *loses*
outright (ring: 1.49 ms/step GPU vs 0.109 ms/step VACASK CPU). The decision-matrix
branch *"GPU loses to CPU even at c6288 on Metal → shelve GPU on Apple"* was
therefore live and untested.

**The raw resistive-eval kernel is fast on Metal (~10× the CPU's full eval+load as a
compute microbenchmark), so Apple Silicon is a viable target and the branch does not
fire.** But "viable at the eval layer" is *not* "proven end-to-end". What remains
unproven (these bound the win; none threatens the architecture):

- **The comparison is not apples-to-apples yet** (see the amendment's speedup table):
  the CPU number is full gather+eval+scatter, resistive+reactive; the GPU number is
  resistive-only compute with no gather/scatter. The real eval-replacement speedup is
  unmeasured and below 10×.
- **The reactive (ddt) path is stubbed.** v3's 44% dead-code win came from backward
  liveness against *resistive* outputs, which stripped much of the reactive path. A
  transient-correct kernel is larger and slower.
- **Jacobian scatter into KLU is unmeasured** and is the likeliest thing to eat the
  win. Gather is cheap by construction — only ~36 of each instance's 2090 inputs
  change per NR iteration (the other 2054 are bias-independent cached init values,
  staged once), so it is a ~1.5 MB unified-memory write. **Scatter is the risk**:
  566k f32 Jacobian entries into KLU's compressed-column structure, every iteration.
  Plan WS1 builds and measures this; WS3 should not commit to a backend before the
  number is in.
- **Only gilbert (6 devices) has been exercised for accuracy**; c6288-scale codegen
  is untested.
- **Small circuits remain CPU territory** — consistent with design choice #5's
  threshold, and the ring datapoint above shows the penalty is real.

### Wall-clock envelope: eval offload helps only on long-enough runs (2026-07-15)

A fresh instrumented c6288 run (`print stats`, same `tran stop=2n step=2p` benchmark
as WS0.3) decomposes the **144.3 s wall** as:

| Phase | Time | % of wall | Scales with run length? |
|---|---|---|---|
| Parse | 0.004 s | ~0% | no (one-time) |
| **Elaborate** | **103.1 s** | **71%** | **no (one-time)** |
| NR solve | 39.9 s | 28% | **yes** |
| — PSP103 eval+load | 36.1 s | 25% | yes |
| — refactor + solve | 3.5 s | 2% | yes |

This resolves an apparent contradiction in this decision: the Context table's *"NR = 99.5%
of total"* (decision-time profiling) vs WS0.3's *"NR = 28% of wall"*. **Both are
right at different run lengths.** This benchmark is a 2 ns / ~1000-step toy, so the
fixed 103 s elaboration dominates and NR is only 28%. On this run, even an
*infinitely fast* eval caps wall speedup at **144.3 / (144.3 − 36.1) = 1.33×**.

Elaboration is one-time, so a realistic transient (many ns → ×10–×1000 the steps)
amortises it toward zero; NR then returns to >90% of wall and eval offload approaches
its asymptotic **~5× wall** (eval ≈ 90% of NR, ~10× eval compute → NR ≈ 0.19×).
**Consequence for the project: the GPU-eval win is real only for runs long enough to
amortise elaboration.** The WS0.3 benchmark is the wrong length to demonstrate it and
should not be used as the wall-clock yardstick — the measurement ladder in the plan
fixes this.

### New cost center: elaboration (103 s) is unaccelerated and unaddressed (2026-07-15)

Elaboration of c6288 is **103 s — ~10.2 ms per PSP103 instance × 10,112** (RSS only
138 MB, so not memory-bound). This is almost certainly dominated by per-instance
PSP103 `init` — the bias-independent f64 setup that produces the 2054 cached values —
rather than KLU symbolic analysis of the 25k matrix *(? inferred: `telab` has no
sub-breakdown; confirming needs a timer inside elaborate)*. Two implications:

1. For short/medium runs it is now the **single largest cost**, and Decision 0001's
   surgical eval offload does nothing for it.
2. `init` is embarrassingly parallel across instances — the *same* property that
   motivates GPU eval — and is the f64-required half of the init/eval split. A future
   GPU (or just parallelised CPU) init is a candidate optimisation, but it is **out
   of scope for this decision** and needs its own decision. Flagged here so it is not lost.

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

**Settles amendment item #4's "Metal runs pure f32".** A generated f32 device-eval
kernel is accurate for both currents AND the resistive Jacobian (~1e-5 / ~8e-5 vs
f64); only **init** needs f64/compensated, and the init/eval MIR split is exactly
where that boundary falls. So the GPU path needs **no f64 emulation and no
compensated Jacobian assembly** for the resistive part — init runs f64 on the host
(once per instance, bias-independent), eval runs f32 on the GPU. (Convergence ≠
accuracy still holds as a general caution; here the measured accuracy is good.)
Scope of the measurement: gilbert, resistive Jacobian only — the reactive (ddt)
Jacobian (still stubbed) and c6288 are not yet covered; see *Metal viability* above.
Source: spike Findings 2026-06-03 / 2026-06-04 + Outcome.

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

> **Note (2026-07-15):** "NR % of total" here is the *steady-state* figure, with
> one-time elaboration amortised. On a short run elaboration dominates wall and NR is
> far lower (the WS0.3 2 ns benchmark: NR = 28% of wall, elaboration = 71%). See
> *Wall-clock envelope* in the amendment — this is why eval offload only pays on runs
> long enough to amortise the 103 s c6288 elaboration.

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

2. **OpenVAF MIR as GPU codegen source** (settled 2026-07-15): walk OpenVAF's
   **MIR** — not OSDI, not the compiled LLVM module — and emit GPU kernels for
   batched device evaluation, via a MIR-walking frontend with pluggable backends:
   **MSL for Metal** today, NVPTX/CUDA near-free from the same frontend. **Eval runs
   f32 on the GPU; init runs f64 on the host** (bias-independent, once per instance).
   Reuse OpenVAF's analytical autodiff; never re-differentiate. See the amendment
   and the *Confirmed by the spike* / *Q2 f32-numerics* sections above for why MIR,
   why MSL, and why that precision split.

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
- **Resistive-eval kernel compute is ~10× the CPU's full eval+load** (Metal, M4 Pro,
  sync included) — an *optimistic upper bound*, not an end-to-end result. The real
  eval-replacement speedup (full resistive+reactive, gather+scatter included) is
  unmeasured and below this; and it only converts to wall-clock on runs long enough
  to amortise elaboration (see *Wall-clock envelope*). Supersedes the decision-time
  "10-50× potential" estimate by grounding the upper bound in a measurement.
- Reuses VACASK's proven simulation infrastructure — confirmed, not just hoped: the
  NR loop, pnjlim/fetlim and timestepping all stay in C++.
- Incremental: can ship CPU-only and add GPU backends over time
- OpenVAF GPU codegen benefits the broader OpenVAF ecosystem

### Risks
- ~~OpenVAF GPU codegen is novel engineering (no prior art for VA→GPU)~~ —
  **retired**: psp103 eval generates from MIR and compiles to a `.metallib`, and the
  f32 result is validated against VACASK-f64 ground truth.
- **Jacobian scatter into KLU's compressed-column structure** (566k f32/iteration on
  c6288) is the main unquantified cost and the likeliest thing to erode the eval win.
  Measure it in WS1 before WS3 commits to a backend.
- **Wall-clock vs NR-time gap — triaged 2026-07-15:** the ~104 s between c6288 wall
  (144 s) and NR (40 s) is **one-time elaboration (103 s), not per-step overhead**
  (confirmed by `print stats`). It amortises over run length, so it does *not* cap a
  realistic long run — but it does mean the WS0.3 benchmark itself tops out at ~1.33×
  wall, and it surfaces elaboration as a large separate cost center (see the
  amendment). Not a risk to the architecture; a constraint on which runs show benefit.
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
