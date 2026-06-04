# Handoff — GPU device-eval spike: Q2 literal-kernel-confirmed (currents AND resistive Jacobian f32-safe; init still f64), via VACASK ground-truth dump

**Created:** 2026-06-02 · **Last updated:** 2026-06-04 (literal-kernel session)
**Working tree:** doc updates to spike/ADR/handoff staged-or-uncommitted; untracked tooling (`.claude/`, `.tldr/`, `.tldrignore`). Prior session commits pushed to `gpu-acceleration` (through `661f431`). Spike code lives in the vajax tree (untracked scratch by convention).
**Branch:** gpu-acceleration

<!--
Ephemeral. At resolution, load-bearing pieces migrate to the spike, the plan,
or ADR 0001, then this file is git rm'd. See docs/handoff-discipline.md.
-->

## Goal & next-up

**Goal of this session:** Advance the GPU-accelerated device-eval work — finish
Phase 0 baseline (WS0), then run the OpenVAF-GPU-codegen spike: prove we can
generate Metal device-eval kernels from OpenVAF for c6288's PSP103, and start
validating f32 accuracy and the GPU-resident-NR-loop direction.

**Q2 f32 accuracy is ANSWERED + literal-kernel-confirmed (spike Findings 2026-06-03 /
2026-06-04 + Outcome).** With VACASK as validated ground truth: the *literal*
`psp103_v3.metal` kernel on the M4 Pro gives accurate currents (residual ~1.3e-5)
**and an accurate resistive Jacobian (~7.8e-5 max)** vs VACASK-f64 over 126 gilbert
op-points. The earlier proxy "Jacobian f32-lossy ~14%" was a **JAX-CPU-f32 artifact**
(corrected by the kernel run). Full f32 incl. init is still catastrophic (~5e5) →
init must be f64/compensated. Validation gate passed (JAX-f64 ≈ VACASK-f64 ~5e-6).

**Next session should pick up (Q2 tail now closed):**
1. ~~Literal MSL-kernel f32 run~~ **DONE** (2026-06-04). Bridge wires the 439 unnamed
   cache slots from the JAX f64 init cache via `get_cache_mapping()`
   (`(init_value_id, eval_param_position)`; `cache[i]` → `mir_func.params[eval_param]`),
   named hidden_state zeroed. Spike code: `spikes/msl-codegen/{build_kernel_inputs,
   compare_kernel_out,diag_kernel}.py`. Result revised Q2 (resistive Jacobian f32-safe).
2. **Q3 resident-loop slice** (the main remaining open spike question): eval+assembly+
   Sprux on Metal, ms/step vs c6288 CPU baseline (40.95 s NR). Depends on follow-ups
   2–3 below (per-instance/converged-X validation; reactive ddt).
3. **Reactive (ddt) Jacobian** — v3 emits resistive only; the f32-safe Jacobian result
   covers the resistive part. Add the ddt path before claiming transient-correct f32.
4. ~~Fold Q2 into ADR 0001~~ **DONE** + revised this session: ADR Q2 addendum now
   records the literal-kernel result (resistive Jacobian f32-safe; only init needs f64).
Harnesses: `spikes/msl-codegen/compare_jax_vacask.py` (the f32 measurement, JAXMODE
= f64|f32|f32eval), `compare_vacask.py` (run_init_eval cross-check, superseded).
Verified-good: voltage mapping (13 OSDI inputs → `V(GP,SI)`..`V(NOI)`, absolutes
from the `A` line), UPPERCASE param names, JAX init honors the card via
`OpenVAFToJAX`+`translate_eval(propagate_constants=False)`.

**Verification command:**

```sh
# 1. VACASK dump still produces the gilbert operating points (VACASK repo):
#    binary: ../build.VACASK/Release/simulator/vacask  (rebuild target `sim` if stale)
cd /tmp/claude/gilbert_dump
VACASK_EVAL_DUMP=/tmp/claude/gilbert_dump/eval_dump.txt VACASK_EVAL_DUMP_STRIDE=10 \
  ~/Code/build.VACASK/Release/simulator/vacask --skip-embed --skip-postprocess gilbert.sim
# expect: 767-line eval_dump.txt, 126 PSP103 INST blocks over 21 timepoints

# 2. The f32-accuracy measurement (spike code in the vajax tree):
cd ~/Code/ChipFlow/vajax/spikes/msl-codegen
PY=~/Code/ChipFlow/vajax/.venv/bin/python
JAXMODE=f64     $PY compare_jax_vacask.py /tmp/claude/gilbert_dump/eval_dump.txt 999  # gate: resid ~5e-6, jac ~9e-5
JAXMODE=f32eval $PY compare_jax_vacask.py /tmp/claude/gilbert_dump/eval_dump.txt 999  # answer: resid 1.1e-5, jac max 0.14
JAXMODE=f32     $PY compare_jax_vacask.py /tmp/claude/gilbert_dump/eval_dump.txt 999  # catastrophic ~5e5 (init in f32)
```
If `/tmp/claude/gilbert_dump/` is gone (it's `/tmp`): copy `demo/gilbert/{gilbert.sim,models.inc}`
there, strip the control block to just `elaborate circuit("sintest")` + `analysis tran1 tran
step=1n stop=4u maxstep=20n`. OSDI resolves from the staged `lib/vacask/mod/`.

## Done this session (literal-kernel, 2026-06-04)

Ran the **literal `psp103_v3.metal` kernel** on the M4 Pro vs the VACASK-f64 dump
(126 gilbert op-points) — the Q2 confirmation the prior session left open. Built
`build_kernel_inputs.py` (assembles the 2090-element `input_order` vector: voltages
from the dump, 439 unnamed cache slots from the JAX f64 init cache via
`get_cache_mapping()`, named hidden_state zeroed), `compare_kernel_out.py`,
`diag_kernel.py`; rebuilt `harness`. **Result revised Q2:** kernel-f32 vs VACASK-f64
= residual 1.3e-5, resistive Jacobian **7.8e-5** — the proxy's "14% lossy Jacobian"
was a JAX-CPU-f32 artifact (on the worst proxy entries, kernel ≈ f64 to ~1e-6).
Folded into the spike (Finding + Outcome 2026-06-04) and ADR 0001 (Q2 addendum
revised). Spike code is untracked vajax scratch; this repo's change is docs-only.

### Prior session (VACASK ground-truth dump, 2026-06-04 earlier — commits in `git log`)

Pivoted f32-accuracy validation from a vajax/JAX oracle to **VACASK as ground truth**,
built the VACASK f64 eval-dump, validated the operating-point reconstruction, and
**answered Q2** (proxy-level). Commits on `gpu-acceleration` (pushed):

| Commit | Subject | Notes |
|---|---|---|
| `8e19043` | docs: pivot f32-accuracy to VACASK ground truth | spike Q2 Step 4 + Finding; abandons vajax `build_system_acc.py` |
| `ee8e26e` | feat(spike): VACASK f64 device-eval dump | `EvalDumpSink`, `OsdiInstance::dumpEvalIO`, TranCore post-acceptance pass; flag-gated, hot path = 1 branch; `/simplify`-clean |
| `35152fd` | feat(spike): dump absolute node voltages (`A` line) | PSP103 bias needs absolute internal-node potentials, not just branch diffs |
| `be1f31d` | docs: refine resume routes | literal OSDI cache-dump impractical (2054 values); init→eval kernel recommended |
| `cfc2296` | docs: VACASK dump validated + run_init_eval blocker | recorded before the JAX-init workaround was found |
| `575cfee` | docs: **Q2 f32-accuracy ANSWERED** | currents f32-safe, Jacobian+init lossy; spike Outcome → PARTIAL |

**Headline result** (eval-only f32 = f64 init/cache + f32 eval, the MSL-kernel regime;
126 gilbert operating points vs VACASK-f64): residual rel err max **1.1e-5** (f32-safe);
Jacobian rel err mean **1.1e-2** / max **0.14** (f32-lossy); **full f32 incl. init =
catastrophic ~5e5** (out-of-f32-range init intermediates). Validation gate passed:
JAX-f64 ≈ VACASK-f64 to ~5e-6 residual / ~9e-5 Jacobian.

## Open follow-ups (priority-ordered)

### 1. Q2 f32-accuracy — DONE + literal-kernel-confirmed

Resolved across two sessions. Literal `psp103_v3.metal` kernel vs VACASK-f64 (126
gilbert op-points): currents f32-safe (residual ~1.3e-5) AND **resistive Jacobian
f32-safe (~7.8e-5 max)**; init still catastrophic in f32 (~5e5, must be
f64/compensated). The proxy's "Jacobian lossy ~14%" was a JAX-CPU-f32 artifact,
corrected here. See spike Findings 2026-06-03/2026-06-04 + Outcome (PARTIAL), and
the revised ADR 0001 Q2 addendum. Remaining: reactive (ddt) Jacobian + c6288 (see #3/#4).

### 1b. Push + fold Q2 into ADR 0001 — DONE

Commits pushed through `2968436`; Q2 result folded into ADR 0001's amendment as a
Q2-numerics addendum (commit `2968436`). **Still pending before any merge to main**
(not pre-push): `devacct=true` compiled in (`include/acct.h`) — revert to `false`
or runtime-gate; and the eval-dump scaffolding is removable spike code (gated by
`VACASK_EVAL_DUMP`; cf. `devacct`). These stay in place while the spike is Open.

### 2. Per-instance / converged-X validation (M)

`V_out` from `FullMNAStrategy.run` is the **reduced external-node trajectory**, not
the full `n_unknowns` internal solution (ring: 11 vs 46). To validate v3 MSL
correctness and per-device f32 accuracy at converged points, capture the **full
internal X** from the solver state (not V_out), gather per-PSP103 eval inputs, and
compare v2/v3 kernels + OpenVAF f64. (The random-input v2/v3 cross-check was
inconclusive — NaN on non-physical inputs.)

### 3. Reactive (ddt) part of eval (M)

v2/v3 emit only the **resistive** residual+Jacobian. PSP103's reactive (charge)
part comes via `TimeDerivative`/`NodeDerivative` calls (currently stubbed to 0).
Needed for a transient-correct kernel.

### 4. Resident-loop slice + Sprux solve (L)

The spike's M3: wire eval kernel + assembly + Sprux (Metal sparse) into a minimal
GPU-resident NR iteration for c6288; compare ms/step vs CPU baseline (c6288 NR =
40.95 s). Depends on 1–3.

### 5. Fold spike decisions into ADR 0001 (when spike resolves)

Per the spike's decision matrix — see ADR 0001 §"Amendment 2026-06-01".

## Critical context

- **VACASK eval-dump (this repo, committed):** `OsdiInstance::dumpEvalIO`
  (`lib/osdiinstance.cpp`) + TranCore post-acceptance pass (`lib/coretran.cpp`) +
  `EvalDumpSink` (`include/evaldump.h`). Env-gated: `VACASK_EVAL_DUMP=<path>`,
  `VACASK_EVAL_DUMP_STRIDE=N`. Emits per-PSP103-instance branch voltages,
  absolute node voltages (`A` line), f64 resistive residual + Jacobian at accepted
  tran timepoints. Removable spike scaffolding (cf. `devacct`).
- **Spike code is in the vajax tree**, not this repo: `~/Code/ChipFlow/vajax/spikes/msl-codegen/`
  (`emit_msl3.py`=v3 emitter, `harness.mm`=Metal runner, **`compare_jax_vacask.py`**=the
  f32 measurement [JAXMODE=f64|f32|f32eval], `compare_vacask.py`=run_init_eval cross-check
  [superseded], plus `emit_msl{,2}.py`, `op_points.py`, `xcheck_v2_v3.py`). Uses
  `~/Code/ChipFlow/vajax/.venv` (`openvaf_py`, `osdi_py` built). Whole dir is untracked
  scratch by convention — persists on disk, not committed to vajax.
- **Sub-agents cannot run Bash here** (sandbox) — all execute-and-verify work must
  run inline in the main session. Read-only investigation still delegates fine.
- **f32 accuracy is now MEASURED on the literal kernel** (was the open question):
  currents f32-safe (~1.3e-5) AND resistive Jacobian f32-safe (~7.8e-5) vs VACASK-f64;
  init catastrophic (~5e5). The proxy's "Jacobian lossy ~14%" was a JAX-CPU-f32
  artifact (full eval in f32, overflow-in-cast, no DCE), **not** a property of the
  generated MSL — corrected by the 2026-06-04 literal-kernel run. So a resident f32
  loop needs f64/compensated **init** but **not** f64 resistive-Jacobian assembly.
  Open: reactive (ddt) Jacobian (still stubbed) + c6288 (only gilbert measured).
  "Metal NR is pure f32 and converges" stands; convergence ≠ accuracy, but here
  the measured accuracy is good.
- **JAX-init reconstruction recipe** (the cache-computer; openvaf_py only, not
  shipped): `openvaf_py.compile_va(...)` → pick max-`num_jacobian` module →
  `OpenVAFToJAX(module)` → conftest `CompiledModel`. `run_init_eval(params)` IGNORES
  init params (no `propagate_constants` flag); `translate_eval(propagate_constants=False)`
  reads them at runtime. Param names are **UPPERCASE**; voltages map positionally
  (13 OSDI inputs → eval `V(GP,SI)`..`V(NOI)`; 6 absolutes from `A`). JAX x64 is
  sticky — force `jax.config.update("jax_enable_x64", ...)` AFTER all imports.
- **psp103 eval is a loop-free DAG** (940 blocks, 0 back-edges) → flattens to
  straight-line MSL. v3 reuses `vajax/openvaf_jax/mir/{cfg,ssa,constprop}.py`
  (JAX-free, separable). SCCP kills 0 blocks (eval branches on cached *init* values,
  not model params); the DCE win is **backward liveness from resistive outputs**.
- **15 out-of-f32-range guard constants** in psp103 (`1e-100`,`1e+100`) — `fconst`
  in emit_msl3 clamps them; watch these in accuracy work.
- **`devacct=true`** is compiled into VACASK (`include/acct.h`, commit `6d6378c`) —
  always-on overhead; revert to `false` (or runtime-gate) before any merge.

## References

- [`../spikes/openvaf-gpu-codegen.md`](../spikes/openvaf-gpu-codegen.md) — the living spike record (all findings, decision matrix)
- [`../plans/gpu-device-eval.md`](../plans/gpu-device-eval.md) — workstream state (WS0 done; WS1+ pending)
- [`../adr/0001-gpu-accelerated-device-evaluation.md`](../adr/0001-gpu-accelerated-device-evaluation.md) — decision + the 2026-06-01 amendment (4 choices under revision)

## Migration note

When the spike resolves:
- Follow-ups 1–4 outcomes → the spike's Findings/Outcome, then the validated
  decisions → ADR 0001 (fold the Amendment into Decision/Consequences) and the
  plan's WS1/WS2/WS3.
- "Spike code in vajax tree" + "sub-agents can't Bash" → if the codegen graduates
  into VACASK proper, a design doc `docs/gpu-codegen.md`.
- Then `git rm docs/handoffs/gpu-device-eval-handoff.md` in the migration commit.
