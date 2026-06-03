# Handoff — GPU device-eval spike: OpenVAF→Metal codegen validated; f32-accuracy + correctness blocked on operating points

**Created:** 2026-06-02
**Working tree:** clean (untracked `.claude/`, `.tldr/`, `.tldrignore` are tooling, not work)
**Branch:** gpu-acceleration (pushed to origin)

<!--
Ephemeral. At resolution, load-bearing pieces migrate to the spike, the plan,
or ADR 0001, then this file is git rm'd. See docs/handoff-discipline.md.
-->

## Goal & next-up

**Goal of this session:** Advance the GPU-accelerated device-eval work — finish
Phase 0 baseline (WS0), then run the OpenVAF-GPU-codegen spike: prove we can
generate Metal device-eval kernels from OpenVAF for c6288's PSP103, and start
validating f32 accuracy and the GPU-resident-NR-loop direction.

**Next session should pick up:** **unblock the openvaf-f64 reconstruction so it
matches VACASK-f64 at a gilbert operating point** (see spike Finding 2026-06-03).
The VACASK ground-truth dump is DONE + validated (commits `ee8e26e`, `35152fd`):
gilbert tran emits per-PSP103-instance branch voltages, absolute node voltages
(`A` line), and f64 residual/Jacobian at accepted timepoints. The cross-check
harness `spikes/msl-codegen/compare_vacask.py` is built; the blocker is that
`openvaf_py.run_init_eval(params)` **ignores init-function params** (verified:
overriding VFBO/W changes nothing), so cached values use default card/geometry,
not gilbert's → openvaf-f64 won't match VACASK-f64. Resume options (spike Finding):
1. Find an openvaf_py init-param API, or use `get_cache_mapping()` + a manual init
   pass to compute cached values from the gilbert card; OR dump cached eval inputs
   directly from VACASK's OSDI instance buffer.
2. Once openvaf-f64 ≈ VACASK-f64 at a gilbert point (the validation gate), run
   `psp103_v3.metal` (f32) on the same input vector → report f32-vs-VACASK-f64
   rel err across the 21×6 sampled points. This is the spike's open Q2.
Verified-good so far: voltage mapping (13 OSDI inputs → `V(GP,SI)`..`V(NOI)`;
absolutes from `A`), param case-sensitivity (UPPERCASE), run_init_eval conducts
PSP103 at default card (max|F|=1.16e-4).

**Verification command:**

```sh
# Spike code lives in the vajax tree, NOT the VACASK repo:
cd ~/Code/ChipFlow/vajax/spikes/msl-codegen
PY=~/Code/ChipFlow/vajax/.venv/bin/python

# v3 emitter generates psp103 MSL that the Metal compiler accepts (the M2 result):
$PY emit_msl3.py                                   # expect live_instr ~11415, MSL ~18951 lines
xcrun -sdk macosx metal -std=metal3.0 -c psp103_v3.metal -o /tmp/psp103_v3.air   # expect exit 0

# M1 resistor kernel validated vs OpenVAF f64 (rebuild + rerun harness):
clang++ -std=c++17 -O2 -fobjc-arc harness.mm -o harness -framework Metal -framework Foundation -framework QuartzCore
$PY bench.py                                       # expect max rel err ~4.7e-8

# Operating-point source runs (returns external-node trajectory):
$PY op_points.py                                   # expect ring V_out (timepoints, 11), 59 filled rows
```

## Done this session

| Commit | Subject | Notes |
|---|---|---|
| `6d6378c` | Per-step transient timing instrumentation (Phase 0) | WS0.1+0.2; `devacct=true` still compiled in — **revert before merge** |
| `d5b9b30` | WS0.3 baseline across all benchmarks | c6288 eval = 90.6% of NR; eval dominance tracks device-model complexity |
| `29cf6ed` | Add GPU codegen spike | resident-loop direction; revises ADR 0001 |
| `9f82192` | LLVM→Metal trawl | MIR→MSL confirmed (SPIR-V bridge hard-fails at f64) |
| `3fd78ed` | Q2 M1 resistor MIR→MSL | validated GPU f32 vs f64 = 4.7e-8 rel err |
| `9ef064b` | Q2 M2 psp103 generates+compiles | 940 blocks/20k instrs → `.metallib`, 0 errors |
| `578e1bb` | Correct f32/f64 claim | **vajax Metal NR loop is pure f32**; f64 emulation not a prerequisite |
| `b10e7b3` | Invalid psp103 f32-accuracy attempt | standalone bias non-physical; need solved operating points |
| `40ab79c` | v3 emitter (vajax SSA opts + backward DCE) | psp103 −44% live instrs, −29% MSL, compiles; DCE sound by compile-consistency |
| `4c6be85` | Operating-point source (vajax CircuitEngine) | `FullMNAStrategy.run` |
| `6dfc0f8` | ADR 0001 amendment | flags 4 choices under spike revision |
| `7ce5940` | Correct V_out claim + harness status | V_out is external-node trajectory (not full internal X) |

## Open follow-ups (priority-ordered)

### 1. f32-accuracy via VACASK ground truth (M, dump DONE; blocked on init-cache)

VACASK dump DONE + validated (commits `ee8e26e`, `35152fd`). Remaining: make the
openvaf-f64 reconstruction (`compare_vacask.py`) match VACASK-f64 at a gilbert
point — blocked because `run_init_eval` ignores init params (cached values use
default card, not gilbert's). See spike Finding 2026-06-03 for the precise
diagnosis + resume options. Then run `psp103_v3.metal` f32 and report rel err.
The open Q2 question in the spike.

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

- **Spike code is in the vajax tree**, not this repo: `~/Code/ChipFlow/vajax/spikes/msl-codegen/`
  (`emit_msl.py`=M1 resistor, `emit_msl2.py`=v2 naive flatten, `emit_msl3.py`=v3
  with vajax SSA opts, `harness.mm`=Metal runner, `op_points.py`, `build_system_acc.py`,
  `psp103_acc.py`, `xcheck_v2_v3.py`). It uses `~/Code/ChipFlow/vajax/.venv` (has
  `openvaf_py`, `osdi_py` built).
- **Sub-agents cannot run Bash here** (sandbox) — all execute-and-verify work must
  run inline in the main session. Read-only investigation still delegates fine.
- **Metal NR loop is pure f32** (verified: `vajax/__init__.py:_backend_supports_x64`
  → `jax_enable_x64=False`). f64 iterative refinement is CUDA-only / host-side Sprux,
  not the Metal loop. So **f64 emulation (Ozaki etc.) is NOT a prerequisite** — but
  f32 *accuracy* (≠ convergence) is still unmeasured at real operating points.
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
