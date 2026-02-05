# VACASK Simulation Options Reference

Comprehensive reference for all VACASK simulation options and analysis parameters.
Based on source code analysis of VACASK (options.cpp, coretran.cpp, coreop.cpp, etc.).

## Table of Contents

- [Global Options](#global-options)
  - [Temperature & Scaling](#temperature--scaling)
  - [Circuit Conductance](#circuit-conductance)
  - [Tolerances](#tolerances)
  - [Tolerance Reference Modes](#tolerance-reference-modes)
  - [Matrix & Solution Checking](#matrix--solution-checking)
  - [Output & File Options](#output--file-options)
  - [Sweep & Cosimulation](#sweep--cosimulation)
- [Newton-Raphson Solver](#newton-raphson-solver)
- [Operating Point Analysis](#operating-point-analysis)
  - [OP Options](#op-options)
  - [OP Analysis Parameters](#op-analysis-parameters)
  - [Homotopy](#homotopy)
- [Transient Analysis](#transient-analysis)
  - [Transient Options](#transient-options)
  - [Transient Analysis Parameters](#transient-analysis-parameters)
  - [Integration Methods](#integration-methods)
  - [LTE Calculation](#lte-calculation)
  - [Timestep Control](#timestep-control)
  - [Discontinuity Handling](#discontinuity-handling)
- [Small-Signal Analysis](#small-signal-analysis)
  - [AC Analysis Parameters](#ac-analysis-parameters)
  - [Noise Analysis Parameters](#noise-analysis-parameters)
  - [DC Transfer Function (DCxf) Parameters](#dc-transfer-function-dcxf-parameters)
  - [DC Incremental (DCinc) Parameters](#dc-incremental-dcinc-parameters)
  - [AC Transfer Function (ACxf) Parameters](#ac-transfer-function-acxf-parameters)
- [Harmonic Balance Analysis](#harmonic-balance-analysis)
  - [HB Options](#hb-options)
  - [HB Analysis Parameters](#hb-analysis-parameters)
- [Device-Level Limiting](#device-level-limiting)
- [Benchmark Examples](#benchmark-examples)

---

## Global Options

Options are set with the `options` keyword in `.sim` files:

```
options reltol=1e-3 vntol=1e-6 tran_method="gear2"
```

### Temperature & Scaling

| Option | Default | Description |
|--------|---------|-------------|
| `temp` | 27 | Ambient simulation temperature (°C). Mapped to `$temperature` in Verilog-A. |
| `tnom` | 27 | Device parameter measurement temperature (°C). Mapped to `$simparam("tnom")`. |
| `scale` | 1.0 | Global circuit scaling factor (>0). Applied to device geometry parameters. |

Both `temp` and `scale` affect parameterized expressions. Changing them triggers
circuit re-parameterization.

#### How Temperature Flows Through the Simulator

1. The user sets `temp` in °C (default 27°C).

2. During device setup, VACASK converts to Kelvin for the OSDI interface
   (`osdidevice.cpp:302`):
   ```cpp
   model->setupWrapper(circuit, sp, opt.temp + 273.15, ...);
   ```
   This value is accessible in Verilog-A as `$temperature` (in Kelvin).

3. `tnom` is passed through as-is in °C via `$simparam("tnom")`
   (`osdidevice.cpp:524-526`). Devices use this to compute parameter adjustments
   for temperature scaling (e.g., mobility vs temperature curves measured at
   `tnom`).

4. Inside device models, temperature-dependent quantities are computed:
   ```
   Vt = k × T / q ≈ 0.02585V at 300K (27°C)
   ```
   Where `k` = Boltzmann constant, `T` = absolute temperature, `q` = electron
   charge. This thermal voltage `Vt` appears in:
   - PN junction current: `I = Is × (exp(V / (n×Vt)) - 1)`
   - Device limiting: `DEVpnjlim` uses device-specific `Vt` and `Vcrit`
   - Saturation current temperature scaling: `Is(T) = Is(Tnom) × f(T, Tnom)`

5. All simulator options (`reltol`, `vntol`, `abstol`, `gmin`, etc.) are passed to
   devices as `$simparam` values, allowing models to query them at runtime.

The full set of `$simparam` values available to devices:
`iniLim`, `gmin`, `gdev`, `tnom`, `minr`, `scale`, `iteration`,
`simulatorVersion`, `simulatorSubversion`, `sourceScaleFactor`,
`reltol`, `vntol`, `abstol`, `chgtol`, `fluxtol`.

String parameters: `analysis_name`, `analysis_type`, `cwd`.

### Circuit Conductance

| Option | Default | Description |
|--------|---------|-------------|
| `gmin` | 1e-12 | Minimum shunt conductance (S) added to device nodes (≥0). Helps convergence by preventing singular matrices. |
| `gshunt` | 0.0 | Shunt conductance (S) from all potential nodes to ground (0=off, >0 active). |
| `minr` | 0.0 | Minimum series resistance (Ω) (≥0). |

`gmin` is added to device terminals and is the value manipulated during gmin-stepping
homotopy. `gshunt` adds a parallel conductance from every potential node to ground,
independent of gmin.

### Tolerances

| Option | Default | Unit | Description |
|--------|---------|------|-------------|
| `tolmode` | `"spice"` | — | Tolerance assignment method: `"spice"`, `"va"`, `"mixed"` |
| `tolscale` | 1.0 | — | Global scaling factor for all absolute tolerances (>0) |
| `reltol` | 1e-3 | — | Relative tolerance (0 < x < 1) |
| `abstol` | 1e-12 | A | Absolute current tolerance |
| `vntol` | 1e-6 | V | Absolute voltage tolerance |
| `chgtol` | 1e-15 | As | Charge tolerance (default = 1mV across 1pF) |
| `fluxtol` | 1e-14 | Vs | Flux/inductance tolerance (default = 1µA across 10nH) |

#### Tolerance Modes

- **`spice`**: SPICE-style assignment:
  - `vntol` for non-flow unknowns, `abstol` for flow unknowns
  - `abstol` for resistive residual of non-flow unknowns
  - `chgtol` for reactive residual of non-flow unknowns
  - `vntol` for resistive residual of flow unknowns
  - `fluxtol` for reactive residual of flow unknowns

- **`va`**: Use Verilog-A natures and disciplines where available; otherwise
  do not enforce tolerances.

- **`mixed`**: Use Verilog-A natures and disciplines where available; otherwise
  fall back to SPICE tolerances.

#### How Tolerances Are Used

Tolerances serve three roles in VACASK, each applied differently:

**1. NR Solution Delta Convergence** (`coreopnr.cpp:966-1050`)

At each NR iteration, the change in each unknown is checked:

```
tol = max(|tolref| × reltol, unknown_abstol[i])
converged if |delta[i]| < tol
```

Where `unknown_abstol[i]` is the per-unknown absolute tolerance assigned from
`vntol` or `abstol` depending on the nature of the unknown and the tolerance mode.
The `tolref` is determined by the `relrefsol` option (see Tolerance Reference Modes).

**2. Instance Residual Convergence** (`osdiinstance.cpp:896-908`)

Each device instance checks whether its residual has changed significantly:

```
tol = max(|f_prev| × reltol, residual_abstol[u]) × nr_convtol
converged if |f_new - f_prev| < tol
```

Here `residual_abstol[u]` is the per-equation tolerance (e.g., `abstol` for KCL
current-sum equations at potential nodes). The `nr_convtol` factor (default 0.01)
makes this check 100× stricter than the absolute tolerance alone.

**3. LTE Timestep Control** (`coretran.cpp:1168`)

LTE tolerance for each unknown:

```
tol = max(|tolref| × reltol, unknown_abstol[i]) × tran_lteratio
```

The `tran_lteratio` factor (default 3.5) loosens the tolerance for LTE, allowing
larger timesteps. The `tolref` is determined by the `relreflte` option.

#### Practical Impact of `abstol`

The default `abstol=1e-12` (1 pA) is the convergence floor for current quantities.
For circuits with very small currents (e.g., leakage-dominated states), this
determines when the solver considers currents "converged." Increasing `abstol`
loosens convergence (faster but less accurate); decreasing it tightens convergence
(slower but more accurate).

`tolscale` provides a convenient way to uniformly scale all absolute tolerances
without changing their relative proportions. Setting `tolscale=10` multiplies all
of `abstol`, `vntol`, `chgtol`, `fluxtol` by 10.

### Tolerance Reference Modes

These control how the reference value for relative tolerance is computed.

| Option | Default | Description |
|--------|---------|-------------|
| `relrefsol` | `"relref"` | Reference for solution delta relative tolerance |
| `relrefres` | `"relref"` | Reference for residual relative tolerance |
| `relreflte` | `"relref"` | Reference for LTE relative tolerance |
| `relref` | `"alllocal"` | Master reference mode when above are set to `"relref"` |

Each `relrefsol`, `relrefres`, `relreflte` can be set independently to:

| Value | Description |
|-------|-------------|
| `"pointlocal"` | At current time, for each unknown separately |
| `"local"` | Maximum over past time, for each unknown separately |
| `"pointglobal"` | At current time, maximum over all unknowns |
| `"global"` | Maximum over past time, maximum over all unknowns |
| `"relref"` | Defer to the `relref` master option |

The `relref` master option values:

| Value | Description |
|-------|-------------|
| `"alllocal"` | Per-unknown, maximum over past time (default) |
| `"sigglobal"` | Signal-global reference |
| `"allglobal"` | Global maximum over all unknowns and past time |

### Matrix & Solution Checking

| Option | Default | Description |
|--------|---------|-------------|
| `matrixcheck` | 0 | Check matrix for inf/nan (0=off) |
| `rhscheck` | 1 | Check RHS vector for inf/nan (1=on) |
| `solutioncheck` | 1 | Check solution vector for inf/nan (1=on) |
| `rcondcheck` | 0 | Matrix condition number check. If >0, fails when rcond < value. Used in small-signal analyses. |

### Output & File Options

| Option | Default | Description |
|--------|---------|-------------|
| `rawfile` | `"binary"` | Output file format: `"ascii"` or `"binary"` |
| `strictoutput` | 2 | Output file handling on error: 0=leave files, 1=delete on first error, 2=delete before analysis |
| `strictsave` | 1 | Binding failure: 0=bind to 0, 1=error on first binding failure only, 2=always error |
| `strictforce` | 1 | Nodeset/IC conflict: 0=warn and ignore, 1=abort |
| `accounting` | 0 | Timing accounting: 0=outside analyses only, 1=inside and outside |

### Sweep & Cosimulation

| Option | Default | Description |
|--------|---------|-------------|
| `sweep_pointmarker` | 0 | Yield at each sweep point for cosimulation (1=enable). Stops analysis and returns `SweepPoint` state. |
| `sweep_debug` | 0 | Sweep debugging (1=debug, ≥2=print details) |

---

## Newton-Raphson Solver

### NR Options

| Option | Default | Description |
|--------|---------|-------------|
| `nr_damping` | 1.0 | Global NR damping factor (0 < x ≤ 1). Scales every Newton delta by this factor. 1.0 = no damping. |
| `nr_convtol` | 0.01 | Convergence tolerance factor for instance residuals. Should be <1. Multiplied against `abstol`/`vntol`. |
| `nr_bypasstol` | 0.01 | Bypass tolerance factor for instance input check. Should be <1. |
| `nr_bypass` | 0 | Enable instance bypass (1=on). Skip device re-evaluation when inputs barely changed. |
| `nr_contbypass` | 1 | Allow forced bypass in first NR iteration of continuation mode. |
| `nr_conviter` | 1 | Number of consecutive convergent iterations required before declaring convergence. |
| `nr_residualcheck` | 1 | Check residual (not just solution delta) for convergence (1=on). |
| `nr_force` | 1e5 | Forcing factor for nodesets and initial conditions. |
| `nr_debug` | 0 | Debug level: ≥1=messages, ≥2=new solution, ≥3=old solution, ≥4=linear system |

### System-Level Damping

VACASK applies damping as a static scalar on the Newton update (`nrsolver.cpp:363-365`):

```cpp
// Compute new solution with static damping
for (i = 1; i <= n; i++) {
    xnew[i] = xprev[i] - xdelta[i] * settings.dampingFactor;
}
```

The default `nr_damping=1.0` means no system-level damping is active. There is no
adaptive damping (line search, backtracking) at the system level. Convergence
assistance comes from device-level `$limit` callbacks (see
[Device-Level Limiting](#device-level-limiting)).

### Convergence Criteria

Convergence is checked at two levels:

1. **Instance level**: Each device instance must have its residuals below
   `nr_convtol × abstol` (or `vntol`). When `nr_bypass=1`, instances that converge
   and whose inputs change by less than `nr_bypasstol × tolerance` are skipped in
   subsequent iterations.

2. **System level**: Solution delta must be below tolerance, and (when
   `nr_residualcheck=1`) the global residual must also be below tolerance.

Convergence requires `nr_conviter` consecutive iterations meeting both criteria.

---

## Operating Point Analysis

### OP Options

| Option | Default | Description |
|--------|---------|-------------|
| `op_itl` | 100 | Max iterations in non-continuation mode |
| `op_itlcont` | 50 | Max iterations in continuation mode |
| `op_skipinitial` | 0 | Skip initial plain OP, go straight to homotopy (1=yes) |
| `op_homotopy` | `["gdev", "gshunt", "src"]` | Homotopy algorithms to try, in order |
| `op_srchomotopy` | `["gdev", "gshunt"]` | Fallback homotopy when source stepping fails at factor=0 |
| `op_nsiter` | 1 | Number of iterations during which nodesets are applied |
| `op_debug` | 0 | OP debugging: ≥1=iteration/homotopy/convergence, ≥2=continuation info |
| `q_debug` | 0 | Reactive residual debugging: ≥1=dump Q after OP, ≥2=at each NR iteration |

### OP Analysis Parameters

```
analysis op1 op [nodeset=<value>] [store=<slot>] [write=<0|1>]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nodeset` | — | Nodesets to guide convergence (string or list format) |
| `store` | — | Slot name to store solution for later use |
| `write` | 1 | Write results to output file |

### Homotopy

#### Algorithm Sequence

The OP solver follows this sequence (`coreop.cpp:407-495`):

1. **Plain OP** (unless `op_skipinitial=1`) — try direct Newton-Raphson
2. **If failed**, iterate through `op_homotopy` list — stop at first success

#### Available Algorithms

| Name | Strategy |
|------|----------|
| `gdev` | Add conductance to device terminals, reduce dynamically |
| `gshunt` | Add shunt conductance to ground, reduce dynamically |
| `spice3gmin` | SPICE3-style fixed-schedule gmin reduction |
| `src` | Scale all sources 0→1 with adaptive step sizing |
| `spice3src` | SPICE3-style fixed-schedule source ramping |

#### Gmin Stepping Options

| Option | Default | Description |
|--------|---------|-------------|
| `homotopy_gminsteps` | 100 | Max steps (≤0 disables gmin stepping) |
| `homotopy_gminfactor` | 10.0 | Initial gmin stepping factor |
| `homotopy_maxgminfactor` | 10.0 | Maximum gmin stepping factor |
| `homotopy_mingminfactor` | 1.00005 | Give up below this factor |
| `homotopy_startgmin` | 1e-3 | Starting gmin value |
| `homotopy_maxgmin` | 1e2 | Fail if can't solve above this gmin |
| `homotopy_mingmin` | 1e-15 | Target gmin value (stopping point) |

#### Source Stepping Options

| Option | Default | Description |
|--------|---------|-------------|
| `homotopy_srcsteps` | 100 | Max steps (≤0 disables source stepping) |
| `homotopy_srcstep` | 0.01 | Initial source step size |
| `homotopy_srcscale` | 3.0 | Multiply on success, divide on failure |
| `homotopy_minsrcstep` | 1e-7 | Give up below this step size |
| `homotopy_sourcefactor` | 1.0 | Manual homotopy scaling (sourcefactor × srcscalefactor). Allows DC sweep to control homotopy. |

The dynamic variants (`gdev`, `gshunt`, `src`) are adaptive — they scale the step
on success/failure. The `spice3*` variants use a fixed schedule.

---

## Transient Analysis

### Transient Options

#### Integration Method

| Option | Default | Description |
|--------|---------|-------------|
| `tran_method` | `"trap"` | Integration method (see [Integration Methods](#integration-methods)) |
| `tran_maxord` | 2 | Maximum integration order (ignored for fixed-order methods) |
| `tran_xmu` | 0.5 | Trapezoidal mixing: 0=pure Euler, 0.5=pure trapezoidal, 1.0=backward Euler |

#### Timestep Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `tran_fs` | 0.25 | Initial timestep fraction: actual first dt = `step × tran_fs` |
| `tran_ffmax` | 0.25 | Fraction of max excitation frequency period for initial timestep limit (0=off) |
| `tran_fbr` | 0.2501 | Breakpoint fraction: limits dt to `fbr × distance_to_breakpoint`. Guarantees ≥4 points between breakpoints. |
| `tran_rmax` | 0.0 | Ratio of step to timestep upper limit (<1 = no upper limit based on step) |
| `tran_minpts` | 50 | Minimum output points from start to stop (<1 = no limit). Caps max_dt to `(stop-start)/minpts`. |
| `tran_ft` | 0.25 | Timestep cut factor when NR fails (0 < x < 1). Timestep is multiplied by this. |

#### NR and LTE

| Option | Default | Description |
|--------|---------|-------------|
| `tran_itl` | 10 | Max NR iterations per transient timepoint |
| `tran_predictor` | 0 | Use polynomial predictor for initial guess (1=yes, 0=use previous solution) |
| `tran_lteratio` | 3.5 | LTE overestimation factor (>1). Higher = more permissive = larger timesteps. |
| `tran_redofactor` | 2.5 | Rejection threshold. If `dt / lte_dt > redofactor`, reject timepoint. 0=never reject. |
| `tran_spicelte` | 0 | LTE mode: 0=correct, 1=SPICE-like (incorrect but compatible) |
| `tran_trapltefilter` | 1 | Trap ringing filter for predictor and LTE (applies only to Adams-Moulton order 2) |

#### Debugging

| Option | Default | Description |
|--------|---------|-------------|
| `tran_debug` | 0 | 1=accept/reject steps, 2=LTE calculation details, 3=per-node LTE |

### Transient Analysis Parameters

```
analysis tran1 tran step=<dt> stop=<tstop> [start=<tstart>] [maxstep=<hmax>] [icmode="uic"|"op"] [ic=<value>]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `step` | — | Nominal timestep. Used as initial dt (scaled by `tran_fs`). |
| `stop` | — | Simulation end time. |
| `start` | 0.0 | Time at which to start recording output. Simulation still runs from t=0. |
| `maxstep` | — | Maximum allowed timestep. Overrides `tran_minpts` calculation if smaller. |
| `icmode` | (none) | Initial condition mode: `"op"` (run DC OP with IC forcing), `"uic"` (use ICs directly, skip OP) |
| `ic` | — | Initial conditions (string slot or list). |
| `store` | — | Slot name to store solution. |
| `write` | 1 | Write results to output file. |

### Integration Methods

| Method | Alias | Order | Description |
|--------|-------|-------|-------------|
| `am` | — | up to `tran_maxord` | Adams-Moulton (variable order) |
| `bdf` / `gear` | — | up to `tran_maxord` | Backward Differentiation Formula (variable order) |
| `euler` | — | 1 | Backward Euler (= AM order 1) |
| `trap` / `am2` | — | 2 | Trapezoidal (= AM order 2, default) |
| `gear2` / `bdf2` | — | 2 | Gear order 2 (= BDF order 2) |

**Trapezoidal** (`trap`): Adams-Moulton order 2. Default method. The `tran_xmu`
parameter controls mixing (0.5=pure trapezoidal). `tran_trapltefilter` enables
filtered history for LTE estimation to avoid trap oscillation. Known to exhibit
numerical ringing on stiff problems.

**Gear2** (`gear2`): BDF order 2. More stable for stiff problems (e.g., diode
switching). Slightly more dissipative than trapezoidal.

### LTE Calculation

#### Formula

The Local Truncation Error is computed per unknown at each accepted timepoint
(`coretran.cpp:1074-1191`):

```
factor = integErrCoeff / (integErrCoeff - predErrCoeff)
lte = factor × (solution[i] - predicted[i])
tol = max(|tolref| × reltol, abstol) × lteratio
ratio = |lte| / tol
```

Where:
- `integErrCoeff` = error coefficient of the integration method
- `predErrCoeff` = error coefficient of the polynomial predictor
- `factor` is always < 1 (typically 0.18 for Gear2 order 2)
- `tolref` depends on `relreflte` option (default: per-unknown maximum over past time)

#### Error Coefficient Values (Gear2)

| Property | Warmup (order 1) | Steady state (order 2) |
|----------|------------------|------------------------|
| `integErrCoeff` | 1.0 | 1.333 |
| `predErrCoeff` | -2.0 | -6.0 |
| `factor` | 0.333 | 0.182 |

#### Timestep Adjustment

```
hkFactor = ratio^(-1/(order+1))    // For order=2: ratio^(-1/3)
dt_new = min(2 × dt, dt × hkFactor)
```

- `ratio < 1` → `hkFactor > 1` → timestep grows (capped at 2×)
- `ratio > 1` → `hkFactor < 1` → timestep shrinks
- `ratio = 1` → timestep unchanged (steady state)

#### SPICE vs Correct LTE

When `tran_spicelte=1`, VACASK uses the SPICE3 formula which divides tolerance by
`(order+1)!` and uses `ratio^(-1/order)` instead of `ratio^(-1/(order+1))`. This
is considered incorrect but maintained for compatibility.

#### Rejection

If `dt / dt_new > tran_redofactor` (default 2.5), the timepoint is rejected:
- Solution is discarded
- Timestep is cut to `dt_new`
- Solver retries from the previous accepted point

Setting `tran_redofactor=0` disables LTE-based rejection entirely.

### Timestep Control

#### Maximum Timestep (hmax)

`hmax` is the minimum of:
1. `(stop - start) / tran_minpts` — ensures minimum output density
2. `tran_fbr × breakpoint_distance` — limits dt near breakpoints
3. `maxstep` — user-specified limit from analysis line
4. Device-requested `$bound_step` — via OSDI callback

#### Minimum Timestep

VACASK has no explicit minimum timestep. The only floor is floating-point precision:

```cpp
const double timeRelativeTolerance = std::numeric_limits<double>::epsilon() * 8;
// ≈ 1.78e-15

if (hk < tSolve * timeRelativeTolerance) {
    setError(TranError::TimestepTooSmall);
    // Abort simulation
}
```

At t=100µs, this means min_dt ≈ 1.78e-19 seconds — effectively no fixed limit.

#### Timestep on NR Failure

When NR fails to converge within `tran_itl` iterations:
1. Timepoint is rejected
2. Integration order drops to 1
3. Timestep is cut by `tran_ft` (default 0.25, i.e., quartered)
4. `pointsSinceLastDiscontinuity` resets

#### Step Budget

VACASK has no limit on total number of transient steps. The simulation runs until
`t >= stop` or a fatal error occurs. Long simulations (e.g., Graetz at 1s) may
require millions of steps.

### Discontinuity Handling

When a breakpoint or `$discontinuity` event is encountered:
1. Integration order drops to 1
2. `pointsSinceLastDiscontinuity` resets to 1
3. LTE warmup period restarts (order stays at 1 until enough history accumulates)

This prevents LTE over-estimation after fast transients by discarding stale history.

---

## Small-Signal Analysis

All small-signal analyses share:

| Option | Default | Description |
|--------|---------|-------------|
| `smsig_debug` | 0 | Debug: 1=iteration/homotopy, ≥100=print linear system, ≥101=print before checks |

### AC Analysis Parameters

```
analysis ac1 ac from=<f1> to=<f2> mode=<sweep> points=<n>
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `from` | 0 | Start frequency (Hz) |
| `to` | 0 | Stop frequency (Hz) |
| `step` | 0 | Step size for linear sweep |
| `mode` | — | Sweep mode: `"lin"`, `"dec"`, `"oct"`, or empty (use step) |
| `points` | 0 | Points per decade (dec/oct) or total points (lin) |
| `values` | — | Vector of explicit frequency values |
| `writeop` | 0 | Dump OP to `<name>.op.raw` (1=yes) |
| `write` | 1 | Write results to output file |

### Noise Analysis Parameters

```
analysis noise1 noise out=<node> in=<source> from=<f1> to=<f2> mode="dec" points=<n>
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `out` | — | Output node(s) for noise measurement |
| `in` | — | Input source for input-referred noise |
| `from` | 0 | Start frequency (Hz) |
| `to` | 0 | Stop frequency (Hz) |
| `step` | 0 | Step size for linear sweep |
| `mode` | — | Sweep mode: `"lin"`, `"dec"`, `"oct"` |
| `points` | 0 | Points per decade (dec/oct) or total points (lin) |
| `values` | — | Vector of explicit frequency values |
| `writeop` | 0 | Dump OP to `<name>.op.raw` (1=yes) |
| `write` | 1 | Write results to output file |

### DC Transfer Function (DCxf) Parameters

```
analysis dcxf1 dcxf out=<node>
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `out` | — | Output node(s) for transfer function |
| `writeop` | 0 | Dump OP to `<name>.op.raw` (1=yes) |
| `write` | 1 | Write results to output file |

### DC Incremental (DCinc) Parameters

```
analysis dcinc1 dcinc
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `writeop` | 0 | Dump OP to `<name>.op.raw` (1=yes) |
| `write` | 1 | Write results to output file |

### AC Transfer Function (ACxf) Parameters

```
analysis acxf1 acxf out=<node> from=<f1> to=<f2> mode=<sweep> points=<n>
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `out` | — | Output node(s) for transfer function |
| `from` | 0 | Start frequency (Hz) |
| `to` | 0 | Stop frequency (Hz) |
| `step` | 0 | Step size for linear sweep |
| `mode` | — | Sweep mode: `"lin"`, `"dec"`, `"oct"` |
| `points` | 0 | Points per decade (dec/oct) or total points (lin) |
| `values` | — | Vector of explicit frequency values |
| `writeop` | 0 | Dump OP to `<name>.op.raw` (1=yes) |
| `write` | 1 | Write results to output file |

---

## Harmonic Balance Analysis

### HB Options

| Option | Default | Description |
|--------|---------|-------------|
| `hb_itl` | 100 | Max HB iterations in non-continuation mode |
| `hb_itlcont` | 50 | Max HB iterations in continuation mode |
| `hb_skipinitial` | 1 | Skip initial HB, go straight to homotopy (1=yes) |
| `hb_homotopy` | `["src"]` | Homotopy algorithms for HB analysis |
| `hb_debug` | 0 | Debug: ≥1=iteration/homotopy, ≥2=continuation, ≥3=spectrum construction |

### HB Analysis Parameters

```
analysis hb1 hb freq=[<f1>, <f2>] nharm=<n>
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `freq` | — | Fundamental frequencies (vector) |
| `nharm` | 4 | Harmonics per frequency (scalar or vector) |
| `immax` | 0 | Max intermodulation product order (≤0 defaults to max nharm) |
| `truncate` | `"diamond"` | Truncation scheme: `"raw"`, `"box"`, `"diamond"` |
| `samplefac` | 2 | Sampling factor in time domain (≥1) |
| `nper` | 3 | Number of periods for colocation points |
| `sample` | `"random"` | Sampling mode: `"uniform"`, `"random"` |
| `harmonic` | — | Annotation: which frequencies are harmonics (raw scheme only) |
| `imorder` | — | Annotation: intermodulation product order (raw scheme only) |
| `store` | — | Slot name to store solution |
| `nodeset` | — | Slot name to read nodesets |
| `write` | 1 | Write results to output file |

---

## Device-Level Limiting

Device-level limiting is the primary convergence mechanism in VACASK. Unlike
system-level damping (`nr_damping`), limiting is applied per-device through the
`$limit` callback in Verilog-A models.

### Available Limiting Functions

VACASK provides five built-in limiting functions via OSDI callbacks
(`osdicallback.cpp`, `limitfunctions.cpp`):

| Function | OSDI callback | Init value | Purpose |
|----------|--------------|------------|---------|
| `DEVpnjlim` | `osdiPnjlim` | `vcrit` | PN junction voltage limiting |
| `DEVpnjlim` (typed) | `osdiTypedpnjlim` | `vcrit` | PN junction with NMOS/PMOS type sign |
| `DEVfetlim` | `osdiFetlim` | `vto + 0.1` | FET gate-source voltage limiting |
| `DEVlimvds` | `osdiLimvds` | `0.1` | FET drain-source voltage limiting |
| `DEVlimitlog` | `osdiLimitlog` | `0.0` | Logarithmic damping (e.g., thermal) |

All limiting functions have an **initialization mode**: on the first NR iteration
of an analysis, the `init` flag is true and the function returns a safe starting
value instead of applying limiting. This gives NR a reasonable initial operating
point.

### Vt and Vcrit

`Vt` (thermal voltage) and `Vcrit` (critical voltage) are not global simulator
options. They are **computed per-device** from the device's physical parameters and
operating temperature, then passed to the limiting functions.

#### Vt — Thermal Voltage

```
Vt = kT/q ≈ 0.02585V at 300K (27°C)
```

For PN junctions, the effective thermal voltage includes the emission coefficient:

```
Vte = n × Vt
```

Where `n` is the ideality/emission factor of the junction (typically 1.0-2.0).
A diode with `n=1.8` at 300K has `Vte = 1.8 × 0.02585 = 0.04653V`.

`Vt` determines the sensitivity of limiting:
- DEVpnjlim triggers when `|vnew - vold| > 2 × Vt`
- Logarithmic compression uses `Vt` as the step-size unit
- Smaller `Vt` means tighter limiting (more compression per volt of change)

#### Vcrit — Critical Voltage

```
Vcrit = Vte × ln(Vte / (√2 × Is))
```

Where `Is` is the junction saturation current. This is the voltage above which
the exponential junction current becomes numerically dangerous.

Example values:
- 1N4007 (Is=1e-14 A, n=1.0): `Vcrit ≈ 0.62V`
- Low-Is diode (Is=1e-18 A, n=1.0): `Vcrit ≈ 0.86V`
- High-n junction (Is=1e-14 A, n=2.0): `Vcrit ≈ 1.28V`

`Vcrit` is the threshold for enabling limiting. When `vnew > Vcrit`, the
exponential current `Is × exp(V/Vte)` is large enough that uncontrolled Newton
steps can cause overflow or massive residual changes.

### SPICE Diode Model (`sp_diode`)

The SPICE diode uses `$limit` with `DEVpnjlim` for PN junction voltage limiting.
Three mechanisms are active:

#### 1. DEVpnjlim (forward bias)

Applied when `vnew > vcrit` and `|vnew - vold| > 2×Vt` (`limitfunctions.cpp:66-98`):

```cpp
if ((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt))) {
    if (vold > 0) {
        arg = (vnew - vold) / vt;
        if (arg > 0)
            vnew = vold + vt * log(1 + arg);   // Logarithmic compression
        else
            vnew = vold - vt * log(1 - arg);
    } else {
        vnew = vt * log(vnew / vt);            // Reverse recovery
    }
} else if (vnew < 0) {
    // Negative voltage clamping (ngspice extension)
    if (vold > 0)
        arg = -1 * vold - 1;
    else
        arg = 2 * vold - 1;
    vnew = max(vnew, arg);
}
```

Note: This is an improved version of the SPICE3 algorithm by A. Buermen, fixing
the `log(arg-2)` issue in ngspice where `|vnew-vold| < 3×Vt` produces negative
log arguments. VACASK uses `log(1+arg)` instead.

After DEVpnjlim, a **last-resort voltage limiter** is applied (`voltagelimit`):
soft tanh limiting at a configurable voltage ceiling. This catches cases where
even the logarithmic compression produces excessive voltages.

#### 2. Reverse breakdown limiting

When `bv` is given and diode is deep in reverse (`vd < -bv + 10×Vt_brk`):
- Shifts voltage relative to breakdown: `vdtemp = -(vd + bv)`
- Applies `DEVpnjlim` to the shifted voltage with `vtebrk` parameters

#### 3. Initialize limiting

At the first NR iteration of DC analysis, forces `vd = Vcrit` to give NR a
good starting point.

#### $limit Callback Mechanism

The limiting works through two `$limit` calls:
1. `$limit(V, DEVlimitOldGet)` — retrieves voltage from previous NR iteration
2. `$limit(V, DEVlimitNewSet, limited_vd, flag)` — stores limited value, calls
   `$discontinuity(-1)` if limiting was applied

### DEVfetlim — FET Gate-Source Limiting

`DEVfetlim` (`limitfunctions.cpp:157-216`) limits the per-iteration change of FET
Vgs based on operating region. The algorithm is region-aware:

- **On-region** (Vgs ≥ Vto, well above threshold):
  - Max step up: `2 × |Vgs_old - Vto| + 2`
  - Max step down: `|Vgs_old - Vto| + 1`
- **Middle region** (near threshold):
  - Decreasing: floor at `Vto - 0.5`
  - Increasing: ceiling at `Vto + 4`
- **Off-region** (Vgs < Vto):
  - Max step up/down: `2 × |Vgs_old - Vto| + 2`
  - Ceiling at `Vto + 0.5` when turning on

The key insight is that larger steps are allowed when the device is firmly in the
on-region (where the model is less sensitive), but tight limiting is applied near
the threshold transition.

On initialization, returns `Vto + 0.1` to start near the threshold.

### DEVlimvds — FET Drain-Source Limiting

`DEVlimvds` (`limitfunctions.cpp:37-55`) limits Vds changes:

- **High Vds** (Vds_old ≥ 3.5V): increasing capped at `3 × Vds_old + 2`,
  decreasing floored at 2V
- **Low Vds** (Vds_old < 3.5V): increasing capped at 4V, decreasing floored at
  -0.5V

On initialization, returns `0.1V`.

### DEVlimitlog — Logarithmic Damping

`DEVlimitlog` (`limitfunctions.cpp:223-244`) provides general-purpose logarithmic
damping for any quantity (used by electrothermal models for temperature limiting):

```cpp
if (deltemp > deltemp_old + LIM_TOL)
    deltemp = deltemp_old + LIM_TOL + log10((deltemp - deltemp_old) / LIM_TOL);
else if (deltemp < deltemp_old - LIM_TOL)
    deltemp = deltemp_old - LIM_TOL - log10((deltemp_old - deltemp) / LIM_TOL);
```

On initialization, returns `0.0`.

### PSP103 Model

PSP103 does **not** use `$limit`. Instead it relies on internal numerical conditioning:

- **`expl` macro**: Safe exponential with 3rd-order polynomial extrapolation for
  `|x| > 230.26`, preventing overflow
- **`CLIP_LOW`/`CLIP_HIGH`/`CLIP_BOTH`**: Parameter clamping macros

For circuits using PSP103 (e.g., the c6288 multiplier), convergence relies
entirely on the intrinsic smoothness of the model. No external pnjlim/fetlim
is applied.

---

## Benchmark Examples

### Graetz (Diode Bridge Rectifier)

Full-wave diode bridge: 4× 1N4007 diodes, 100µF smoothing cap, 1kΩ load,
20V 50Hz sine input.

#### Configuration

```
options tran_method="gear2" nr_bypass=0 nr_contbypass=1
analysis tran1 tran step=1u stop=1 maxstep=1u
```

Uses `spice/sn/diode.osdi` — the "simplified noise" variant of `sp_diode`,
identical for DC/transient behavior.

#### MNA Structure

9 unknowns:
- 4 external node voltages (inp, inn, outp, outn)
- 4 internal diode nodes (a_int for each d1-d4, from series resistance rs=42mΩ)
- 1 voltage source branch current

#### Performance

| Metric | Value |
|--------|-------|
| Accepted timepoints | 1,000,004 |
| Rejected timepoints | **0** |
| NR iterations | 2,000,275 |
| Average NR/step | **2.0** |
| Wall time | 1.3s |
| Smallest dt | ~2.2ns |

#### Active Limiting

| Mechanism | Active? |
|-----------|---------|
| DEVpnjlim (forward) | Yes, 4 diodes |
| Breakdown limiting (bv=1kV) | Yes, deep reverse only |
| initialize_limiting | Yes, first OP iteration |
| Global nr_damping | 1.0 (disabled) |

#### LTE Behaviour

Three distinct phases observed:

**Startup (steps 2-7)**: dt grows from 250ns → 6µs as LTE ratio is very small
(0.003-0.01).

**Steady shrinking (steps 15-24)**: LTE ratio increases from 0.16 → 1.45 as the
sine wave ramps up. Timestep gradually shrinks from 6µs → 3.35µs. All points
accepted — ratio never exceeds `tran_redofactor` (2.5).

**Equilibrium (steps 35+)**: LTE ratio hovers at ~1.0, timestep stabilizes at
~3.35µs. Self-regulating: ratio slightly > 1 → dt shrinks fractionally; ratio
slightly < 1 → dt grows fractionally.

### c6288 (CMOS Multiplier)

16×16 CMOS multiplier with PSP103 transistors.

#### Configuration

```
options nr_convtol=1 nr_bypasstol=1 nr_bypass=0 nr_contbypass=1 tran_lteratio=1.5
analysis tran stop=2n step=2p icmode="uic"
```

Note: `nr_convtol=1` and `nr_bypasstol=1` relax per-instance convergence checks
(matching full tolerances rather than the default 1% of tolerance).

#### Performance

| Metric | Value |
|--------|-------|
| Accepted timepoints | 1,024 |
| Rejected timepoints | 8 |
| Total NR iterations | ~3,482 |
| Average NR/step | **3.4** |
| Wall time | 58.8s |
| Smallest dt | ~2ps |

#### Active Limiting

| Mechanism | Active? |
|-----------|---------|
| DEVpnjlim | **Not used** (PSP103 has no `$limit`) |
| PSP103 `expl` macro | Yes (safe exponential) |
| PSP103 CLIP macros | Yes (parameter clamping) |
| Global nr_damping | 1.0 (disabled) |

#### hmax

`hmax = stop / tran_minpts = 2ns / 50 = 40ps`. During fast switching, LTE is the
dominant timestep constraint.
