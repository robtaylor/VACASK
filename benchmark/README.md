# Simulators

Benchmarks were performed with [Xyce 7.9](https://xyce.sandia.gov/), [Gnucap 20240220](http://gnucap.org), [Ngspice pre-master 45](https://ngspice.sourceforge.io/), and [VACASK](https://codeberg.org/arpadbuermen/VACASK). 

All the tests measure single threaded performance. Xyce, Ngspice, and VACASK all used KLU as the linear solver. Gnucap used its own solver. 

Xyce was compiled without parallel computing support. Trilinos 16.0.0 was used. Since it has an obsession with taking time measurements this can result in significantly slowdowns for small circuits. To disable time measurements the file `Xyce/src/UtilityPKG/N_UTL_StatMetricTraits.C` was modified by replacing `return Xyce::cpu_time();` and `return Xyce::wall_time();` with `return 0;` (thank goes to Eric Keiter for suggesting this). This modified version is named Xyce-fast. 

Gnucap was compiled with default options. Version 20240220 was the last stable release at the time of testing. We had convergence problems with version 20250301-dev so we stuck with the stable release. 

Ngspice was configured with 
```
configure --with-x --with-readline=yes --enable-osdi --disable-debug --disable-openmp
```
and then compiled. 

VACASK was compiled with default options. 

# Methodology

Each case was run 6 times. The first run was not timed. We assumed that the caches are filled after the first run so the performance variability in subsequent runs is smaller. No results were stored in files to prevent disk IO operations from affecting the results. Timing includes the whole run - loading of the input files, parsing, elaboration, and simulation. The tables list the average values across 5 runs. The standard deviation of timings was below 2% of the average. You can review the methodology by looking at [benchmark.py](benchmark.py). 

Ngspice, Xyce and Gnucap were using builtin models of circuit elements. VACASK used Ngspice models converted into Verilog-A by [Verilog-A Distiller](https://codeberg.org/arpadbuermen/VADistiller). The CMOS test problems used the PSP103.4 MOSFET model. 

For each test problem the number of timepoints and the number of Newton-Raphson iterations (linear solves) is listed. In some cases the number of rejected timepoints is also given. 

# Results

## RC circuit excited by a pulse train (rc)
This is a simple circuit with 3 elements (voltage source, resistor, capacitor). The timestep was forced artificially to a small value so that roughly 1000000 time steps are computed. 

|Simulator  |Time (s)        |Timepoints       |Linear solves|
|-----------|----------------|-----------------|-------------|
|Xyce       |9.39            |1011527          |1016041
|Xyce-fast  |4.12            |1011527          |1016041      |
|Gnucap     |8.54            |1006982          |3025043      |
|Ngspice    |1.31            |1006013          |2012031      |
|VACASK     |0.95            |1005006          |2010014      |

## Full wave rectifier with smoothing and load (graetz)

## Diode voltage multiplier (mul)

## 9 stage CMOS ring oscillator (ring)

## 16x16 CMOS multiplier (c6288)

