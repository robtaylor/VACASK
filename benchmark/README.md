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

Ngspice, Xyce and Gnucap were using builtin models of circuit elements. VACASK used Ngspice models converted into Verilog-A by [Verilog-A Distiller](https://codeberg.org/arpadbuermen/VADistiller). The simplified noise model was used (`simplified_noise` set to `True`) and opvars that induce extra internal nodes were left out (`opvars_intnodes` set to `False`) to make models as close as possible to the ones provided by Ngspice in terms of complexity. The CMOS test problems used the PSP103.4 MOSFET model. 

For each test problem the number of timepoints and the number of Newton-Raphson iterations (residual evaluations for Xyce) is listed. In some cases the number of rejected timepoints is also given. 

# Results

## RC circuit excited by a pulse train (rc)
This is a simple circuit with 2 elements (resistor, capacitor) excited by a voltage pulse train. The timestep was forced artificially to a small value so that roughly 1000000 time steps are computed. This is a linear circuit with no rejected timepoints. 

|Simulator  |Time (s)        |Timepoints |Iterations |
|-----------|----------------|-----------|-----------|
|Xyce       | 9.39           |1011527    |2029571    |
|Xyce-fast  | 4.12           |1011527    |2029571    |
|Gnucap     | 8.54           |1006982    |3025043    |
|Ngspice    | 1.31           |1006013    |2012031    |
|VACASK     | 0.95           |1005006    |2010014    |

## Full wave rectifier with smoothing and load (graetz)
A full wave rectifier (4 diodes) with a capacitor and a resistor as load excited by a sinusoidal voltage. Two 1GOhm resistors are used for setting up a DC path to ground. The timestep was forced artificially to a small value so that roughly 1000000 time steps are computed. Second order Gear integration was used. The diode model has transit time (tt) set to 0 because Gnucap exhibits convergence problem when reverse recovery effects are modelled. 

|Simulator  |Time (s)        |Timepoints |Rejected |Iterations |
|-----------|----------------|-----------|---------|-----------|
|Xyce       |10.60           |1000002    |0        |2000014    |
|Xyce-fast  | 5.34           |1000002    |12       |2000014    |
|Gnucap     |15.16           |1000026    |739      |4459517    |
|Ngspice    | 2.21           |1000008    |0        |2000024    |
|VACASK     | 1.88           |1000003    |0        |2000277    |

## Diode voltage multiplier (mul)
A voltage multiplier (4 diodes, 4 capacitors) with a series resistor at its input excited by a sinusoidal voltage. The timestep was forced artificially to a small value so that roughly 500000 time steps are computed. Second order Gear integration was used. The diode model has transit time (tt) set to 0 because Gnucap exhibits convergence problem when reverse recovery effects are modelled. 

|Simulator  |Time (s)        |Timepoints |Rejected |Iterations |
|-----------|----------------|-----------|---------|-----------|
|Xyce       | 5.51           |502341     |1270     |1041833    |
|Xyce-fast  | 2.78           |502341     |1270     |1041833    |
|Gnucap     | 9.94           |520797     |739      |2822527    |
|Ngspice    | 1.16           |500467     |957      |1019733    |
|VACASK     | 0.96           |500056     |3        |1001217    |

## 9 stage CMOS ring oscillator (ring)
This is a ring oscillator with 9 CMOS inverters (18 transistors) powered by 1.2V. The timestep was limited to 50ps. Xyce timeint reltol parameter was set to 5e-3 to make the number of computed timepoints roughly equal to that of VACASK. 

|Simulator  |Time (s)        |Timepoints |Rejected |Iterations |
|-----------|----------------|-----------|---------|-----------|
|Xyce       | 3.33           |27310      |0        |95462      |
|Xyce-fast  | 3.10           |27310      |0        |95462      |
|Ngspice    | 1.56           |20483      |890      |79407      |
|VACASK     | 1.19           |26066      |0        |81875      |

## 16x16 CMOS multiplier (c6288)
A medium size digital circuit with 10112 transistors and 25380 nodes. It computes the product of 0xFFFF with itself. The timestep was not limited. The simulation is purely analog (i.e. transistor-level). Due to its size we do not test it with Xyce-fast since the timing overheads are small compared to simulation time. VACASK tran_lteratio parameter was set to 1.5 to roughly match the number of computed timepoints to those of Ngspice. For the same reason Xyce timeint reltol parameter was set 2.5e-4. 

|Simulator  |Time (s)        |Timepoints |Rejected |Iterations |
|-----------|----------------|-----------|---------|-----------|
|Xyce       | 151.57         |1013       |37       |3559       |
|Ngspice    |  73.16         |1020       |1        |3474       |
|VACASK     |  61.52         |1021       |7        |3487       |



### Effect of continuation bypass 

### Effect of residual tolerance checks

