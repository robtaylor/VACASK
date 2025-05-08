# Simulators

Benchmarks were performed with [Xyce 7.9](https://xyce.sandia.gov/), [Gnucap 20240220](http://gnucap.org), [Ngspice pre-master 45](https://ngspice.sourceforge.io/), and [VACASK](https://codeberg.org/arpadbuermen/VACASK). 

All the tests measure single threaded performance. Xyce, Ngspice, and VACASK all used KLU as the linear solver. Gnucap used its own solver. 

Xyce was compiled without parallel computing support. Trilinos 16.0.0 was used. Since it has an obsession with taking time measurements this can result in significantly slowdowns for small circuits. To disable time measurements the file `Xyce/src/UtilityPKG/N_UTL_StatMetricTraits.C` was modified by replacing `return Xyce::cpu_time();` and `return Xyce::wall_time();` with `return 0;` (thanks to Eric Keiter for a tip). This modified version is named Xyce-fast. 

Gnucap was compiled with default options. Version 20240220 was the last stable release at the time of testing. We had convergence problems with version 20250301-dev so we stuck with the stable release. 

Ngspice was configured with 
```
configure --with-x --with-readline=yes --enable-osdi --disable-debug --disable-openmp
```
and then compiled. 

VACASK was compiled with default options. 

