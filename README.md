![cask](./cask.png) 

# About VACASK
VACASK (Verilog-A Circuit Analysis Kernel) is an analog circuit simulator. VACASK uses the [OpenVAF-reloaded Verilog-A compiler](https://github.com/arpadbuermen/OpenVAF) for building the device models as shared libraries. The compiled device models are loaded by the simulator on demand at runtime. The simulator communicates with the models via the [OSDI API](https://openvaf.semimod.de/docs/details/osdi/). Currently OSDI API 0.4 is used which is supported only by OpenVAF-reloaded. Of course, you can also create device models using VACASK's APIs in C++ and link them statically with the simulator. 

VACASK is not SPICE (although one could write a SPICE-compatible netlist parser for it with little effort). SPICE3 is more than 30 years old, written in C, and the code is hard to maintain. In some respect SPICE looks more like a proof of concept one writes before building the real thing. The way circuit equations are handled in SPICE makes it hard to extend the simulator with new algorithms. VACASK's goal is to be better than SPICE, not only in terms of what it offers, but also in terms of extensibility and ease of maintenance. 

# Wait... is it any good? 

The following benchmark results were obtained on the [C6288 16x16 multiplier circuit](benchmark/c6288) multiplying 0xFFFF with itself (simulated as an analog circuit). This is a medium size circuit with 10112 transistors and 25380 nodes. The transistors were modelled with the PSP103.4 model. Ngspice and Xyce tolerances were set so that the number of output timepoints is roughly equal to that of VACASK. OpenMP support was disabled to make the comparison fair. All simulators used KLU as the linear solver. Results output was disabled so that the impact of disk operations was minimized. The benchmark was run on a computer with an AMD Threadripper 7970 processor. 

|Simulator  |Time (s)        |Timepoints |Rejected |Iterations |
|-----------|----------------|-----------|---------|-----------|
|Xyce       | 151.57         |1013       |37       |3559       |
|Ngspice    |  71.81         |1020       |1        |3474       |
|VACASK     |  57.98         |1021       |7        |3487       |

If you want to find out more, there is [a page dedicated to benchmarks](benchmark).

# What does VACASK offer? 

- user defined global and ground nodes
- fully parameterized hierarchical circuit description
- RPN interpreter for parameterized expression evaluation
- integer, real, and string data types
- vectors (homogeneous containers) and lists (heterogeneous containers)
- a library of [built-in functions and constants](lib/context.cpp) for use in parameterized expressions
- [operating point](include/coreop.h), [DC small-signal](include/coredcinc.h), [DC transfer function](include/coredcxf.h), [AC small-signal](include/coreac.h), [AC transfer function](include/coreacxf.h), [noise](include/corenoise.h), [transient](include/coretran.h), and [(multitone) harmonic balance](include/corehb.h) analyses 
- [options](lib/options.cpp) for fine tuning the simulator
- selection of what should be saved during simulation (save directives)
- collection of auxiliary values (opvars) computed by device models
- [parametric sweep](include/answeep.h) of any analysis with arbitrary depth
- almost anything can be swept (instance, model, and subcircuit parameters, options, and circuit variables)
- anything that can be swept can also be modified without reloading the circuit (no need to build a new netlist and restart the simulator)
- automatic partial circuit elaboration when circuit topology changes due to a change in parameters, options, or variables
- a different topology can be elaborated without restarting the simulator (circuits with multiple testbenches can be simulated without restarting the simulator)
- residual-based convergence test (improved accuracy of the nonlinear solver)
- inactive instances bypass in nonlinear solver (disabled by default, see the `nr_bypass`, `nr_convtol`, and `nr_bypasstol` options)
- instance evaluation bypass in the first iteration of nonlinear solver running in continuation mode (enabled by default, see the `nr_contbypass` option)
- several homotopy algorithms (SPICE3/adaptive gmin stepping and source stepping) for finding the operating point of problematic circuits
- nodesets for improving convergence speed and selecting the operating point
- analysis results can be used as nodesets for subsequent analyses (combined with circuit variable sweeps this feature can be used for implementing custom arc length homotopy algorithms at the netlist level)
- setting initial conditions (Spectre style and legacy SPICE3 style)
- backward Euler, trapezoidal, and Gear integration algorithms
- predictor-corrector local truncation error control in transient analysis
- numerical solvers based on the KLU sparse matrix library
- SPICE ASCII/binary raw file output
- embedded postprocessing scripts and Verilog-A models in the netlist
- postprocessing of results with external tools (some basic Python/Numpy scripts are provided)
- a growing [library of Verilog-A models](devices/) (diode, BSIM3, BSIM4, BSIMBULK, ...)
- simulator library that can be linked to 3rd party software
- circuit can be built programmatically or read by a custom parser
- netlist parser with Spectre-like syntax
 
Certain devices (independent voltage and current sources, linear controlled sources, and inductive coupling) are implemented as builtin devices because certain features needed by these devices are not available in OpenVAF-reloaded or even Verilog-A. 

VACASK is being developed by Árpád Bűrmen at the EDA Laboratory, University of Ljubljana, Slovenia. It is written in C++20 and is free software released under the [GNU Affero General Public License 3.0](LICENSE). 

# What about device models? 

The following device models are supplied with VACASK. 

|Builtin device                   |Name   |
|---------------------------------|-------|
|Independent voltage source       |vsource|
|Independent current source       |isource|
|Voltage-controlled voltage source|vcvs   |
|Voltage-controlled current source|vccs   |
|Current-controlled voltage source|ccvs   |
|Current-controlled current source|cccs   |
|Inductive coupling               |mutual |

|Verilog-A device          |File               |Module   |
|--------------------------|-------------------|---------|
|Linear resistor           |resistor.va        |resistor |
|Linear capacitor          |capacitor.va       |capacitor|
|Linear inductor           |inductor.va        |inductor |
|SPICE diode               |diode.va           |diode    |
|BSIM3v3 MOSFET (Cogenda)  |bsim3v3.va         |bsim3    |
|BSIM4v8 MOSFET (Cogenda)  |bsim4v8.va         |bsim4    |
|PSP103.4 MOSFET           |psp103v4/psp103.va |PSP103VA |
|BSIMBULK MOSFET 106.2.0   |bsimbulk.va        |bsimbulk |

All Verilog-A models supplied with VACASK are located in [devices](devices). You can find several models at [www.mos-ak.org](https://www.mos-ak.org/open_dir/). All recent models developed by the [BSIM group at UC Berkeley](https://bsim.berkeley.edu/) are released in Verilog-A. Also take a look at [The Designer's Guide community](https://designers-guide.org/index.html) where various models are available in the [Verilog AMS section](https://designers-guide.org/verilog-ams/index.html). 

Recently a project ([Verilog-A Distiller](https://codeberg.org/arpadbuermen/VADistiller)) whose aim is to create a converter from SPICE3 C model format to Verilog-A took off. At this point the following converted models are available in VACASK. 

|Verilog-A device (SPICE)  |File        |Module       |
|--------------------------|------------|-------------|
|Linear resistor           |resistor.va |sp_resistor  |
|Linear capacitor          |capacitor.va|sp_capacitor |
|Linear inductor           |inductor.va |sp_inductor  |
|Diode (levels 1 and 3)    |diode.va    |sp_diode     |
|Gummel-Poon BJT           |bjt.va      |sp_bjt       |
|JFET level 1 (Schichman-Hodges)     |jfet1.va    |sp_jfet1     |
|JFET level 2 (Parker-Skellern) *    |jfet2.va    |sp_jfet2     |
|MESFET level 1 (Statz et. al.) *    |mes1.va     |sp_mes1      |
|MOSFET level 1 (Schichman-Hodges) * |mos1.va     |sp_mos1      |
|MOSFET level 2 (Grove-Frohman) *    |mos2.va     |sp_mos2      |
|MOSFET level 3 (empirical) *        |mos3.va     |sp_mos3      |
|MOSFET level 6 (Sakurai-Newton) *   |mos6.va     |sp_mos6      |
|MOSFET level 9 (modified level 3) * |mos9.va     |sp_mos9      |
|BSIM3 3.3.0                         |bsim3v30.va |sp_bsim3v30  |

Devices marked with an asterisk (*) do not conserve charge because of the modeling approach chosen by their respective authors. 

The converted SPICE models can be found in the [devices/spice](devices/spice) directory.

Most devices provide several model variants. The `sn` variant ([devices/spice/sn](devices/spice/sn) directory) models have a simplified noise model and do not expose opvars that would introduce extra internal nodes. These models are the fastest, but they cannot be used in advanced noise analyses (they give correct noise values in the ordinary small-signal noise analysis only). 

The `full` variant  ([devices/spice/full](devices/spice/full) directory) of a model exposes all opvars and the noise model is appropriate for all types of analysis. 

The `default` variant of models can be found in the [devices/spice](devices/spice) directory. These models do not expose opvars that introduce extra internal nodes. The noise model, however, is appropriate for all types of noise analysis. If a particular device does not have a `sn` or a `full` variant then that variant is equal to the `default` variant. For more information consult the [Verilog-A Distiller repository](https://codeberg.org/arpadbuermen/VADistiller). 

Examples of SPICE3 model usage are in [demo/spice](demo/spice). 


# Installation from pre-built packages
[Pre-built packages](https://codeberg.org/arpadbuermen/VACASK/releases) for Linux (based on the stable version of Debian) and Windows are available. The OpenVAF-reloaded compiler is included in all binary packages. Linux users can choose between a .tgz archive and a .deb package. The Windows package is a .zip file that you can unpack wherever you want. It is recommended to add the `bin` directory to the system path. 

A new version of VACASK is released every now and then. Between releases [(not quite) nightly builds](https://fides.fe.uni-lj.si/vacask/download/) are released. These are great if you want to try VACASK with latest bugfixes. 

# Getting started
There are some examples available in the [`demo`](demo) directory. You can try the simulation of a Miller OTA by running
```
vacask demo/bsim3-ptm-amp/toplevel.sim
```

If you have Python 3, NumPy, and [Matplotlib](https://matplotlib.org/) installed the results will be plotted by the postprocessor script. 

You can learn about the netlist syntax by studying the demos in the [`demo`](demo) directory and the system tests in the [`test`](test) directory. Documentation is planned for the future. :)

If VACASK fails to find something, first check all the paths by typing
```
vacask -dp
```

If you specify the `-df` option VACASK will print the paths to the files it is loading, dumping, or compiling. 

VACASK detects the Python 3 interpreter and sets the `PYTHON` circuit variable to the interpreter's full path. This variable can then be used in the netlist for launching Python to postprocess the simulation results without having to specify its full path. For VACASK to find the Python interpreter the interpreter's binary directory must be in the system path. VACASK supplements the `PYTHONPATH` variable with the directory holding the supplied Python scripts (`<vacask library directory>/python`). These scripts can be used for loading binary .raw files. They depend on the [NumPy library](https://numpy.org/). 

When a file is included with the `include` netlist directive and the given path is absolute VACASK loads it based on the given absolute path. If the path is relative VACASK first looks for the file in the directory where the netlist that invoked the `include` directive resides, then in the current working directory, and finally in the include files path. The include files path is by default set to`<vacask library directory>/inc`. You can override it by setting the `SIM_INCLUDE_PATH` environmental variable. The directories in the list must be separated by colons (in Windows they must be separated by semicolons). 

Models are loaded with the `load` netlist directive. If the given path is absolute VACASK looks for the model only at the given path. If, however, it is relative VACASK first searches for the model in the directory where the netlist invoking the `load` directive is located, followed by the current working directory, and the modules search path. The modules search path is by default set to `<vacask library directory>/mod`. You can override it with the `SIM_MODULE_PATH` environmental variable (same syntax as for `SIM_INCLUDE_PATH`). 

VACASK can compile Verilog-A files on the fly. For that purpose VACASK looks for the OpenVAF-reloaded compiler in the directory where the VACASK binary is installed and in the system path. You can override this by specifying the path to the OpenVAF-reloaded compiler in the `SIM_OPENVAF` envirnonmental variable. If a `load` directive specifies a raw Verilog-A file (ending in .va), VACASK will try to compile it. The compiled model is placed in the current working directory and then loaded. 

# Building VACASK
VACASK has only a few dependencies. You will need a C++20 compiler with an implementation of the standard C++ library, the Boost library (use version 1.84, versions up to and including 1.87 reportedly work, 1.88 does not), and the KLU library (SuiteSparse). All these components come as pre-built packages for [Debian](https://www.debian.org) (and other Linux distributions). You will also need a working Python3 installation (for the system tests). 

You will also have to install the OpenVAF-reloaded compiler. The latest development version of OpenVAF-reloaded can be found at [https://fides.fe.uni-lj.si/openvaf/download](https://fides.fe.uni-lj.si/openvaf/download/). Make sure you download the OSDI 0.4 version. Of course, you can also take the OpenVAF-reloaded binary from the VACASK binary packages (.deb and .tar.gz for Linux, .zip for Windows). Make sure the OpenVAF binary you pick up is named `openvaf-r` (this indicates the compiler produces models with the OSDI 0.4 interface). If you decide to build the compiler yourself, git-clone the [OpenVAF-reloaded repository](https://github.com/arpadbuermen/OpenVAF). Instructions for building can be found in the [README.md](https://github.com/arpadbuermen/OpenVAF/blob/master/README.md) file. 

## Linux
Install gcc, Boost, and KLU. You will also need CMake and GNU make or Ninja for building. 

First, create a `build` directory and create the build system
```
cmake -G Ninja  -S <sources directory> -B <build directory> -DCMAKE_BUILD_TYPE=Release -DOPENVAF_DIR=<path to the OpenVAF-reloaded compiler>
```

To build with GNU make, replace `-G Ninja` with `-G "Unix Makefiles"`. The build process is started by typing
```
cmake --build <build directory>
```

After the build process is finished the binary can be found in `<build directory>/simulator`. 

The .tar.gz and .deb packages can be built by changing the current directory to `<build directory>` and typing 
```
cpack
```
The packages are created in the `<build directory>`. 

## Windows
Building for Windows is performed with the Mingw64 compiler. Unfortunately you will have to build all of the prerequisites manually. 

### Building the prerequisites
First, install the compiler and the tools. We will assume everything will be unpacked and built in `e:\`. Adjust the paths accordingly if you are going to use a different drive. Download Msys2 from [https://www.msys2.org/](https://www.msys2.org/) and install it to `e:\msys64`. Start the MSYS prompt and install bison, flex, and binutils. 
```
pacman -S bison
pacman -S flex
pacman -S binutils
```

Download MinGW64 from [https://winlibs.com/](https://winlibs.com/). Get the posix-seh-ucrt version without LLVM/Clang/LLD/LLDB (we don't need them). Unpack it to `e:\mingw64`. This version of MinGW64 comes with CMake and Ninja (required for building). 

Add MSYS2 and MinGW64 to the system path (Control Panel/System/Advanced System Settings/...). Make sure you add MinGW (`e:\mingw64\bin`) before MSYS2 (`e:\msys64\usr\bin`). 

Start the command prompt (`cmd.exe`). Create a directory for all the stuff we are going to build.
```
e:
cd \
mkdir build
```

Next, download the Boost library from [https://www.boost.org/users/download/](https://www.boost.org/users/download/). Use Boost version 1.84.0. Get the Windows version and unpack it in `e:\build`. In the Boost sources directory type
```
cd tools\build
bootstrap mingw
cd ..\..
tools\build\b2 --with-filesystem link=static toolset=gcc
```
After a short time the required part of Boost is built and placed in the `stage` subdirectory.

Prepare the toolchain file (`e:\build\mingw.cmake`) by putting the following definitions in the file. 
```
set(CMAKE_SYSTEM_NAME Windows)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_RC_COMPILER windres)

set(CMAKE_FIND_ROOT_PATH e:/mingw64/bin)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
```
If needed, replace `e:/mingw64/bin` with your path. 

Now, download OpenBLAS from [https://www.openblas.net/](https://www.openblas.net/). There is a link to [GitHub where the released sources can be found](https://github.com/OpenMathLib/OpenBLAS/releases). Take the latest source (at the time of writing 0.3.29). Unpack it in `e:\build`. Enter the directory with the OpenBLAS sources and type
```
mkdir build
cd build
cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE=e:\build\mingw.cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC_LIBS=ON -DUSE_THREADS=0 -DUSE_LOCKING=1 -DDYNAMIC_LIST="CORE2;NEHALEM;BARCELONA;SANDYBRIDGE;BULLDOZER;PILEDRIVER;STEAMROLLER;EXCAVATOR;HASWELL;ZEN;SKYLAKEX;COOPERLAKE;SAPPHIRERAPIDS" -DTARGET=NEHALEM
cmake --build . -j 8
cmake --install . --prefix e:/build/installation
```

The `-j 8` option enables parallel building with 8 processors. Since OpenBLAS is big, this will save you some time. In the end OpenBLAS will be installed in `e:\build\installation`. 

Finally, download [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) from [GitHub](https://github.com/DrTimothyAldenDavis/SuiteSparse/releases). Get the latest release source code (at the time of writing 7.10.3). Unpack it in `e:\build`. In the sources directory type
```
mkdir build
cd build
cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE=e:\build\mingw.cmake -DSUITESPARSE_ENABLE_PROJECTS="klu" -DCMAKE_BUILD_TYPE=Release -DBLAS_LIBRARIES=e:\build\installation\lib\libopenblas.a -DLAPACK_INCLUDE_DIRS=e:\buid\installation\include\openblas -DBLAS_INCLUDE_DIRS=e:\build\installation\include\openblas
cmake --build . -j 8
cmake --install . --prefix e:/build/installation
```

Replace the `e:\...` paths with your own, if needed. In the end OpenBLAS will be installed in `e:\build\installation`. 

### Building the simulator
Unpack the sources, create a build directory, and type. 
```
cmake .. -G Ninja -S <sources directory> -B <build directory> -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=e:\build\mingw.cmake -DOPENVAF_DIR=<path to the OpenVAF-reloaded compiler> -DBoost_ROOT=e:/build/boost_1_84_0/stage -DSuiteSparse_DIR=e:/build/installation
cmake --build <build directory>
```
Replace the `e:\...` paths with your own, if needed. In the end the simulator can be found in `<build directory>/simulator`. To create a package (.zip), go to the `<build directory>` and type. 
```
cpack
```
The created packages are located in the `<build directory>`. 

# Visual Studio Code project for developers
A [Visual Studio Code](https://code.visualstudio.com/) setup is available in the [`.vscode`](.vscode) subdirectory of the sources. Files [`settings-linux.json`](.vscode/settings-linux.json) and [`settings-windows.json`](.vscode/settings-windows.json) are the settings templates for Linux and Windows. Depending on your platform copy one of these two to `settings.json` and edit it to reflect your configuration. File `settings.json` is not tracked by git so editing it won't result in any changes that need committing. 

Install the following Visual Studio Code extensions: C/C++, CMake, CMake Tools, VSCode-YACC, and VSCode-YYLEX. Assuming the Visual Studio Code binary is in the system path you can start the IDE by entering the simulator source directory and typing
```
code .
```

In Windows select the MinGW64 toolchain. In Linux select GCC. Configure the project with Ctrl+Shift+P 'CMake: Delete Cache and Reconfigure', followed by building with Ctrl+Shift+P 'CMake: Build'. A full debugging setup is available in [`launch.json`](.vscode/launch.json). System tests are located in [`test`](test) and can be run via CMake/CTest. The path to the built debug version (relative to the sources) is `../build.VACASK/Debug`. The release version is built under `../build.VACASK/Release`. 
