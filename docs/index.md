# VACASK Documentation

**VACASK** (Verilog-A Circuit Analysis Kernel) is an analog circuit simulator written in C++20. It uses the OpenVAF-reloaded Verilog-A compiler to build device models as shared libraries, which are loaded at runtime via the OSDI 0.4 API.

## Features

- Multiple analysis types: Operating Point, Transient, AC, Noise, Harmonic Balance
- Verilog-A device models via OpenVAF-reloaded compiler
- OSDI 0.4 API for device model interface
- KLU sparse matrix solver
- Variable-order integration methods (Adams-Moulton, BDF/Gear)
- Adaptive timestep control with LTE estimation
- Homotopy methods for convergence assistance

## Quick Start

```bash
# Run a simulation
vacask netlist.sim

# Print search paths
vacask -dp

# Print file paths being loaded
vacask -df netlist.sim
```

## Analysis Types

VACASK supports the following analysis types:

| Analysis | Command | Description |
|----------|---------|-------------|
| **op** | `op` | Operating Point (DC bias point) |
| **tran** | `tran` | Transient (time-domain) |
| **ac** | `ac` | AC small-signal frequency response |
| **noise** | `noise` | Small-signal noise analysis |
| **hb** | `hb` | Harmonic Balance (steady-state periodic) |
| **dcinc** | `dcinc` | DC incremental (small-signal at DC) |
| **dcxf** | `dcxf` | DC transfer function |
| **acxf** | `acxf` | AC transfer function |

See [Analysis Types](analysis_types.md) for detailed information.

## Documentation

- [Analysis Types Overview](analysis_types.md) - Overview of all supported analyses
- [Operating Point Analysis](operating_point.md) - Detailed OP algorithm description
- [Transient Analysis](transient.md) - Detailed transient algorithm description

## Building

### Prerequisites

- C++20 compiler (GCC/Clang)
- CMake 3.18+
- Boost 1.88 (filesystem, process, system components)
- SuiteSparse (KLU library)
- toml++ library (version 3.4)
- Bison and Flex
- OpenVAF-reloaded compiler (openvaf-r)

### Build Commands

```bash
# Configure with Ninja
cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DOPENVAF_DIR=<path-to-openvaf-r> \
    -DBoost_ROOT=<boost-directory>/stage

# Build
cmake --build build

# Run tests
cd build && ctest
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `SIM_INCLUDE_PATH` | Override include files search path |
| `SIM_MODULE_PATH` | Override compiled modules search path |
| `SIM_OPENVAF` | Override OpenVAF-reloaded compiler path |

## Configuration

TOML configuration files are read in order:

1. `/etc/vacask/vacaskrc.toml` (Linux) or `<install>/lib/vacaskrc.toml` (Windows)
2. `~/.vacaskrc.toml`
3. `.vacaskrc.toml` in startup directory
4. `.vacaskrc.toml` in netlist directory

## License

See the repository for license information.
