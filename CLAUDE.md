# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the Python component of VACASK (Verilog-A Circuit Analysis Kernel), specifically focused on multi-dialect SPICE netlist parsing and conversion to VACASK format. The Python code lives in the `python/` directory.

## Build and Test Commands

```bash
# All commands run from python/ directory
cd python

# Install dependencies (uses uv)
uv sync --extra dev

# Run all tests
uv run pytest -v

# Run a single test
uv run pytest tests/test_parser.py::test_function_name -v

# Lint (new parser code only)
uv run ruff check spiceparser/ vcwriter/ tests/

# Lint with fixes
uv run ruff check --fix spiceparser/ vcwriter/ tests/
```

## Architecture

### Multi-Dialect SPICE Parser (`spiceparser/`)

The parser supports three SPICE dialects: **Ngspice**, **HSPICE**, and **LTSpice**.

**Core Components:**
- `dialect.py` - Abstract `SpiceDialect` base class with `@register_dialect` decorator for auto-registration
- `parser.py` - `NetlistParser` class that uses dialect instances to handle syntax differences
- `netlist.py` - Dialect-independent data structures: `Netlist`, `ModelDef`, `Subcircuit`, `Instance`
- `elements.py` - `DEVICE_TYPES` and `OSDI_MODULES` registries for device type info and OSDI module mappings

**Dialect Implementations** (`dialects/`):
- `ngspice.py` - Reference implementation based on original netlist_converter
- `hspice.py` - HSPICE with `.if/.endif` conditionals, `.ALTER`, single-quoted expressions
- `ltspice.py` - LTSpice with `Rser`/`Lser`/`Rpar`/`Cpar` parasitic parameters

**Key Pattern:** Dialects are auto-registered via decorator and retrieved by name:
```python
from spiceparser import parse_netlist, detect_dialect_from_file

dialect = detect_dialect_from_file("circuit.sp")  # Returns "ngspice", "hspice", or "ltspice"
netlist = parse_netlist("circuit.sp", dialect=dialect)
```

### VACASK Writer (`vcwriter/`)

- `writer.py` - `VacaskWriter` class converts parsed netlists to VACASK format (Spectre-like syntax)
- Generates `load` directives for OSDI modules and converts model/subcircuit definitions

### Legacy Code (`netlist_converter/`)

Original Ngspice-specific converter being refactored into the multi-dialect architecture. The new `spiceparser` and `vcwriter` packages are the modern replacement.

### Converter Scripts

- `netlist_converter.py` - Ngspice to VACASK converter (uses legacy netlist_converter)
- `hspice2vc.py` - HSPICE to VACASK converter
- `ltspice2vc.py` - LTSpice to VACASK converter
- `sg13g2tovc.py` - IHP SG13G2 PDK converter

## Key Data Structures

**ModelDef** - Represents `.model` with device_type, level, version, and parameters

**Subcircuit** - Contains ports, parameters, instances, and nested models

**Netlist** - Top-level container with models, subcircuits, library sections, and includes

**OsdiModuleInfo** - Maps SPICE device family/level/version to OSDI file and module name
