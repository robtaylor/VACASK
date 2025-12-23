# VACASK Multi-Dialect SPICE Parser

This package provides a multi-dialect SPICE netlist parser for VACASK, supporting:

- **Ngspice** - Open-source SPICE3 derivative (reference implementation)
- **HSPICE** - Synopsys HSPICE (priority)
- **LTSpice** - Analog Devices LTSpice

## Structure

- `spiceparser/` - Unified SPICE parsing library
  - `dialect.py` - Base dialect class and registry
  - `dialects/` - Dialect-specific implementations
  - `parser.py` - Main parser
  - `netlist.py` - Data structures
- `vcwriter/` - VACASK output generation
- `ng2vclib/` - Original ng2vc library (being refactored)

## Usage

```python
from spiceparser import parse_netlist
from spiceparser.dialects import NgspiceDialect

netlist = parse_netlist("circuit.sp", dialect=NgspiceDialect())

for model in netlist.models:
    print(f"Model: {model.name} ({model.device_type})")
```

## Testing

```bash
uv run pytest tests/ -v
```
