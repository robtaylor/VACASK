# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""VACASK output file writer.

Generates VACASK-format files from parsed netlist objects.
"""

import re
from io import StringIO
from pathlib import Path
from typing import TextIO

from spiceparser.netlist import Instance, ModelDef, Netlist, Subcircuit


# SI prefix conversion for VACASK output
_SI_PREFIX_PATTERN = re.compile(r"(\d+\.?\d*|\.\d+)(meg|g|t|mil)", re.IGNORECASE)
_SI_PREFIX_MAP = {
    "meg": "M",
    "g": "G",
    "t": "T",
}


def _convert_si_prefix(match: re.Match) -> str:
    """Convert SI prefix to VACASK format."""
    num = match.group(1)
    prefix = match.group(2).lower()
    if prefix == "mil":
        return f"({num}*25.4e-6)"
    return num + _SI_PREFIX_MAP.get(prefix, prefix)


def format_value(value: str) -> str:
    """Format a parameter value for VACASK output."""
    if isinstance(value, str):
        # Remove curly braces
        if value.startswith("{") and value.endswith("}"):
            value = value[1:-1]
        # Convert SI prefixes
        return _SI_PREFIX_PATTERN.sub(_convert_si_prefix, value)
    return str(value)


def write_vacask(
    netlist: Netlist,
    output: str | Path | TextIO,
    signature: str = "// Converted by ng2vc converter\n",
    columns: int = 80,
) -> None:
    """Write a netlist in VACASK format.

    Args:
        netlist: Parsed netlist to write
        output: Output file path or file object
        signature: Signature comment to add at top
        columns: Maximum line width for formatting
    """
    writer = VacaskWriter(signature=signature, columns=columns)
    content = writer.write(netlist)

    if isinstance(output, (str, Path)):
        Path(output).write_text(content)
    else:
        output.write(content)


class VacaskWriter:
    """VACASK output file writer.

    Converts parsed SPICE netlists to VACASK format, which uses:
    - Spectre-like syntax with parentheses around node lists
    - `load` directives for OSDI modules
    - `model` statements for model definitions
    - `subckt`/`ends` blocks for subcircuits
    """

    def __init__(
        self,
        signature: str = "// Converted by ng2vc converter\n",
        columns: int = 80,
    ):
        self.signature = signature
        self.columns = columns

        # OSDI module mapping (family, level, version) -> (osdi_file, module_name)
        # This should be loaded from configuration
        self._osdi_modules: set[str] = set()

    def write(self, netlist: Netlist) -> str:
        """Write netlist to VACASK format string."""
        buf = StringIO()

        # Write signature
        buf.write(self.signature)
        buf.write("\n")

        # Collect required OSDI modules
        self._collect_osdi_modules(netlist)

        # Write load statements
        for module in sorted(self._osdi_modules):
            buf.write(f'load "{module}"\n')

        if self._osdi_modules:
            buf.write("\n")

        # Write models
        for model in netlist.models:
            self._write_model(buf, model)

        # Write subcircuits
        for subckt in netlist.subcircuits:
            self._write_subcircuit(buf, subckt)

        return buf.getvalue()

    def _collect_osdi_modules(self, netlist: Netlist) -> None:
        """Collect required OSDI modules from netlist."""
        # For now, just track unique device types
        # Real implementation would map to specific OSDI files
        pass

    def _write_model(self, buf: TextIO, model: ModelDef, indent: int = 0) -> None:
        """Write a model definition."""
        prefix = " " * indent

        # Model statement: model <module_name> <model_name>
        # For now, use device_type as module name
        module_name = f"sp_{model.device_type}"
        buf.write(f"{prefix}model {module_name} {model.name}\n")

        # Parameters
        if model.parameters:
            params_str = self._format_params(model.parameters)
            buf.write(f"{prefix}    parameters {params_str}\n")

        buf.write("\n")

    def _write_subcircuit(
        self, buf: TextIO, subckt: Subcircuit, indent: int = 0
    ) -> None:
        """Write a subcircuit definition."""
        prefix = " " * indent

        # Subcircuit header: subckt name (port1 port2 ...)
        ports_str = " ".join(subckt.ports)
        buf.write(f"{prefix}subckt {subckt.name} ({ports_str})\n")

        # Default parameters
        if subckt.parameters:
            for name, value in subckt.parameters.items():
                buf.write(f"{prefix}    parameters {name}={format_value(value)}\n")

        # Nested models
        for model in subckt.models:
            self._write_model(buf, model, indent + 4)

        # Instances
        for inst in subckt.instances:
            self._write_instance(buf, inst, indent + 4)

        buf.write(f"{prefix}ends {subckt.name}\n\n")

    def _write_instance(self, buf: TextIO, inst: Instance, indent: int = 0) -> None:
        """Write a device instance."""
        prefix = " " * indent

        # Instance format: name (node1 node2 ...) model_or_type params...
        nodes_str = " ".join(inst.nodes)

        # Determine type string
        type_str = inst.model_name if inst.model_name else inst.device_type

        # Build instance line
        line = f"{prefix}{inst.name} ({nodes_str}) {type_str}"

        # Add parameters
        if inst.parameters:
            params_str = self._format_params(inst.parameters)
            line += f" {params_str}"

        buf.write(line + "\n")

    def _format_params(self, params: dict[str, any]) -> str:
        """Format parameters for output."""
        parts = []
        for name, value in params.items():
            formatted = format_value(str(value)) if value is not None else ""
            parts.append(f"{name}={formatted}")
        return " ".join(parts)
