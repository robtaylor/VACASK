# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""SPICE netlist parser with multi-dialect support.

This package provides a unified interface for parsing SPICE netlists from
different simulators (Ngspice, HSPICE, LTSpice), extracting models and
subcircuits for conversion to Verilog-A.

Supported dialects:
    - ngspice: Ngspice simulator netlists
    - hspice: Synopsys HSPICE netlists (priority)
    - ltspice: LTSpice netlists

Usage:
    from spiceparser import parse_netlist
    from spiceparser.dialects import NgspiceDialect

    netlist = parse_netlist("circuit.sp", dialect=NgspiceDialect())
"""

from spiceparser.dialect import SpiceDialect, get_dialect, register_dialect
from spiceparser.netlist import ModelDef, Subcircuit, Instance, Netlist
from spiceparser.parser import parse_netlist

# Import dialects to register them
from spiceparser.dialects import ngspice as _ngspice  # noqa: F401

__all__ = [
    "SpiceDialect",
    "get_dialect",
    "register_dialect",
    "ModelDef",
    "Subcircuit",
    "Instance",
    "Netlist",
    "parse_netlist",
]
