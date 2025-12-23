# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""SPICE dialect implementations.

This package contains dialect-specific parsing rules for different
SPICE simulators:

    - ngspice: Ngspice simulator (reference implementation)
    - hspice: Synopsys HSPICE (priority dialect)
    - ltspice: Analog Devices LTSpice
    - spectre: Cadence Spectre (industry standard)
    - spectre-spice: Spectre-SPICE mode for mixed-dialect files
"""

from spiceparser.dialects.hspice import HspiceDialect
from spiceparser.dialects.ltspice import LtspiceDialect
from spiceparser.dialects.ngspice import NgspiceDialect
from spiceparser.dialects.spectre import SpectreDialect, SpectreSpiceDialect

__all__ = [
    "NgspiceDialect",
    "HspiceDialect",
    "LtspiceDialect",
    "SpectreDialect",
    "SpectreSpiceDialect",
]
