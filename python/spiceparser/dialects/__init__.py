# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""SPICE dialect implementations.

This package contains dialect-specific parsing rules for different
SPICE simulators:

    - ngspice: Ngspice simulator (reference implementation)
    - hspice: Synopsys HSPICE (priority dialect)
    - ltspice: Analog Devices LTSpice
"""

from spiceparser.dialects.ngspice import NgspiceDialect
from spiceparser.dialects.hspice import HspiceDialect

# Import dialects to register them
# from spiceparser.dialects.ltspice import LtspiceDialect  # TODO

__all__ = [
    "NgspiceDialect",
    "HspiceDialect",
    # "LtspiceDialect",
]
