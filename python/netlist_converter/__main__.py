# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Entry point for running netlist_converter as a module.

Usage:
    python -m netlist_converter [<args>] <input file> [<output file>]
"""

import sys
from pathlib import Path

# Add parent directory to path to import netlist_converter.py script
_parent = Path(__file__).parent.parent
if str(_parent) not in sys.path:
    sys.path.insert(0, str(_parent))

# Import main from the script (netlist_converter.py in parent directory)
# This avoids circular import issues
if __name__ == "__main__":
    # Import here to avoid issues with module loading
    import importlib.util

    script_path = _parent / "netlist_converter.py"
    spec = importlib.util.spec_from_file_location("netlist_converter_script", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.exit(module.main())
