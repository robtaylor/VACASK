#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""LTSpice to VACASK netlist converter.

Converts LTSpice-format netlists to VACASK format using the multi-dialect
SPICE parser library.

Usage:
    python3 -m ltspice2vc input.net [output.vc]
    python3 -m ltspice2vc --help

Example:
    python3 -m ltspice2vc LT1001_TA05.net models.vc
"""

import argparse
import sys
from pathlib import Path

from spiceparser import parse_netlist
from vcwriter import write_vacask


def main() -> int:
    """Main entry point for ltspice2vc converter."""
    parser = argparse.ArgumentParser(
        prog="ltspice2vc",
        description="Convert LTSpice netlists to VACASK format.",
        epilog="If no output file is provided, output is written to stdout.",
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input LTSpice netlist file (.net, .cir, .lib)",
    )

    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=None,
        help="Output VACASK file (default: stdout)",
    )

    parser.add_argument(
        "-sp",
        "--sourcepath",
        action="append",
        default=[],
        help="Add a directory to search for .include and .lib files",
    )

    parser.add_argument(
        "--ltspice-lib",
        type=Path,
        default=None,
        help="Path to LTSpice library directory (for resolving LTC.lib, etc.)",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output",
    )

    args = parser.parse_args()

    # Check input file exists
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    # Build search paths
    search_paths = [args.input.parent, Path(".")]
    search_paths.extend(Path(p) for p in args.sourcepath)

    # Add LTSpice library path if provided
    if args.ltspice_lib:
        search_paths.append(args.ltspice_lib)

    try:
        # Parse the LTSpice netlist
        if args.verbose:
            print(f"Parsing {args.input}...", file=sys.stderr)

        netlist = parse_netlist(
            args.input,
            dialect="ltspice",
            search_paths=search_paths,
        )

        if args.verbose:
            print(f"  Found {len(netlist.models)} models", file=sys.stderr)
            print(f"  Found {len(netlist.subcircuits)} subcircuits", file=sys.stderr)
            print(f"  Found {len(netlist.instances)} instances", file=sys.stderr)

        # Generate signature
        signature = f"// Converted from LTSpice by ltspice2vc\n// Source: {args.input.name}\n"

        # Write output
        if args.output:
            # Ensure output directory exists
            args.output.parent.mkdir(parents=True, exist_ok=True)

            write_vacask(netlist, args.output, signature=signature)

            if args.verbose:
                print(f"Wrote output to {args.output}", file=sys.stderr)
        else:
            # Write to stdout
            write_vacask(netlist, sys.stdout, signature=signature)

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            raise
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
