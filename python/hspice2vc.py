#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""HSPICE to VACASK netlist converter.

Converts HSPICE-format netlists to VACASK format using the multi-dialect
SPICE parser library.

Usage:
    python3 -m hspice2vc input.hspice [output.vc]
    python3 -m hspice2vc --help

Example:
    python3 -m hspice2vc 130BCD_wrapper.hspice models.vc
"""

import argparse
import sys
from pathlib import Path

from spiceparser import parse_netlist
from vcwriter import write_vacask


def main() -> int:
    """Main entry point for hspice2vc converter."""
    parser = argparse.ArgumentParser(
        prog="hspice2vc",
        description="Convert HSPICE netlists to VACASK format.",
        epilog="If no output file is provided, output is written to stdout.",
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input HSPICE netlist file",
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
        "-s",
        "--section",
        default=None,
        help="Library section to extract (if input is a library file)",
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

    try:
        # Parse the HSPICE netlist
        if args.verbose:
            print(f"Parsing {args.input}...", file=sys.stderr)

        netlist = parse_netlist(
            args.input,
            dialect="hspice",
            search_paths=search_paths,
        )

        if args.verbose:
            print(f"  Found {len(netlist.models)} models", file=sys.stderr)
            print(f"  Found {len(netlist.subcircuits)} subcircuits", file=sys.stderr)
            print(
                f"  Found {len(netlist.library_sections)} library sections",
                file=sys.stderr,
            )

        # Generate signature
        signature = f"// Converted from HSPICE by hspice2vc\n// Source: {args.input.name}\n"

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
