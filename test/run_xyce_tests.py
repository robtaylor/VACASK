#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 ChipFlow
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Run Xyce regression tests with VACASK.

This script converts Xyce test netlists to VACASK format using netlist_converter,
runs VACASK on them, and compares results with Xyce reference output.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

# Add python directory to path for imports
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent
PYTHON_DIR = REPO_ROOT / "python"
sys.path.insert(0, str(PYTHON_DIR))


@dataclass
class TestResult:
    """Result of a single test."""

    name: str
    passed: bool
    error: str | None = None
    max_diff: float | None = None


# Device models supported by VACASK (based on available OSDI modules)
SUPPORTED_MODELS = {
    "bsim4",
    "bsim4_v4p7",
    "bsim4_v4p82",
    # Add more as VACASK supports them
}

# Xyce test directories that use supported models
SUPPORTED_TEST_DIRS = [
    "BSIM4",
    "BSIM4_v4p7",
    "BSIM4_v4p82",
]


def find_vacask_binary() -> Path | None:
    """Find the VACASK binary."""
    # Check in build directory
    build_vacask = REPO_ROOT / "build" / "simulator" / "vacask"
    if build_vacask.exists():
        return build_vacask

    # Check if in PATH
    result = subprocess.run(
        ["which", "vacask"], capture_output=True, text=True, check=False
    )
    if result.returncode == 0:
        return Path(result.stdout.strip())

    return None


def find_xyce_tests(test_dirs: list[str] | None = None) -> list[Path]:
    """Find Xyce test netlists in specified directories."""
    xyce_regression = SCRIPT_DIR / "Xyce_Regression" / "Netlists"

    if test_dirs is None:
        test_dirs = SUPPORTED_TEST_DIRS

    tests = []
    for test_dir in test_dirs:
        dir_path = xyce_regression / test_dir
        if dir_path.exists():
            for cir_file in dir_path.glob("*.cir"):
                tests.append(cir_file)

    return sorted(tests)


def convert_netlist(input_path: Path, output_path: Path) -> bool:
    """Convert Xyce netlist to VACASK format using netlist_converter."""
    try:
        # Use the netlist_converter module
        from netlist_converter import Converter

        converter = Converter()
        converter.convert(str(input_path), str(output_path))
        return True
    except Exception as e:
        print(f"  Conversion error: {e}")
        return False


def run_vacask(vacask_binary: Path, netlist: Path, work_dir: Path) -> tuple[bool, str]:
    """Run VACASK on a netlist."""
    try:
        result = subprocess.run(
            [str(vacask_binary), str(netlist)],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if result.returncode != 0:
            return False, result.stderr or result.stdout
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


def parse_output_file(path: Path) -> dict[str, list[float]]:
    """Parse VACASK or Xyce output file into columns."""
    columns: dict[str, list[float]] = {}

    with open(path) as f:
        lines = f.readlines()

    # Find header line (contains column names)
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Index") or "v-sweep" in line.lower():
            header_idx = i
            break

    if header_idx is None:
        return columns

    # Parse header
    header_line = lines[header_idx]
    col_names = header_line.split()

    # Initialize columns
    for name in col_names:
        columns[name] = []

    # Parse data lines
    for line in lines[header_idx + 2 :]:  # Skip header and separator
        line = line.strip()
        if not line or line.startswith("-"):
            continue

        values = line.split()
        if len(values) >= len(col_names):
            for i, name in enumerate(col_names):
                try:
                    columns[name].append(float(values[i]))
                except ValueError:
                    pass

    return columns


def compare_outputs(
    vacask_output: Path, reference_output: Path, rel_tol: float = 1e-3
) -> tuple[bool, float]:
    """Compare VACASK output with Xyce reference output."""
    vacask_data = parse_output_file(vacask_output)
    ref_data = parse_output_file(reference_output)

    if not vacask_data or not ref_data:
        return False, float("inf")

    max_diff = 0.0

    # Compare common columns
    for col_name in vacask_data:
        if col_name in ref_data:
            vacask_vals = vacask_data[col_name]
            ref_vals = ref_data[col_name]

            if len(vacask_vals) != len(ref_vals):
                continue

            for v, r in zip(vacask_vals, ref_vals):
                if r != 0:
                    diff = abs((v - r) / r)
                else:
                    diff = abs(v)
                max_diff = max(max_diff, diff)

    passed = max_diff <= rel_tol
    return passed, max_diff


def run_test(
    test_path: Path, vacask_binary: Path, work_dir: Path, verbose: bool = False
) -> TestResult:
    """Run a single Xyce test with VACASK."""
    test_name = f"{test_path.parent.name}/{test_path.name}"

    if verbose:
        print(f"Running: {test_name}")

    # Create test work directory
    test_work_dir = work_dir / test_path.parent.name / test_path.stem
    test_work_dir.mkdir(parents=True, exist_ok=True)

    # Copy any include files
    for inc_file in test_path.parent.glob("*.nmos"):
        (test_work_dir / inc_file.name).write_text(inc_file.read_text())
    for inc_file in test_path.parent.glob("*.pmos"):
        (test_work_dir / inc_file.name).write_text(inc_file.read_text())

    # Convert netlist
    converted_path = test_work_dir / f"{test_path.stem}.vc"
    if not convert_netlist(test_path, converted_path):
        return TestResult(test_name, False, "Conversion failed")

    # Run VACASK
    success, error = run_vacask(vacask_binary, converted_path, test_work_dir)
    if not success:
        return TestResult(test_name, False, f"VACASK failed: {error}")

    # Find reference output
    ref_output = (
        SCRIPT_DIR
        / "Xyce_Regression"
        / "OutputData"
        / test_path.parent.name
        / f"{test_path.name}.prn"
    )
    if not ref_output.exists():
        return TestResult(test_name, False, "No reference output")

    # Find VACASK output (TODO: determine actual output filename)
    vacask_output = test_work_dir / f"{test_path.stem}.prn"
    if not vacask_output.exists():
        # Try other common output names
        for pattern in ["*.prn", "*.csv", "*.out"]:
            outputs = list(test_work_dir.glob(pattern))
            if outputs:
                vacask_output = outputs[0]
                break

    if not vacask_output.exists():
        return TestResult(test_name, False, "No VACASK output")

    # Compare outputs
    passed, max_diff = compare_outputs(vacask_output, ref_output)

    if passed:
        return TestResult(test_name, True, max_diff=max_diff)
    else:
        return TestResult(test_name, False, f"Max diff: {max_diff:.2e}", max_diff)


def main():
    parser = argparse.ArgumentParser(description="Run Xyce regression tests with VACASK")
    parser.add_argument(
        "--test-dir",
        action="append",
        help="Specific test directory to run (e.g., BSIM4)",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=SCRIPT_DIR / "xyce_test_results",
        help="Working directory for test outputs",
    )
    parser.add_argument(
        "--vacask",
        type=Path,
        help="Path to VACASK binary",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available tests without running",
    )

    args = parser.parse_args()

    # Find VACASK binary
    vacask_binary = args.vacask or find_vacask_binary()
    if not args.list and not vacask_binary:
        print("Error: VACASK binary not found. Use --vacask to specify path.")
        sys.exit(1)

    # Find tests
    test_dirs = args.test_dir if args.test_dir else None
    tests = find_xyce_tests(test_dirs)

    if args.list:
        print(f"Found {len(tests)} tests:")
        for test in tests:
            print(f"  {test.parent.name}/{test.name}")
        return

    print(f"Found {len(tests)} tests")
    print(f"VACASK binary: {vacask_binary}")
    print(f"Work directory: {args.work_dir}")
    print()

    # Create work directory
    args.work_dir.mkdir(parents=True, exist_ok=True)

    # Run tests
    results: list[TestResult] = []
    passed = 0
    failed = 0

    for test in tests:
        result = run_test(test, vacask_binary, args.work_dir, args.verbose)
        results.append(result)

        if result.passed:
            passed += 1
            status = "PASS"
        else:
            failed += 1
            status = "FAIL"

        if args.verbose or not result.passed:
            print(f"  [{status}] {result.name}")
            if result.error:
                print(f"         {result.error}")

    print()
    print(f"Results: {passed} passed, {failed} failed out of {len(tests)}")


if __name__ == "__main__":
    main()
