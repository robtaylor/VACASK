#!/usr/bin/env python3
"""Convert and run tb_dp512x8 test.

This script:
1. Converts the original ngspice netlist to VACASK format using ng2vc
2. Adds PSP model load and control block
3. Runs the simulation (optional, very slow in debug builds)
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Paths
TEST_DIR = Path(__file__).parent
PYTHON_DIR = TEST_DIR.parent.parent / "python"
SPICE_FILE = TEST_DIR / "tb_dp512x8_klu.spi"
SIM_FILE = TEST_DIR / "tb_dp512x8.sim"
PSP_MODEL_PATH = "./psp103_ihp.osdi"

# IHP PDK path - check common locations
IHP_PDK_CANDIDATES = [
    Path(os.environ.get("IHP_PDK_ROOT", "")) / "libs.tech/ngspice/models",
    Path.home() / "Code/ChipFlow/PDK/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models",
    Path("/foss/pdks/ihp-sg13g2/libs.tech/ngspice/models"),
]


def find_ihp_pdk():
    """Find the IHP PDK models directory."""
    for candidate in IHP_PDK_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def convert_netlist():
    """Convert ngspice netlist to VACASK format."""
    print(f"Converting {SPICE_FILE.name} to VACASK format...")

    ihp_models = find_ihp_pdk()
    if not ihp_models:
        print("Warning: IHP PDK not found. Trying conversion anyway...")
        ihp_models = IHP_PDK_CANDIDATES[0]

    # Create temp file with fixed paths
    with tempfile.NamedTemporaryFile(mode='w', suffix='.spi', delete=False) as tmp:
        content = SPICE_FILE.read_text()
        # Replace container path with local PDK path
        content = content.replace(
            "/foss/pdks/ihp-sg13g2/libs.tech/ngspice/models",
            str(ihp_models)
        )
        tmp.write(content)
        tmp_path = tmp.name

    try:
        cmd = [
            sys.executable, "-m", "ng2vc",
            "-sp", str(ihp_models),
            tmp_path,
            str(SIM_FILE),
        ]

        result = subprocess.run(cmd, cwd=PYTHON_DIR, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Conversion failed: {result.stderr}")
            return False

        print(f"Generated {SIM_FILE.name}")
        return True
    finally:
        Path(tmp_path).unlink()


def add_psp_load_and_control():
    """Add PSP OSDI load and control block to the converted netlist."""
    print("Adding PSP load and control block...")

    content = SIM_FILE.read_text()

    # Add PSP load before capacitor load
    content = content.replace(
        'load"spice/capacitor.osdi"',
        f'load"{PSP_MODEL_PATH}"\nload"spice/capacitor.osdi"'
    )

    # Add control block at the end
    control_block = """
// Control block
control
  analysis tran1 tran step=1e-10 stop=3e-8
  save all
endc
"""
    content += control_block

    SIM_FILE.write_text(content)
    print("Added PSP load and control block")
    return True


def main():
    if not SPICE_FILE.exists():
        print(f"Error: {SPICE_FILE} not found")
        return 1

    ihp_pdk = find_ihp_pdk()
    if not ihp_pdk:
        print("Warning: IHP PDK models not found")
        print("Searched locations:")
        for loc in IHP_PDK_CANDIDATES:
            print(f"  - {loc}")
        print("Set IHP_PDK_ROOT environment variable or adjust paths in script")

    if not convert_netlist():
        return 1

    if not add_psp_load_and_control():
        return 1

    print(f"\nSuccess! Generated {SIM_FILE}")
    print("\nTo run the simulation:")
    print(f"  cd {TEST_DIR}")
    print(f"  vacask {SIM_FILE.name}")
    print("\nNote: Simulation is very slow in debug builds due to circuit complexity.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
