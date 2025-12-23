#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "matplotlib",
# ]
# ///
"""Plot dual-port SRAM simulation results.

This script reads VACASK .raw output files and generates plots showing
the behavior of the dual-port SRAM during read/write operations.

Usage:
    uv run plot.py tb_dp512x8_klu.raw          # Plot from raw file
    uv run plot.py tb_dp512x8_klu.raw --save   # Save plots to PNG
    uv run plot.py tb_dp512x8_klu.raw --port 1 # Plot only port 1
    uv run plot.py --list tb_dp512x8_klu.raw   # List available signals
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def rawread(fname):
    """Read a VACASK/SPICE raw file.

    Returns a dict with:
        - 'names': list of signal names
        - 'data': 2D numpy array (npoints x nvars)
        - 'title': simulation title
        - 'plotname': analysis type
    """
    BSIZE_SP = 512
    header_entries = {
        b'title', b'date', b'plotname', b'flags', b'no. variables',
        b'no. points', b'dimensions', b'command', b'option'
    }

    with open(fname, 'rb') as fp:
        plot = {}
        while True:
            try:
                line = fp.readline(BSIZE_SP)
                if not line:
                    break
                split_line = line.split(b':', maxsplit=1)
            except Exception:
                raise RuntimeError("Failed to read a line from file.")

            if len(split_line) == 2:
                key = split_line[0].lower()
                if key in header_entries:
                    plot[key.decode('ascii')] = split_line[1].strip()

                if key == b'variables':
                    nvars = int(plot['no. variables'])
                    npoints = int(plot['no. points'])
                    plot['varnames'] = []
                    plot['varunits'] = []
                    for ii in range(nvars):
                        txt = fp.readline(BSIZE_SP).strip().decode('ascii')
                        var_desc = txt.split(maxsplit=3)
                        assert ii == int(var_desc[0])
                        plot['varnames'].append(var_desc[1])
                        if len(var_desc) > 2:
                            plot['varunits'].append(var_desc[2])
                        else:
                            plot['varunits'].append('')

                if key == b'binary':
                    flags = plot.get('flags', b'')
                    if isinstance(flags, bytes):
                        flags = flags.decode('ascii')
                    dtype = np.complex128 if 'complex' in flags else np.float64
                    arr = np.fromfile(
                        fp, dtype=dtype, count=npoints * nvars
                    ).reshape((npoints, nvars))
                    plot['data'] = arr
                    break
            else:
                break

    return {
        'names': plot.get('varnames', []),
        'units': plot.get('varunits', []),
        'data': plot.get('data', np.array([])),
        'title': plot.get('title', b'').decode('ascii') if isinstance(plot.get('title'), bytes) else plot.get('title', ''),
        'plotname': plot.get('plotname', b'').decode('ascii') if isinstance(plot.get('plotname'), bytes) else plot.get('plotname', ''),
    }


def find_signal(names, pattern):
    """Find signals matching a pattern (case-insensitive)."""
    regex = re.compile(pattern, re.IGNORECASE)
    matches = []
    for i, name in enumerate(names):
        if regex.search(name):
            matches.append((i, name))
    return matches


def get_signal(raw, name):
    """Get signal data by name.

    Handles both direct names and v(name) format.
    """
    # Try direct match first
    if name in raw['names']:
        idx = raw['names'].index(name)
        return raw['data'][:, idx]

    # Try with v() wrapper
    v_name = f'v({name})'
    if v_name in raw['names']:
        idx = raw['names'].index(v_name)
        return raw['data'][:, idx]

    # Try lowercase variants
    for variant in [name, v_name, name.lower(), f'v({name.lower()})']:
        for i, n in enumerate(raw['names']):
            if n.lower() == variant.lower():
                return raw['data'][:, i]

    return None


def get_time(raw):
    """Get the time vector."""
    # Try common names for time/sweep variable
    for name in ['time', 'TIME', 'Time', 'sweep']:
        sig = get_signal(raw, name)
        if sig is not None:
            return sig.real if np.iscomplexobj(sig) else sig
    # Fall back to first column
    return raw['data'][:, 0].real if np.iscomplexobj(raw['data'][:, 0]) else raw['data'][:, 0]


def plot_port_activity(raw, port: int, ax, title_prefix=""):
    """Plot activity for one port of the dual-port SRAM."""
    time = get_time(raw) * 1e9  # Convert to ns
    names = raw['names']

    # Find clock signal - try various formats
    clk = get_signal(raw, f'clk{port}')
    clk_name = f'clk{port}'

    # Find address signals (a1[0] .. a1[8] or a2[0] .. a2[8])
    addr_signals = find_signal(names, f'a{port}\\[')

    # Find data input signals (d1[0] .. d1[7] or d2[0] .. d2[7])
    data_signals = find_signal(names, f'd{port}\\[')

    # Find output signals (q1[0] .. q1[7] or q2[0] .. q2[7])
    q_signals = find_signal(names, f'q{port}\\[')

    # Find write enable
    we_signals = find_signal(names, f'we{port}')

    # Plot
    plot_idx = 0
    signals_plotted = []

    if clk is not None:
        data = clk.real if np.iscomplexobj(clk) else clk
        ax.plot(time, data + plot_idx * 2, label=clk_name, linewidth=1)
        signals_plotted.append(clk_name)
        plot_idx += 1

    # Plot first 2 address bits
    for idx, name in addr_signals[:2]:
        data = raw['data'][:, idx]
        data = data.real if np.iscomplexobj(data) else data
        # Extract short name from v(a1[0]) -> a1[0]
        short_name = name.replace('v(', '').replace(')', '')
        ax.plot(time, data + plot_idx * 2, label=short_name, linewidth=1)
        signals_plotted.append(short_name)
        plot_idx += 1

    # Plot write enable
    for idx, name in we_signals[:1]:
        data = raw['data'][:, idx]
        data = data.real if np.iscomplexobj(data) else data
        short_name = name.replace('v(', '').replace(')', '')
        ax.plot(time, data + plot_idx * 2, label=short_name, linewidth=1)
        signals_plotted.append(short_name)
        plot_idx += 1

    # Plot first 2 data bits
    for idx, name in data_signals[:2]:
        data = raw['data'][:, idx]
        data = data.real if np.iscomplexobj(data) else data
        short_name = name.replace('v(', '').replace(')', '')
        ax.plot(time, data + plot_idx * 2, label=short_name, linewidth=1)
        signals_plotted.append(short_name)
        plot_idx += 1

    # Plot first 2 output bits
    for idx, name in q_signals[:2]:
        data = raw['data'][:, idx]
        data = data.real if np.iscomplexobj(data) else data
        short_name = name.replace('v(', '').replace(')', '')
        ax.plot(time, data + plot_idx * 2, label=short_name, linewidth=1)
        signals_plotted.append(short_name)
        plot_idx += 1

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Voltage (V) + offset')
    ax.set_title(f'{title_prefix}Port {port} Activity')
    if signals_plotted:
        ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_overview(raw, ax_list, title=""):
    """Create overview plots for both ports."""
    for i, port in enumerate([1, 2]):
        if i < len(ax_list):
            plot_port_activity(raw, port, ax_list[i], title_prefix=f"{title} - " if title else "")


def plot_zoom(raw, ax, port: int, t_start: float, t_end: float):
    """Plot zoomed view of read operation."""
    time = get_time(raw) * 1e9
    names = raw['names']

    # Filter to time range
    mask = (time >= t_start) & (time <= t_end)
    time_zoom = time[mask]

    # Find key signals for zoom
    signals_to_plot = []

    # Clock
    clk = get_signal(raw, f'clk{port}')
    if clk is not None:
        signals_to_plot.append((f'clk{port}', clk[mask]))

    # Address signals
    for match in find_signal(names, f'a{port}\\[')[:3]:
        idx, name = match
        short_name = name.replace('v(', '').replace(')', '')
        signals_to_plot.append((short_name, raw['data'][mask, idx]))

    # Write enable
    for match in find_signal(names, f'we{port}')[:1]:
        idx, name = match
        short_name = name.replace('v(', '').replace(')', '')
        signals_to_plot.append((short_name, raw['data'][mask, idx]))

    # Data inputs
    for match in find_signal(names, f'd{port}\\[')[:2]:
        idx, name = match
        short_name = name.replace('v(', '').replace(')', '')
        signals_to_plot.append((short_name, raw['data'][mask, idx]))

    # Output signals
    for match in find_signal(names, f'q{port}\\[')[:2]:
        idx, name = match
        short_name = name.replace('v(', '').replace(')', '')
        signals_to_plot.append((short_name, raw['data'][mask, idx]))

    # Plot with offsets
    for i, (name, data) in enumerate(signals_to_plot):
        data = data.real if np.iscomplexobj(data) else data
        ax.plot(time_zoom, data + i * 2, label=name, linewidth=1)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Voltage (V) + offset')
    ax.set_title(f'Port {port} Read Operation Detail ({t_start:.0f}-{t_end:.0f} ns)')
    if signals_to_plot:
        ax.legend(loc='upper right', fontsize=7)
    ax.grid(True, alpha=0.3)


def list_signals(raw):
    """List all available signals."""
    print(f"\nTitle: {raw['title']}")
    print(f"Analysis: {raw['plotname']}")
    print(f"\nAvailable signals ({len(raw['names'])}):")
    for i, (name, unit) in enumerate(zip(raw['names'], raw['units'])):
        print(f"  {i:4d}: {name} [{unit}]")


def main():
    parser = argparse.ArgumentParser(
        description="Plot dual-port SRAM simulation results"
    )
    parser.add_argument(
        "rawfile",
        help="Path to .raw simulation output file",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save plots to PNG files instead of displaying",
    )
    parser.add_argument(
        "--port",
        type=int,
        choices=[1, 2],
        help="Plot only specified port (1 or 2)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available signals and exit",
    )
    parser.add_argument(
        "--zoom",
        nargs=2,
        type=float,
        metavar=('START', 'END'),
        help="Zoom to time range (in ns)",
    )
    parser.add_argument(
        "--signals",
        nargs='+',
        help="Plot specific signals by name or regex pattern",
    )

    args = parser.parse_args()

    rawfile = Path(args.rawfile)
    if not rawfile.exists():
        print(f"Error: {rawfile} not found")
        return 1

    print(f"Reading {rawfile}...")
    raw = rawread(str(rawfile))

    if not raw['names']:
        print("Error: No data found in raw file")
        return 1

    print(f"Loaded {len(raw['names'])} signals, {raw['data'].shape[0]} points")

    if args.list:
        list_signals(raw)
        return 0

    # Create output filename base
    basename = rawfile.stem

    # Custom signal plotting
    if args.signals:
        fig, ax = plt.subplots(figsize=(12, 8))
        time = get_time(raw) * 1e9

        for i, pattern in enumerate(args.signals):
            matches = find_signal(raw['names'], pattern)
            for idx, name in matches:
                data = raw['data'][:, idx]
                data = data.real if np.iscomplexobj(data) else data
                ax.plot(time, data + i * 2, label=name, linewidth=1)

        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Voltage (V) + offset')
        ax.set_title(f'{basename} - Custom Signals')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)

        if args.save:
            outfile = rawfile.parent / f"{basename}_custom.png"
            plt.savefig(outfile, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile}")
        else:
            plt.show()
        return 0

    # Zoom plot
    if args.zoom:
        t_start, t_end = args.zoom
        port = args.port or 1

        fig, ax = plt.subplots(figsize=(12, 8))
        plot_zoom(raw, ax, port, t_start, t_end)

        if args.save:
            outfile = rawfile.parent / f"{basename}_port{port}_zoom.png"
            plt.savefig(outfile, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile}")
        else:
            plt.show()
        return 0

    # Standard port plots
    if args.port:
        fig, ax = plt.subplots(figsize=(12, 6))
        plot_port_activity(raw, args.port, ax, title_prefix=f"{basename} - ")

        if args.save:
            outfile = rawfile.parent / f"{basename}_port{args.port}.png"
            plt.savefig(outfile, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile}")
        else:
            plt.show()
    else:
        # Plot both ports
        fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
        fig.suptitle(f'{basename} - Dual-Port SRAM Simulation', fontsize=14)

        plot_port_activity(raw, 1, axes[0])
        plot_port_activity(raw, 2, axes[1])

        plt.tight_layout()

        if args.save:
            outfile = rawfile.parent / f"{basename}_overview.png"
            plt.savefig(outfile, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile}")

            # Also create zoom plots for read operations
            # Port 1: read happens around 20-25ns
            fig2, ax2 = plt.subplots(figsize=(12, 8))
            plot_zoom(raw, ax2, 1, 20, 25)
            outfile2 = rawfile.parent / f"{basename}_port1_readzoom.png"
            plt.savefig(outfile2, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile2}")

            # Port 2: similar timing
            fig3, ax3 = plt.subplots(figsize=(12, 8))
            plot_zoom(raw, ax3, 2, 20, 25)
            outfile3 = rawfile.parent / f"{basename}_port2_readzoom.png"
            plt.savefig(outfile3, dpi=150, bbox_inches='tight')
            print(f"Saved {outfile3}")
        else:
            plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
