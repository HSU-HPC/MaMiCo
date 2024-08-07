#! /usr/bin/env python3

"""
Simple but versitile script to plot couette flow profiles

(Replaces https://github.com/HSU-HPC/MaMiCo/commit/42ad244c75640a692ae1b70d56c7060431fdab0d)
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

args = None


def couette_analytic(z, t):
    """Analytic Couette startup equation"""
    v = args.viscosity / args.density
    k_sum = 0
    for k in range(1, 1001):
        k_sum += (
            (1.0 / k)
            * np.sin(k * np.pi * z / args.channel_height)
            * np.exp(
                (
                    -1.0
                    * (k * k)
                    * (np.pi * np.pi)
                    * v
                    * t
                    / (args.channel_height * args.channel_height)
                )
            )
        )
    k_sum *= 2.0 * args.wall_velocity / np.pi
    return args.wall_velocity * (1.0 - (z / args.channel_height)) - k_sum


def load_avg_ux_from_csv(csv_file):
    """Get CSV data from one cycle and compute
    the average velocity per layer of cells in z-direction
    """
    # load data in pandas DataFrame
    df = pd.read_csv(csv_file)
    # get Avg x velocity per z layer
    avg_ux = []
    idx_col = "I01_z"
    for i in range(4, 10):
        avg = 0
        mass = 0
        j = 0
        for _, row in df[df[idx_col] == i].iterrows():
            avg += row["vel_x"]
        if df[df[idx_col] == i].shape[0] > 0:
            avg_ux.append(avg / df[df[idx_col] == i].shape[0])
    return avg_ux


def plot_couette_profile(coupling_cycle, color, ax=plt.gca()):
    """Plot the flow profile for a single cycle"""
    csv_path = Path(f"CouetteAvgMultiMDCells_0_0_{coupling_cycle}.csv")
    if not csv_path.exists():
        print(f"File {csv_path} does not exist!", file=sys.stderr)
        return
    data = load_avg_ux_from_csv(csv_path)
    z = np.linspace(0, 50, num=21)
    y = couette_analytic(
        z, coupling_cycle / 4
    )  # Multiply MD timestep by number of MD per coupling, multiply that factor with coupling_cycle here
    ax.plot(z, y, "-", color=color)
    x_start = (
        (args.overlap_size * args.coupling_cell_size)
        + (args.coupling_cell_size / 2)
        + args.offset
    )
    x_stop = (
        (args.overlap_size * args.coupling_cell_size)
        + (args.coupling_cell_size / 2)
        + (args.coupling_cell_size * (args.coupling_cells - 1))
        + args.offset
    )
    ax.plot(
        np.linspace(x_start, x_stop, num=args.coupling_cells), data, "o", color=color
    )
    ax.fill_between([], [], [], color=color, label=f"Coupling cycle #{coupling_cycle}")


def parse_args(argv=sys.argv[1:]):
    """Parse all command line arguments and make them available globally"""
    global args
    arg_parser = argparse.ArgumentParser()
    # Scenario parameters
    arg_parser.add_argument("--offset", default=2.5)
    arg_parser.add_argument("--wall-velocity", default=0.5)
    arg_parser.add_argument("--channel-height", default=50.0)
    arg_parser.add_argument("--density", default=0.813037037)
    arg_parser.add_argument("--viscosity", default=2.4)
    arg_parser.add_argument("--coupling-cells", default=6)
    arg_parser.add_argument(
        "--overlap-size", default=3, help="In number of coupling cells"
    )
    arg_parser.add_argument("--coupling-cell-size", default=2.5)
    # Script parameters
    arg_parser.add_argument("--workdir", default=Path(__file__).parent.parent)
    arg_parser.add_argument(
        "--coupling-cycles",
        default="",
        type=str,
        help="Comma separated coupling cycle indices",
    )
    arg_parser.add_argument("--output", default=None, type=Path)
    args = arg_parser.parse_args(argv)


if __name__ == "__main__":
    plt.style.use("ggplot")
    parse_args()
    os.chdir(Path(args.workdir))
    coupling_cycles = [
        int(s.strip()) for s in args.coupling_cycles.strip().split(",") if len(s) > 0
    ]
    if len(coupling_cycles) == 0:
        filename_prefix = "CouetteAvgMultiMDCells_0_0_"
        for path in Path(".").glob(f"{filename_prefix}*.csv"):
            coupling_cycles.append(int(path.stem[len(filename_prefix) :]))
    coupling_cycles.sort()
    for i, cc in enumerate(coupling_cycles):
        plot_couette_profile(cc, f"C{i}")
    plt.plot([], [], "-", color="grey", label="Analytical")
    plt.plot([], [], "o", color="grey", label="Sampled")
    plt.ylabel("$u_x$")
    plt.xlabel("$z$")
    plt.grid(True)
    plt.legend()
    if args.output is not None:
        plt.savefig(args.output, dpi=300)
        print("Saved figure to", Path(args.output).absolute())
    else:
        plt.show()
