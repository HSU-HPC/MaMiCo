#! /usr/bin/env python3

"""
Simple but versatile script to plot couette flow profiles based on the CSV output.

(Replaces https://github.com/HSU-HPC/MaMiCo/commit/42ad244c75640a692ae1b70d56c7060431fdab0d)
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from xml.parsers.expat import errors as xmlErrors
import xml.etree.cElementTree as ET


def couette_analytic(z, t, args):
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
    for i in range(df[idx_col].min(), df[idx_col].max() + 1):
        avg = 0
        for _, row in df[df[idx_col] == i].iterrows():
            avg += row["vel_x"]
        if df[df[idx_col] == i].shape[0] > 0:
            avg_ux.append(avg / df[df[idx_col] == i].shape[0])
    return avg_ux


def plot_couette_profile(coupling_cycle, color, args, ax=plt.gca()):
    """Plot the flow profile for a single cycle"""
    csv_path = Path(f"CouetteAvgMultiMDCells_d0_r0_c{coupling_cycle}.csv")
    if not csv_path.exists():
        print(f"File {csv_path} does not exist!", file=sys.stderr)
        return
    data = load_avg_ux_from_csv(csv_path)
    z = np.linspace(
        0,
        args.channel_height,
        num=int(args.channel_height / args.coupling_cell_size) + 1,
    )
    y = couette_analytic(
        z, coupling_cycle * args.md_ts_per_cc * args.md_ts_length, args
    )
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


def parse_xml(args, defaults):
    # Load file
    try:
        tree = ET.parse(args.couette_xml)
        tree = tree.getroot()
    except ET.ParseError as err:
        if (
            err.code == xmlErrors.codes[xmlErrors.XML_ERROR_JUNK_AFTER_DOC_ELEMENT]
        ):  # legacy style xml
            with open(args.couette_xml) as file:
                xmlList = file.readlines()
                # add a root node
                xmlList.insert(1, "<scenario-configuration>")
                xmlList.append("</scenario-configuration>")
                tree = ET.fromstringlist(xmlList)
        else:
            print(err)
            exit(1)

    # Update all needed variables
    if args.offset == defaults.offset:
        args.offset = float(
            tree.find("molecular-dynamics/domain-configuration")
            .attrib["domain-offset"]
            .split(";")[2]
        )
    if args.wall_velocity == defaults.wall_velocity:
        args.wall_velocity = float(
            tree.find("couette-test/domain").attrib["wall-velocity"].split(";")[0]
        )
    if args.channel_height == defaults.channel_height:
        args.channel_height = float(
            tree.find("couette-test/domain").attrib["channelheight"]
        )
    if args.density == defaults.density:
        args.density = float(
            tree.find("couette-test/microscopic-solver").attrib["density"]
        )
    if args.viscosity == defaults.viscosity:
        args.viscosity = float(
            tree.find("couette-test/macroscopic-solver").attrib["viscosity"]
        )
    if args.overlap_size == defaults.overlap_size:
        args.overlap_size = int(
            tree.find("mamico/momentum-insertion").attrib["innermost-overlap-layer"]
        )
    if args.coupling_cell_size == defaults.coupling_cell_size:
        args.coupling_cell_size = float(
            tree.find("mamico/coupling-cell-configuration")
            .attrib["cell-size"]
            .split(";")[0]
        )
    if args.md_ts_per_cc == defaults.md_ts_per_cc:
        args.md_ts_per_cc = int(
            tree.find("molecular-dynamics/simulation-configuration").attrib[
                "number-of-timesteps"
            ]
        )
    if args.md_ts_length == defaults.md_ts_length:
        args.md_ts_length = float(
            tree.find("molecular-dynamics/simulation-configuration").attrib["dt"]
        )
    if args.coupling_cells == defaults.coupling_cells:
        temp_mdSize = float(
            tree.find("molecular-dynamics/domain-configuration")
            .attrib["domain-size"]
            .split(";")[2]
        )
        temp_numCellsInMD = int(temp_mdSize / args.coupling_cell_size)
        args.coupling_cells = temp_numCellsInMD - (2 * args.overlap_size)
    return args


def parse_args(argv=sys.argv[1:]):
    """Parse all command line arguments"""
    arg_parser = argparse.ArgumentParser()
    # Scenario parameters
    arg_parser.add_argument("--offset", default=2.5, type=float)
    arg_parser.add_argument("--wall-velocity", default=0.5, type=float)
    arg_parser.add_argument("--channel-height", default=50.0, type=float)
    arg_parser.add_argument("--density", default=0.813037037, type=float)
    arg_parser.add_argument("--viscosity", default=2.14, type=float)
    arg_parser.add_argument(
        "--coupling-cells",
        default=6,
        type=int,
        help="Number of coupling cells in the z direction in the CSV file (MD cells minus overlap and ghost)",
    )
    arg_parser.add_argument(
        "--overlap-size", default=3, help="In number of coupling cells", type=int
    )
    arg_parser.add_argument("--coupling-cell-size", default=2.5, type=float)
    arg_parser.add_argument(
        "--md-ts-per-cc",
        default=50,
        type=int,
        help="Number of MD timesteps per coupling step",
    )
    arg_parser.add_argument(
        "--md-ts-length", default=0.005, type=float, help="Length of MD timestep"
    )
    # Script parameters
    arg_parser.add_argument(
        "--couette-xml", default="", type=Path, help="Path to couette.xml file"
    )
    arg_parser.add_argument("--workdir", default=Path(), type=Path)
    arg_parser.add_argument(
        "--coupling-cycles",
        default="",
        type=str,
        help="Comma separated coupling cycle indices",
    )
    arg_parser.add_argument("--output", default=None, type=Path)
    args = arg_parser.parse_args(argv)
    defaults = arg_parser.parse_args([])
    return args, defaults


if __name__ == "__main__":
    plt.style.use("ggplot")
    args, defaults = parse_args()
    if args.couette_xml != "":
        args = parse_xml(args, defaults)
    os.chdir(Path(args.workdir))
    coupling_cycles = [
        int(s.strip()) for s in args.coupling_cycles.strip().split(",") if len(s) > 0
    ]
    if len(coupling_cycles) == 0:
        filename_prefix = "CouetteAvgMultiMDCells_d0_r0_c"
        for path in Path(".").glob(f"{filename_prefix}*.csv"):
            coupling_cycles.append(int(path.stem[len(filename_prefix) :]))
    coupling_cycles.sort()
    for i, cc in enumerate(coupling_cycles):
        plot_couette_profile(cc, f"C{i}", args)
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
