#! /usr/bin/env python3

"""
Numerically and visually assess the difference between the sampled and analytical couette flow profile.
"""

import math
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except:
    print("Could not import matplotlib.pyplot", file=sys.stderr)
    plt = None

if __name__ != "__main__":
    raise ImportError(f"{__file__} cannot be imported as a module!")

density = None
viscosity = None
wall_velocity = None
channel_height = None
overlap_size = 3
coupling_cell_size = None
coupling_cells = 6
offset = 2.5


def load_values_from_xml():
    """Loads the uninitialized values above from the couette.xml file."""
    global density, viscosity, wall_velocity, channel_height, coupling_cell_size
    tree = ET.parse("couette.xml")
    domain = tree.find("./couette-test/domain")
    channel_height = int(domain.get("channelheight"))
    wall_velocity = float(domain.get("wall-velocity").split(";")[0].strip())
    viscosity = float(tree.find("./couette-test/macroscopic-solver").get("viscosity"))
    density = float(tree.find("./couette-test/microscopic-solver").get("density"))
    cc_conf = tree.find("./mamico/coupling-cell-configuration")
    coupling_cell_size = float(cc_conf.get("cell-size").split(";")[0].strip())


def couette_analytic(z, t):
    """Computes the analytic Couette startup equation at time t and position z."""
    v = viscosity / density
    k_sum = [0] * len(z)
    for k in range(1, 1001):
        for i in range(len(z)):
            k_sum[i] += (
                (1.0 / k)
                * math.sin(k * math.pi * z[i] / channel_height)
                * math.exp(
                    (
                        -1.0
                        * (k * k)
                        * (math.pi * math.pi)
                        * v
                        * t
                        / (channel_height * channel_height)
                    )
                )
            )
    k_sum = [_k_sum * 2.0 * wall_velocity / math.pi for _k_sum in k_sum]
    return [
        wall_velocity * (1.0 - (_z / channel_height)) - _k_sum
        for (_z, _k_sum) in zip(z, k_sum)
    ]


def load_avg_ux_from_csv(csv_file):
    """Get CSV data from one cycle and compute the average velocity per layer of cells in z-direction."""
    csv = []
    with open(csv_file, "r") as file:
        lines = file.readlines()
        if len(lines) >= 2:
            lines = lines[1:]

        def try_parse(s):
            try:
                return float(s)
            except:
                return 0

        csv = [[try_parse(s) for s in line.split(",")] for line in lines]
    # get Avg x velocity per z layer
    col_idx_z = 2
    col_vel_x = 3
    avg_ux = []
    for i in range(4, 10):
        avg = 0
        j = 0
        for row in csv:
            if row[col_idx_z] != i:
                continue
            avg += row[col_vel_x]
            j += 1
        if j > 0:
            avg_ux.append(avg / j)
    return avg_ux


def get_rmspe(expected, actual):
    """Computes the root mean squared percentage error."""
    residual = residual = [e - a for e, a in zip(expected, actual)]
    pe = [100 * r / e for r, e in zip(residual, expected)]
    spe = [e**2 for e in pe]
    mspe = sum(spe) / len(expected)
    rmspe = mspe**0.5
    return rmspe


def compute_couette_flow_profile_match():
    """Computes the sampling error and plots the flow profile for the last cycle."""
    coupling_cycles = []
    for path in Path().glob("CouetteAvgMultiMDCells_0_0_*.csv"):
        coupling_cycle = int(path.stem.split("_")[-1])
        coupling_cycles.append(coupling_cycle)
    coupling_cycle = max(coupling_cycles)
    csv_path = Path(f"CouetteAvgMultiMDCells_0_0_{coupling_cycle}.csv")
    if not csv_path.exists():
        print(f"File {csv_path} does not exist!", file=sys.stderr)
        return
    sampled = load_avg_ux_from_csv(csv_path)
    z = [channel_height / 21 * i for i in range(21)]
    z_full = z
    analytical = couette_analytic(
        z, coupling_cycle / 4
    )  # Multiply MD timestep by number of MD per coupling, multiply that factor with coupling_cycle here
    z_start = (overlap_size * coupling_cell_size) + (coupling_cell_size / 2) + offset
    z_stop = (
        (overlap_size * coupling_cell_size)
        + (coupling_cell_size / 2)
        + (coupling_cell_size * (coupling_cells - 1))
        + offset
    )
    z_step = (z_stop - z_start) / coupling_cells
    z = [z_start + z_step * i for i in range(coupling_cells)]
    analytical = couette_analytic(z_full, coupling_cycle / 4)
    rmspe = get_rmspe(analytical, sampled)
    print(f"=== Couette flow profile from {csv_path} ===")
    print(f"max. coupling cycles = {coupling_cycle}")
    print(f"analytical vs. sampled RMSPE = {rmspe:.3f}%")
    if plt is not None:
        plt.style.use("ggplot")
        plt.plot(z_full, analytical, "C0-", label="Analytical")
        plt.plot(z, sampled, "C0o", label="Sampled")
        plt.ylabel("$u_x$")
        plt.xlabel("$z$")
        plt.grid(True)
        plt.legend()
        plt.title(
            f"Couette startup flow\n(coupling cycle={coupling_cycle}, RMSPE={rmspe:.3f}%)"
        )
        plt.savefig("couette-flow-profile.png")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <directory>", file=sys.stderr)
        exit(1)
    os.chdir(sys.argv[1])
    load_values_from_xml()
    compute_couette_flow_profile_match()
