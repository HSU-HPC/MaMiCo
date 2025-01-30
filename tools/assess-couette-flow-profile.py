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
        max(0, wall_velocity * (1.0 - (_z / channel_height)) - _k_sum)
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

def plot_compute_couette_flow_mpl(analytical, sampled, coupling_cycle=None, rmspe=None):
    if plt is None:
        return False
    plt.style.use("ggplot")
    plt.plot(analytical[0], analytical[1], "C0-", label="Analytical")
    plt.plot(sampled[0], sampled[1], "C0o", label="Sampled")
    plt.ylabel("$u_x$")
    plt.xlabel("$z$")
    plt.grid(True)
    plt.legend()
    info_text = ""
    if coupling_cycle is not None or rmspe is not None:
        if coupling_cycle is None:
            coupling_cycle = "NaN"
        if rmspe is None:
            rmspe = "NaN"
        else:
            rmspe = f"{rmspe:.3f}%"
        info_text = f"\n(coupling cycle={coupling_cycle}, RMSPE={rmspe})"
    plt.title(
        f"Couette startup flow{info_text}"
    )
    plt.ylim(0, wall_velocity)
    plt.savefig("couette-flow-profile.png")
    return True

def plot_compute_couette_flow_mmd(analytical, sampled_z_start_stop, sampled_u_x, coupling_cycle=None, rmspe=None):
    with open("couette-flow-profile.mmd", "w") as file:
        print("xychart-beta", file=file)
        info_text = ""
        if coupling_cycle is not None or rmspe is not None:
            if coupling_cycle is None:
                coupling_cycle = "NaN"
            if rmspe is None:
                rmspe = "NaN"
            else:
                rmspe = f"{rmspe:.3f}%"
            info_text = f" (coupling cycle={coupling_cycle}, RMSPE={rmspe})"
        print(f"\ttitle \"Couette startup flow{info_text}\"", file=file)
        print(f"\tx-axis \"z\" 0 --> {channel_height}", file=file)
        print(f"\ty-axis \"u_x\" 0 --> {wall_velocity}", file=file)
        print("\tbar [", end="", file=file)
        i = 0
        for j, z in enumerate(analytical[0]):
            u_x = -1
            if z >= sampled_z_start_stop[0] and z <= sampled_z_start_stop[1]:
                u_x = sampled_u_x[i]
                i += 1
            if j > 0:
                print(", ", end="", file=file)
            print(u_x, end="", file=file)
        print("]", file=file)
        print("\tline", analytical[1], file=file)
    return True

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
    z = [channel_height / 20 * i for i in range(21)]
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
    print(f"### Couette flow profile from `{csv_path.name}`")
    print(f"max. coupling cycles = {coupling_cycle}")
    print(f"analytical vs. sampled RMSPE = {rmspe:.3f}%")
    plot_compute_couette_flow_mpl([z_full, analytical], [z, sampled], coupling_cycle, rmspe)
    plot_compute_couette_flow_mmd([z_full, analytical], [z_start, z_stop], sampled, coupling_cycle, rmspe)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <couette output directory>", file=sys.stderr)
        exit(1)
    os.chdir(sys.argv[1])
    load_values_from_xml()
    compute_couette_flow_profile_match()
