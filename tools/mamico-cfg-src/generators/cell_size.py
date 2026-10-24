from generators.domain_size import get_domain_sizes
from utils import MD_TIMESTEP, get_domain_size

# Other values scale linearly with cell size
_default_cell_size = 2.5

# Coupling time step in MD time units. Wittenberg & Neumann, Computation 9 (2021) 128
# use 0.25 for the validated two-way OpenFOAM setup at cell size 5, and keep it fixed
# across the resolutions they tried. It deliberately does not scale with cell size:
# the Nie upper bound dx^2/nu grows with dx, but the binding constraint is resolving
# the Couette start-up transient, and their sim-2 shows that 4x too large leaves a
# systematic lag the run never recovers from.
_coupling_timestep = 0.25  # default of the coupling_dt option

# CFD solvers that discretise the channel on a grid of the coupling cell size.
# "precice" is absent: its OpenFOAM mesh is generated to fit whatever cell size is chosen.
_grid_cfd_solvers = ("lb", "fd", "foam")


def get_molecules_per_direction(get_config_value) -> float:
    domain_size = get_domain_size(get_config_value)
    size_factor = domain_size
    if domain_size == 3:
        size_factor = 4
    return 28 * size_factor


def get_timesteps_per_coupling_cycle(get_config_value) -> int:
    return round(float(get_config_value("coupling_dt")) / MD_TIMESTEP)


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    The MD domain must be an integer number of coupling cells, and grid-based CFD
    solvers additionally need the channel to be an integer number of coupling cells.
    """
    cell_size = get_config_value("cell_size")
    md_domain_size, cfd_domain_size = get_domain_sizes(get_config_value)
    if md_domain_size % cell_size != 0:
        return f"MD domain of {md_domain_size} is not an integer multiple of cell size {cell_size}."
    if (
        get_config_value("solver_cfd") in _grid_cfd_solvers
        and cfd_domain_size % cell_size != 0
    ):
        return f"Channel height of {cfd_domain_size} is not an integer multiple of cell size {cell_size}."
    return None


def apply(partial_xml, get_config_value) -> None:
    cell_size = get_config_value("cell_size")
    timesteps_per_coupling_cycle = get_timesteps_per_coupling_cycle(get_config_value)
    linked_cells_per_coupling_cell = int(1 * cell_size / _default_cell_size)
    molecules_per_direction = get_molecules_per_direction(get_config_value)
    partial_xml.substitute("cell-size", cell_size)
    partial_xml.substitute("timesteps-per-coupling-cycle", timesteps_per_coupling_cycle)
    partial_xml.substitute(
        "linked-cells-per-coupling-cell", linked_cells_per_coupling_cell
    )
    partial_xml.substitute("molecules-per-direction", molecules_per_direction)
    print(
        "Substituted cell size, timesteps per coupling cycle, linked cells per coupling cell and molecules per direction"
    )
