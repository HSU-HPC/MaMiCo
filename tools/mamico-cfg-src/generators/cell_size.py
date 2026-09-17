from utils import get_domain_size

# Other values scale linearly with cell size
_default_cell_size = 2.5


def get_molecules_per_direction(get_config_value) -> float:
    domain_size = get_domain_size(get_config_value)
    size_factor = domain_size
    if domain_size == 3:
        size_factor = 4
    return 28 * size_factor


def get_linked_cells_per_coupling_cell(get_config_value) -> float:
    """Number of linked cells along each axis of a coupling cell.

    The LAMMPS adapter embeds exactly one linked cell per coupling cell
    (MamicoLammpsMDSolverInterface::getLinkedCell ignores linkedCellInCouplingCell and
    returns the coupling cell itself), so asking for more makes every per-cell mapping
    traverse a block of neighbouring coupling cells instead of one cell.
    """
    if get_config_value("solver_md") == "lammps-md":
        return 1
    return get_config_value("cell_size") / _default_cell_size


def apply(partial_xml, get_config_value) -> None:
    cell_size = get_config_value("cell_size")
    timesteps_per_coupling_cycle = 50 * cell_size / _default_cell_size
    linked_cells_per_coupling_cell = get_linked_cells_per_coupling_cell(get_config_value)
    # derive the linked cell size from the two values it has to be consistent with,
    # rather than assuming it stays at the default cell size
    linked_cell_size = cell_size / linked_cells_per_coupling_cell
    molecules_per_direction = get_molecules_per_direction(get_config_value)
    partial_xml.substitute("cell-size", cell_size)
    partial_xml.substitute("timesteps-per-coupling-cycle", timesteps_per_coupling_cycle)
    partial_xml.substitute(
        "linked-cells-per-coupling-cell", linked_cells_per_coupling_cell
    )
    partial_xml.substitute("linked-cell-size", linked_cell_size)
    partial_xml.substitute("molecules-per-direction", molecules_per_direction)
    print(
        "Substituted cell size, timesteps per coupling cycle, linked cells per coupling cell, linked cell size and molecules per direction"
    )
