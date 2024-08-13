# Other values scale linearly with cell size
_default_cell_size = 2.5


def get_cell_size(get_config_value) -> float:
    key = __name__.split(".")[-1]
    return get_config_value(key)


def get_molecules_per_direction(get_config_value) -> float:
    return 28 * get_cell_size(get_config_value) / _default_cell_size


def apply(partial_xml, get_config_value) -> None:
    cell_size = get_cell_size(get_config_value)
    timesteps_per_coupling_cycle = 50 * cell_size / _default_cell_size
    linked_cells_per_coupling_cell = 1 * cell_size / _default_cell_size
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
