def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    cell_size = get_config_value(key)
    # Default cell_size is 2.5 (values below scale linearly with it)
    timesteps_per_coupling_cycle = 50 * cell_size / 2.5
    linked_cells_per_coupling_cell = 1 * cell_size / 2.5
    molecules_per_direction = 28 * cell_size / 2.5
    partial_xml.substitute("{cell-size}", cell_size)
    partial_xml.substitute(
        "{timesteps-per-coupling-cycle}", timesteps_per_coupling_cycle
    )
    partial_xml.substitute(
        "{linked-cells-per-coupling-cell}", linked_cells_per_coupling_cell
    )
    partial_xml.substitute("{molecules-per-direction}", molecules_per_direction)
    print(
        "Substituted cell size, timesteps per coupling cycle, linked cells per coupling cell and molecules per direction"
    )
