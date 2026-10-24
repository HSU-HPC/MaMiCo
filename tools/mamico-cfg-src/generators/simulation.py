from generators.domain_size import get_domain_size

# Cycles between full coupling-cell dumps (CouplingCell_macro_*), 0 disables them.
# The CSV that write2csv emits covers only the md2macro region, so it cannot see
# the Nie imposition layers (2-3) where the continuum velocity is actually applied.
# This is the only output that does, which makes it the instrument for diagnosing
# where the MD-side amplification originates.
_plot_every_macroscopic_timestep = 20


def get_coupling_cycles(get_config_value) -> int:
    key = __name__.split(".")[-1]
    simulation_type = get_config_value(key)
    domain_size = get_domain_size(get_config_value)
    if simulation_type == "test":
        return 50
    elif simulation_type == "short":
        return 250 * domain_size**2
    elif simulation_type == "normal":
        return 500 * domain_size**2
    raise ValueError(key + " = " + simulation_type)


def apply(partial_xml, get_config_value) -> None:
    partial_xml.substitute("coupling-cycles", get_coupling_cycles(get_config_value))
    partial_xml.substitute(
        "plot-every-macroscopic-timestep", _plot_every_macroscopic_timestep
    )
    print("Substituted coupling cycles and macroscopic plot interval")
