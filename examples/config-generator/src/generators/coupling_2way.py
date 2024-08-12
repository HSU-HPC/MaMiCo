def validate(get_config_value) -> str:
    key = __name__.split(".")[-1]
    use_2way_coupling = get_config_value(key)
    simulation_type = get_config_value("simulation")
    num_md_sims = get_config_value("multi_md")
    if use_2way_coupling and simulation_type != "test" and num_md_sims < 200:
        return f"Too many coupling cycles or too few MD instances for 2-way coupling."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    value = "yes" if get_config_value(key) else "no"
    partial_xml.substitute("two-way-coupling", value)
    print("Substituted 2-way coupling")
