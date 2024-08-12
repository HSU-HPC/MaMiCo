def validate(get_config_value) -> str:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    use_2way_coupling = get_config_value("coupling_2way")
    if solver == "analytical" and use_2way_coupling:
        return f"Cannot use 2-way coupling with analytical CFD solver."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("cfd-solver", solver)
    print("Substituted CFD solver")
    if solver == "foam":
        print('TODO: Please update "foam-setup-directory" in the generated file!')
