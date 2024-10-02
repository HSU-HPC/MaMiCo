from utils import get_domain_size


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    The analytical CFD solver does not support two-way coupling.
    Only the OpenFOAM solver should be used with two-way coupling.
    """
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    use_2way_coupling = get_config_value("coupling_2way")
    if solver == "foam" and get_domain_size(get_config_value) != 2:
        return f"OpenFOAM can only be used with the medium domain size."
    if solver == "analytical" and use_2way_coupling:
        return f"Cannot use two-way coupling with analytical CFD solver."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("cfd-solver", solver)
    print("Substituted CFD solver")
    if solver == "foam":
        print('TODO: Please update "foam-setup-directory" in the generated file!')
