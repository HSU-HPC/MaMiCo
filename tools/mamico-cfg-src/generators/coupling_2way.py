_min_md_sims_unfiltered = 200
_min_md_sims_filtered = 20


def _instance_shortfall(get_config_value):
    """Return the required MD instance count if there are too few for stable 2-way coupling, else None."""
    key = __name__.split(".")[-1]
    if not get_config_value(key) or get_config_value("simulation") == "test":
        return None
    min_md_sims = (
        _min_md_sims_filtered
        if get_config_value("filtering")
        else _min_md_sims_unfiltered
    )
    if get_config_value("multi_md") < min_md_sims:
        return min_md_sims
    return None


def _message(min_md_sims) -> str:
    return f"Too many coupling cycles or too few MD instances or insufficient filtering for 2-way coupling. ({min_md_sims} MD instances recommended)"


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    Two way coupling is only stable for up to 50 coupling cycles or if averaging over many (>= 200) MD solvers.
    The gaussian filter downgrades this to a warning so that under-sampled runs stay reachable.
    """
    if get_config_value("filtering") == "gauss":
        return None
    min_md_sims = _instance_shortfall(get_config_value)
    return None if min_md_sims is None else _message(min_md_sims)


def warn(get_config_value) -> str:
    if get_config_value("filtering") != "gauss":
        return None
    min_md_sims = _instance_shortfall(get_config_value)
    return None if min_md_sims is None else _message(min_md_sims)


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    # Kept "yes" for preCICE even when one-way: PreciceAdapter only registers MD2macro
    # mesh vertices under twoWayCoupling(), and preCICE rejects an empty provided mesh.
    # The published setup is "technically bidirectional, physically unidirectional"
    # (Neumann & Bian, CPC 220 (2017) 390-402, section 4) -- the OpenFOAM hole patches
    # decide whether the MD data actually reaches the solution, see solver_cfd.py.
    technical = get_config_value(key) or get_config_value("solver_cfd") == "precice"
    partial_xml.substitute("two-way-coupling", "yes" if technical else "no")
    print("Substituted 2-way coupling")
