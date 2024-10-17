def get_ranks_xyz(get_config_value):
    key = __name__.split(".")[-1]
    ranks = get_config_value(key)
    ranks_xyz = [1, 1, 1]
    ax = 0
    while ranks > 1:
        ranks_xyz[ax] *= 2
        ranks /= 2
        ax = (ax + 1) % len(ranks_xyz)
    return ranks_xyz


def apply(partial_xml, get_config_value) -> None:
    for ranks, ax in zip(get_ranks_xyz(get_config_value), ["x", "y", "z"]):
        partial_xml.substitute(f"mpi-size-{ax}", ranks)
    print("Substituted MPI size")
