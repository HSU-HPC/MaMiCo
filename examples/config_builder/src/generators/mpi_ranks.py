def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    ranks = get_config_value(key)
    ranks_xyz = [1, 1, 1]
    ax = 0
    while ranks > 1:
        ranks_xyz *= 2
        ranks /= 2
        ranks = (ranks + 1) % len(ranks_xyz)
    for ranks, ax in zip(ranks_xyz, ["x", "y", "z"]):
        partial_xml.substitute("{mpi-size-" + ax + "}", ranks)
    print("Substituted MPI size")
