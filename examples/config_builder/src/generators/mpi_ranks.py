from generators.utils import get_config_value


def apply(partialXml: object, configs: dict) -> None:
    ranks = get_config_value(configs, __name__.split(".")[-1])
    # TODO multiply by number of MD instances?
    ranks_xyz = [1, 1, 1]
    ax = 0
    while(ranks > 1):
        ranks_xyz *= 2
        ranks /= 2
        ranks = (ranks + 1) % len(ranks_xyz)
    for ranks, ax in zip(ranks_xyz, ["x","y","z"]):
        partialXml.substitute("{mpi-size-" + ax + "}", ranks)
    print("Substituted MPI size")
