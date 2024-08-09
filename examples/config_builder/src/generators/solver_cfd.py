from generators.utils import get_config_value


def apply(partialXml: object, configs: dict) -> None:
    solver = get_config_value(configs, __name__.split(".")[-1])
    partialXml.substitute("{cfd-solver}", solver)
    print("Substituted CFD solver")
    if solver == "foam":
        print("TODO: Please update \"foam-setup-directory\" in the generated file!")
