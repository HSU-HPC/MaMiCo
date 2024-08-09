from generators.utils import get_config_value, get_selected


def apply(partialXml: object, configs: dict) -> None:
    solver = get_config_value(configs, __name__.split(".")[-1])
    partialXml.substitute("{md-solver}", solver)
    wall_velocity = 0.5 if solver == "md" else 1.5
    partialXml.substitute("{wall-velocity}", wall_velocity)
    if solver == "ls1":
        raise NotImplementedError("TODO: Create corresponding ls1config.xml")
    print("Substituted MD solver and wall velocity")
