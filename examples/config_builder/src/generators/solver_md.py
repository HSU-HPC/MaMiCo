def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("{md-solver}", solver)
    wall_velocity = 0.5 if solver == "md" else 1.5
    partial_xml.substitute("{wall-velocity}", wall_velocity)
    if solver == "ls1":
        raise NotImplementedError("TODO: Create corresponding ls1config.xml")
    print("Substituted MD solver and wall velocity")
