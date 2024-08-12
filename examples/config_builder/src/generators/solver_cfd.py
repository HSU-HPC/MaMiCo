def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("{cfd-solver}", solver)
    print("Substituted CFD solver")
    if solver == "foam":
        print('TODO: Please update "foam-setup-directory" in the generated file!')
