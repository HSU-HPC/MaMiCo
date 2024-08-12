def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    simulation_type = get_config_value(key)
    domain_size = get_config_value("domain")
    match simulation_type:
        case "test":
            partial_xml.substitute("{coupling-cycles}", 50)
        case "small":
            partial_xml.substitute("{coupling-cycles}", 250 * domain_size**2)
        case "normal":
            partial_xml.substitute("{coupling-cycles}", 500 * domain_size**2)
        case _:
            raise ValueError(key + " = " + simulation_type)

    print("Substituted coupling cycles")
