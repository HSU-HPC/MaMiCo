from generators.domain_size import get_domain_size


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    simulation_type = get_config_value(key)
    domain_size = get_domain_size(get_config_value)
    if simulation_type == "test":
            partial_xml.substitute("coupling-cycles", 50)
    elif simulation_type == "small":
            partial_xml.substitute("coupling-cycles", 250 * domain_size**2)
    elif simulation_type == "normal":
            partial_xml.substitute("coupling-cycles", 500 * domain_size**2)
    else:
            raise ValueError(key + " = " + simulation_type)

    print("Substituted coupling cycles")
