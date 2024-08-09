from generators.utils import get_config_value


def apply(partialXml: object, configs: dict) -> None:
    key = __name__.split(".")[-1]
    simulation_type = get_config_value(configs, key)
    domain_size = get_config_value(configs, "domain")
    match simulation_type:
        case "test":
            partialXml.substitute("{coupling-cycles}", 50)
        case "small":
            partialXml.substitute("{coupling-cycles}", 250 * domain_size**2)
        case "normal":
            partialXml.substitute("{coupling-cycles}", 500 * domain_size**2)
        case _:
            raise ValueError(key + " = " + simulation_type)
    
    print("Substituted coupling cycles")
