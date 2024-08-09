from generators.utils import get_config_value


def apply(partialXml: object, configs: dict) -> None:
    value = "yes" if  get_config_value(configs, __name__.split(".")[-1]) else "no"
    partialXml.substitute("{two-way-coupling}", value)
    print("Substituted 2-way coupling")
