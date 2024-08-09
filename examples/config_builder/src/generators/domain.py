from generators.utils import get_config_value

def apply(partialXml: object, configs: dict) -> None:
    md_base_size = 30
    cfd_base_size = 50
    equilibration_steps = 10000
    equilibration_steps_max = 20000
    size = get_config_value(configs, __name__.split(".")[-1])
    partialXml.substitute("{md-size}", md_base_size * size)
    partialXml.substitute("{cfd-size}", cfd_base_size * size)
    partialXml.substitute("{equilibration-steps}", min(equilibration_steps_max, equilibration_steps * size))
    print("Substituted domain size for MD and CFD solver and MD equilibration steps")
