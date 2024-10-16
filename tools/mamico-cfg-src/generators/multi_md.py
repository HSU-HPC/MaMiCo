def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    num_md_sims = get_config_value(key)
    partial_xml.substitute("num-md-sims", num_md_sims)
    print("Substituted number of MD simulations")
