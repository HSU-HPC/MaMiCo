def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    value = "yes" if get_config_value(key) else "no"
    partial_xml.substitute("two-way-coupling", value)
    print("Substituted 2-way coupling")
