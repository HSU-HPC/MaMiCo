def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    # ls1 has always run thermostatted everywhere; keep that.
    thermostat = "all" if get_config_value("solver_md") == "ls1" else get_config_value(key)
    partial_xml.substitute("thermostat-type", thermostat)
    print("Substituted thermostat region")
