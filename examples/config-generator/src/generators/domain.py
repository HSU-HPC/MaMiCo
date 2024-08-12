def apply(partial_xml, get_config_value) -> None:
    md_base_size = 30
    cfd_base_size = 50
    domain_offset_xy = (cfd_base_size - md_base_size) / 2  # Centered
    equilibration_steps = 10000
    equilibration_steps_max = 20000
    key = __name__.split(".")[-1]
    size = get_config_value(key)
    partial_xml.substitute("md-size", md_base_size * size)
    partial_xml.substitute("cfd-size", cfd_base_size * size)
    partial_xml.substitute(
        "equilibration-steps",
        min(equilibration_steps_max, equilibration_steps * size),
    )
    partial_xml.substitute("domain-offset-x", domain_offset_xy * size)
    partial_xml.substitute("domain-offset-y", domain_offset_xy * size)
    cell_size = get_config_value(key)
    partial_xml.substitute("domain-offset-z", cell_size)
    print(
        "Substituted domain size and offset for MD, CFD solver domain size, and MD equilibration steps"
    )
