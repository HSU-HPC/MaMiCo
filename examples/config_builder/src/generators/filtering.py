def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    filtering = get_config_value(key)

    xml_2d_gaussian = """
<per-instance output = "gaussXY" >
    <gaussXY filtered-values="macro-mass macro-momentum">
        <gauss dim="0" sigma="1" extrapolation="mirror" />
        <gauss dim="1" sigma="1" extrapolation="mirror" />
    </gaussXY>
</per-instance>
""".strip()

    if filtering == False:
        filter_xml = ""
    elif filtering == "gauss":
        filter_xml = xml_2d_gaussian
    else:
        raise NotImplementedError(key + " = " + filtering)

    partial_xml.substitute("{filtering}", filter_xml)
    print("Substituted filtering pipeline")
