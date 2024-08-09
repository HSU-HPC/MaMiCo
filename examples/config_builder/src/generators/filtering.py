from generators.utils import get_config_value


def apply(partialXml: object, configs: dict) -> None:
    key = __name__.split(".")[-1]
    filtering = get_config_value(configs, key)

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

    partialXml.substitute("{filtering}", filter_xml)
    print("Substituted filtering pipeline")
