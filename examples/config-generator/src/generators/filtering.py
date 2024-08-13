def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    filtering = get_config_value(key)

    xml_2d_gaussian = """
<per-instance output="md">
    <gaussXY filtered-values="macro-mass macro-momentum">
        <gauss dim="0" sigma="1" extrapolation="mirror" />
        <gauss dim="1" sigma="1" extrapolation="mirror" />
    </gaussXY>
</per-instance>
""".strip()

    xml_nlm = """
<per-instance output="nlm-junction">
    <prefilter filtered-values="macro-mass macro-momentum">
        <gauss dim="0" sigma="1" extrapolation="mirror" />
        <gauss dim="1" sigma="1" extrapolation="mirror" />
        <gauss dim="2" sigma="1" extrapolation="mirror" />
    </prefilter>
    <nlm-junction filtered-values="macro-mass macro-momentum" input="md prefilter">
        <NLM time-window-size="5" sigsq="10" hsq="20" sigsq_rel="0.05" hsq_rel="0.1" />
    </nlm-junction>
</per-instance>
""".strip()

    if filtering == False:
        per_instance_filtering_xml = "<per-instance output=\"md\"></per-instance>"
    elif filtering == "gauss":
        per_instance_filtering_xml = xml_2d_gaussian
    elif filtering == "nlm":
        per_instance_filtering_xml = xml_nlm
    else:
        raise NotImplementedError(key + " = " + filtering)

    partial_xml.substitute("per-instance-filtering", per_instance_filtering_xml)
    print("Substituted filtering pipeline")
