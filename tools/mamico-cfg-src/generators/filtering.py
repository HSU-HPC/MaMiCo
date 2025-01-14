def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    filtering = get_config_value(key)

    xml_write_to_file = '<write-to-file location="per-instance-filtering-result.csv" overwrite="false" />'

    xml_2d_gaussian = """
<per-instance output="md">
    <gaussXY filtered-values="macro-mass macro-momentum">
        <gauss dim="0" sigma="1" extrapolation="mirror" />
        <gauss dim="1" sigma="1" extrapolation="mirror" />
    </gaussXY>
    <output filtered-values="macro-mass macro-momentum" input="gaussXY">
        {write-to-file}
    </output>
</per-instance>
""".strip()

    xml_nlm = """
<per-instance output="nlm-junction">
    <prefilter filtered-values="macro-mass macro-momentum">
        <gauss dim="0" sigma="1" extrapolation="mirror" />
        <gauss dim="1" sigma="1" extrapolation="mirror" />
        <gauss dim="2" sigma="1" extrapolation="mirror" />
        <constant dir="1" value="0" />
        <constant dir="2" value="0" />
    </prefilter>
    <nlm-junction filtered-values="macro-mass macro-momentum" input="md prefilter">
        <NLM time-window-size="5" sigsq="10" hsq="20" sigsq_rel="0.05" hsq_rel="0.1" />
    </nlm-junction>
    <output filtered-values="macro-mass macro-momentum" input="nlm-junction">
        {write-to-file}
    </output>
</per-instance>
""".strip()

    if filtering == False:
        per_instance_filtering_xml = '<per-instance output="md"></per-instance>'
    elif filtering == "gauss":
        per_instance_filtering_xml = xml_2d_gaussian
    elif filtering == "nlm":
        per_instance_filtering_xml = xml_nlm
    else:
        raise NotImplementedError(key + " = " + filtering)

    partial_xml.substitute("per-instance-filtering", per_instance_filtering_xml)
    # Output the filter result to a file so it is not lost
    should_write_to_file = filtering != False and not get_config_value("coupling_2way")
    partial_xml.substitute(
        "write-to-file",
        xml_write_to_file
        if should_write_to_file
        else "<!-- No output of filtering result -->",
    )
    print("Substituted filtering pipeline")
