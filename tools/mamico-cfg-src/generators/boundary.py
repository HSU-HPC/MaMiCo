def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    boundary_condition = get_config_value(key)
    boundary_condition_xml = None
    boundary_condition_value = "periodic" if boundary_condition == "periodic" else "reflective"
    boundary_condition_XY_value = "periodic" if "periodic" in boundary_condition else "reflective"
    if boundary_condition == "reflecting":
        boundary_condition_xml = """
<particle-insertion type="usher" maximum-number-of-iterations="100" maximum-number-of-restarts="500" insert-every-timestep="10" tolerance="0.5" />
    <boundary-force type="zhou-boundary-force" west="yes" east="yes" north="yes" south="yes" bottom="yes" top="yes" density="0.81" temperature="1.1" />
""".strip()
    elif boundary_condition == "periodic":
        boundary_condition_xml = """
<particle-insertion type="none" />
    <boundary-force type="none" west="no" east="no" north="no" south="no" bottom="no" top="no" density="0.81" temperature="1.1"/>
""".strip()
    elif boundary_condition == "XY-periodic":
        boundary_condition_xml = """
<particle-insertion type="none" />
    <boundary-force type="zhou-boundary-force" west="no" east="no" north="no" south="no" bottom="yes" top="yes" density="0.81" temperature="1.1" />
""".strip()
    else:
        raise ValueError(f'Invalid boundary condition "boundary_condition"')
    partial_xml.substitute("boundary-condition", boundary_condition_value)
    partial_xml.substitute("boundary-condition-XY", boundary_condition_XY_value)
    partial_xml.substitute("boundary-condition-xml", boundary_condition_xml)
    print("Substituted boundary condition")
