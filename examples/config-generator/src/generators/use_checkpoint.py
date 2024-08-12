def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    use_checkpoint = get_config_value(key)
    checkpoint_xml = '<checkpoint-configuration filename="CheckpointSimpleMD" write-every-timestep="0"></checkpoint-configuration>'
    partial_xml.substitute("checkpoint", checkpoint_xml if use_checkpoint else "")
    print("Substituted loading checkpoint")
    if use_checkpoint:
        print(
            'TODO: Please make sure to provide the file "CheckpointSimpleMD.checkpoint" in the same directory as couette.xml!'
        )
