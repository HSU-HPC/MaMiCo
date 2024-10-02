import shutil
from pathlib import Path

from generators.domain_size import get_domain_size
from utils import check_if_replacing


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    A pre-computed checkpoint is (currently) only available for Simple MD with a 30x30x30 domain.
    """
    key = __name__.split(".")[-1]
    use_checkpoint = get_config_value(key)
    domain_size = get_domain_size(get_config_value)
    solver_md = get_config_value("solver_md")
    if use_checkpoint and (domain_size > 1 or solver_md != "md"):
        return f"Example checkpoint file only provided for small Simple MD simulation."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    use_checkpoint = get_config_value(key)
    checkpoint_xml = '<checkpoint-configuration filename="CheckpointSimpleMD" write-every-timestep="0"></checkpoint-configuration>'
    partial_xml.substitute("checkpoint-out", checkpoint_xml if use_checkpoint else "")
    checkpoint_key = 'init-from-sequential-checkpoint="CheckpointSimpleMD"'
    partial_xml.substitute("checkpoint-in", checkpoint_key if use_checkpoint else "")
    print("Substituted loading checkpoint")
    if use_checkpoint:
        boundary_condition = get_config_value("boundary")
        checkpoint_src_path = (
            Path(__file__).parent.parent.parent
            / "assets"
            / f"CheckpointSimpleMD_10000_{boundary_condition}_0.checkpoint"
        )
        checkpoint_dst_path = (
            Path(get_config_value("output_dir")) / "CheckpointSimpleMD.checkpoint"
        )
        check_if_replacing(checkpoint_dst_path, get_config_value)
        # Avoid reading and writing the contents of the file, because it is rather large
        shutil.copyfile(checkpoint_src_path, checkpoint_dst_path)
        print(
            'Also created file "CheckpointSimpleMD.checkpoint" in the same directory as couette.xml.'
        )
