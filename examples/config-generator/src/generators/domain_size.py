from generators.mpi_ranks import get_ranks_xyz
from utils import get_cell_size, get_domain_size


def get_domain_sizes(get_config_value) -> int:
    size = get_domain_size(get_config_value)
    md_base_size = 30
    cfd_base_size = 50
    md_domain_size = md_base_size * size
    cfd_domain_size = cfd_base_size * size
    if size == 3:
        md_domain_size = 120  # instead of 90
        cfd_domain_size = 200  # instead of 150
    return md_domain_size, cfd_domain_size


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    For each axis x/y/z the MD domain must be divisible by the corresponding number of MPI ranks along the same axis.
    """
    md_domain_size, _ = get_domain_sizes(get_config_value)
    ranks_xyz = get_ranks_xyz(get_config_value)

    for ranks in ranks_xyz:
        if md_domain_size % ranks != 0:
            ranks_xyz_fmtd = "\u00d7".join([str(x) for x in ranks_xyz])
            return f"Cannot decompose MD domain of {md_domain_size}\u00b3 over {ranks_xyz_fmtd} ranks."
    return None


def apply(partial_xml, get_config_value) -> None:
    equilibration_steps = 10000
    equilibration_steps_max = 20000
    size = get_domain_size(get_config_value)
    md_domain_size, cfd_domain_size = get_domain_sizes(get_config_value)
    domain_offset_xy = (cfd_domain_size - md_domain_size) / 2  # Centered
    partial_xml.substitute("md-size", md_domain_size)
    partial_xml.substitute("cfd-size", cfd_domain_size)
    use_checkpoint = get_config_value("use_checkpoint")
    partial_xml.substitute(
        "equilibration-steps",
        min(equilibration_steps_max, equilibration_steps * size)
        * (0 if use_checkpoint else 1),
    )
    partial_xml.substitute("domain-offset-x", domain_offset_xy)
    partial_xml.substitute("domain-offset-y", domain_offset_xy)
    cell_size = get_cell_size(get_config_value)
    partial_xml.substitute("domain-offset-z", cell_size)
    print(
        "Substituted domain size and offset for MD, CFD solver domain size, and MD equilibration steps"
    )
