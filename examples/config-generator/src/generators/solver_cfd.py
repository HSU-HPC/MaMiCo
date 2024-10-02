import shutil
from pathlib import Path

from utils import check_if_replacing, get_domain_size


def _create_foam_setup(get_config_value) -> Path:
    src_path = Path(__file__).parent.parent.parent / "assets" / f"FoamSetup.template"
    dst_path = Path(get_config_value("output_dir")) / "FoamSetup"
    check_if_replacing(dst_path, get_config_value)
    shutil.copytree(src_path, dst_path)
    return dst_path


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    The analytical CFD solver does not support two-way coupling.
    Only the OpenFOAM solver should be used with two-way coupling.
    """
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    use_2way_coupling = get_config_value("coupling_2way")
    if solver == "foam" and get_domain_size(get_config_value) != 2:
        return f"OpenFOAM can only be used with the medium domain size."
    if solver == "analytical" and use_2way_coupling:
        return f"Cannot use two-way coupling with analytical CFD solver."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("cfd-solver", solver)
    print("Substituted CFD solver")
    foam_setup = ""
    if solver == "foam":
        path = _create_foam_setup(get_config_value)
        print(f"Also created files under {path}")
        output_dir = get_config_value("output_dir")
        foam_setup = (
            f'foam-setup-directory="{output_dir}"\n' + 'foam-setup-folder="FoamSetup"'
        )
    partial_xml.substitute("foam-setup", foam_setup)
