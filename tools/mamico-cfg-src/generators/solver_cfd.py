import shutil
from pathlib import Path

from utils import check_if_replacing, get_domain_size
import sys

def _create_foam_setup(get_config_value) -> Path:
    src_path = Path(__file__).parent.parent / "assets" / f"FoamSetup.template"
    dst_path = Path(get_config_value("output_dir")) / "FoamSetup"
    check_if_replacing(dst_path, get_config_value)
    if sys.version_info >= (3, 8):
        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
    else:
        try:
             shutil.copytree(src_path, dst_path)
        except FileExistsError:
             pass
    return dst_path


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    The analytical CFD solver does not support two-way coupling.
    Only the OpenFOAM solver should be used with two-way coupling.
    """
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    use_2way_coupling = get_config_value("coupling_2way")
    cell_size = get_config_value("cell_size")
    if solver == "foam" and get_domain_size(get_config_value) != 2:
        return f"OpenFOAM can only be used with the medium domain size."
    if solver == "foam" and cell_size != 5.0:
        return f"OpenFOAM can only be used with cell size 5.0"
    if solver == "analytical" and use_2way_coupling:
        return f"Cannot use two-way coupling with analytical CFD solver."
    if (
        solver == "lb"
        and get_config_value("simulation") != "test"
        and use_2way_coupling
    ):
        return f"Two-way coupling is not stable with LBCouette solver for more than a quick test."
    if (
        solver == "fd"
        and get_config_value("simulation") != "test"
        and use_2way_coupling
    ):
        return f"Two-way coupling is not stable with FDCouette solver for more than a quick test."
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
