from pathlib import Path

from generators.cell_size import get_molecules_per_direction
from generators.domain_size import get_domain_sizes
from utils import MD_DENSITY, check_if_replacing, get_asset_text
from xml_templating import PartialXml


def get_wall_velocity(get_config_value) -> float:
    return 1.5 if get_config_value("solver_md") == "ls1" else 0.5


def _create_ls1_config(get_config_value) -> Path:
    xml = PartialXml(get_asset_text("ls1config.xml.template"))
    md_domain_size, _ = get_domain_sizes(get_config_value)
    xml.substitute("md-size", md_domain_size)
    molecules_per_direction = get_molecules_per_direction(get_config_value)
    grid_filler_lattice = md_domain_size / molecules_per_direction
    xml.substitute("grid-filler-lattice", grid_filler_lattice)
    boundary_condition = get_config_value("boundary")
    xml.substitute("boundary-condition", boundary_condition)
    path = Path(get_config_value("output_dir")) / "ls1config.xml"
    check_if_replacing(path, get_config_value)
    with open(path, "w") as file:
        file.write(xml.get())
    return path


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("md-solver", solver)
    partial_xml.substitute("wall-velocity", get_wall_velocity(get_config_value))
    partial_xml.substitute("density", MD_DENSITY)
    parallel_topology_type = "zyx" if solver == "ls1" else "xyz"
    partial_xml.substitute("parallel-topology-type", parallel_topology_type)
    print("Substituted MD solver and wall velocity")
    if solver == "ls1":
        path = _create_ls1_config(get_config_value)
        print(f"Also created file {path}")
