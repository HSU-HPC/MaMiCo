from pathlib import Path

from generators.cell_size import get_molecules_per_direction
from generators.domain_size import get_domain_sizes
from utils import check_if_replacing, get_asset_text
from xml_templating import PartialXml


def _create_ls1_config(get_config_value) -> None:
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


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("md-solver", solver)
    wall_velocity = 1.5 if solver == "ls1" else 0.5
    partial_xml.substitute("wall-velocity", wall_velocity)
    thermostat_type = "all" if solver == "ls1" else "outerLayers"
    partial_xml.substitute("thermostat-type", thermostat_type)
    parallel_topology_type = "zyx" if solver == "ls1" else "xyz"
    partial_xml.substitute("parallel-topology-type", parallel_topology_type)
    print("Substituted MD solver and wall velocity")
    if solver == "ls1":
        _create_ls1_config(get_config_value)
        print('Also created file "ls1config.xml" in the same directory as couette.xml.')
