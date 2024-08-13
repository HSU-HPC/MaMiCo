from pathlib import Path

from generators.cell_size import get_molecules_per_direction
from generators.domain import get_domain_sizes
from utils import get_asset_text
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
    filename = Path(get_config_value("output_filename")).parent / "ls1config.xml"
    with open(filename, "w") as file:
        file.write(xml.get())


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("md-solver", solver)
    wall_velocity = 0.5 if solver == "md" else 1.5
    partial_xml.substitute("wall-velocity", wall_velocity)
    print("Substituted MD solver and wall velocity")
    if solver == "ls1":
        _create_ls1_config(get_config_value)
        print('Also created file "ls1config.xml" in the same directory as couette.xml.')
