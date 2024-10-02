"""General utilities used in multiple modules."""

from pathlib import Path

import term


def get_domain_size(get_config_value) -> int:
    """Return the domain size as an ordinal integer where 1=small, 2=medium, and 3=large.

    Keyword arguments:
    get_config_value -- A function to get the configuration values by key
    """
    domain_size_name = get_config_value("domain_size")
    domain_size_names = ["small", "medium", "large"]
    return domain_size_names.index(domain_size_name) + 1


def get_cell_size(get_config_value) -> float:
    """Return the cell size.

    Keyword arguments:
    get_config_value -- A function to get the configuration values by key
    """
    return get_config_value("cell_size")


def check_if_replacing(path: Path, get_config_value) -> None:
    """Exit with an error if a would be replaced and this is not requested explicitly.

    Keyword arguments:
    path -- The path to the file that will be written to
    get_config_value -- A function to get the configuration values by key
    """
    if not get_config_value("replace_existing") and path.exists():
        print(
            term.fmt_red(
                f"The file {path} already exists.\nUse --replace-existing to overwrite it!"
            )
        )
        exit(1)


def get_asset_text(name: str) -> str:
    """Load the text from a file in the assets folder next to the src folder.

    Keyword arguments:
    name -- The filename of the asset in the src folder including extension
    """
    base_path = Path().absolute().parent
    return (base_path / "assets" / name).read_text()
