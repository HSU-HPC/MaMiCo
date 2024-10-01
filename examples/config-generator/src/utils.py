"""General utilities used in multiple modules."""

from pathlib import Path

import term


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
