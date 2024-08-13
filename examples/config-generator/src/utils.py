"""General utilities used in multiple modules."""

from pathlib import Path


def get_asset_text(name: str) -> str:
    """Load the text from a file in the assets folder next to the src folder."""
    base_path = Path().absolute().parent
    return (base_path / "assets" / name).read_text()
