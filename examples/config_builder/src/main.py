#! /usr/bin/env python3

import importlib
import json
import os
from pathlib import Path
import sys

from generators.utils import get_config_value, get_selected
import term

import lxml.etree as etree


def select_option(config: dict) -> None:
    print(config)
    pre_selected = None
    option_labels = []
    for i, option in enumerate(config["options"]):
        if "selected" in option and option["selected"]:
            pre_selected = i
            option["selected"] = False
        label = option["label"]
        if "description" in option:
            label += "\t" + option["description"]
        option_labels.append(label)
    selected = term.select(
        option_labels, title=config["label"], pre_selected=pre_selected
    )
    config["options"][selected]["selected"] = True


def get_text_file(name: str) -> str:
    base_path = Path().absolute().parent
    return (base_path / "assets" / name).read_text()


def generate(configs):
    class PartialXml:
        def __init__(self) -> None:
            self.string = get_text_file("couette.xml.template")
            print("Loaded configuration template")

        def substitute(self, a: str, b: object) -> None:
            self.string = self.string.replace(a, str(b))

    xml = PartialXml()
    for config in configs:
        generator = config["key"]
        try:
            generator = importlib.import_module(f"generators.{generator}")
        except ModuleNotFoundError:
            print(f'Skipping missing generator "{generator}"')
            continue
        generator.apply(xml, configs)
    if "{" in xml.string:
        print("\nWarning: Not all values have been substituted!")
    output_filename = sys.argv[1]
    with open(output_filename, "w") as file:
        file.write(xml.string)
    # Format
    parsed = etree.parse(output_filename)
    with open(output_filename, "wb") as file:
        file.write(etree.tostring(parsed, pretty_print=True))
    print(f"\nWrote configuration to {output_filename}")
    ranks = get_config_value(configs, "mpi_ranks") * get_config_value(configs, "multi_md")
    print(f"\nRun simulation using \"mpirun -n {ranks} $MAMICO_PATH/build/couette\"")


def main() -> None:
    os.chdir(Path(__file__).parent)
    title = "MaMiCo couette.xml Generator\n"
    # 1. Load configuration options
    configs = json.loads(get_text_file("configuration_template.json"))
    # Fill missing labels
    for config in configs:
        if "label" not in config:
            config["label"] = str(config["key"])
        for option in config["options"]:
            if "label" not in option:
                option["label"] = str(option["value"])

    while True:
        # 2. Create the application menu
        main_menu = []
        for config in configs:
            main_menu.append(config["label"] + "\t" + get_selected(config)["label"])
        # 3. Validate
        is_valid = True  # TODO
        if is_valid:
            main_menu[-1] += "\n"
            main_menu.append("Generate")
        # 4. Show and the menu
        selected = term.select(
            main_menu, title=title, pre_selected=-1 if is_valid else None
        )
        if selected == -1:
            exit(1)
        elif is_valid and selected == len(main_menu) - 1:
            # Generate
            generate(configs)
            exit()
        else:
            # Update configuration
            select_option(configs[selected])


if __name__ == "__main__":
    main()
