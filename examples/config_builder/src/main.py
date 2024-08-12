#! /usr/bin/env python3

import importlib
import json
import os
import sys
from pathlib import Path

import term
from xml_templating import PartialXml


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


def get_selected(config: dict) -> object:
    for option in config["options"]:
        if "selected" in option and option["selected"]:
            return option
    raise ValueError(f'No option selected for config with key "{config["key"]}"')


def get_config_value(configs: dict, key: str) -> object:
    for config in configs:
        if config["key"] == key:
            return get_selected(config)["value"]
    return ValueError(f'No config with key "{key}" found')


def generate(configs):
    xml = PartialXml(get_text_file("couette.xml.template"))
    print("Loaded configuration template")
    for config in configs:
        generator = config["key"]
        try:
            generator = importlib.import_module(f"generators.{generator}")
        except ModuleNotFoundError:
            print(f'Missing generator "{generator}"', file=sys.stderr)
            exit(1)
        generator.apply(xml, lambda k: get_config_value(configs, k))
    output_filename = sys.argv[1]
    with open(output_filename, "w") as file:
        file.write(xml.get())
    print(f"\nWrote configuration to {output_filename}")


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
            ranks = get_config_value(configs, "mpi_ranks") * get_config_value(
                configs, "multi_md"
            )
            print(
                f'\nRun simulation using "mpirun -n {ranks} $MAMICO_PATH/build/couette"'
            )
            exit(0)
        else:
            # Update configuration
            select_option(configs[selected])


if __name__ == "__main__":
    main()
