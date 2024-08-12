#! /usr/bin/env python3

"""
Utility script to generate basic couette.xml configuration for MaMiCo.
Do NOT run this script directly! Use ../run instead.
"""

import importlib
import json
import os
import sys
from pathlib import Path

import term
from xml_templating import PartialXml


def select_option(config: dict) -> None:
    """Show a terminal menu for a single option of the configuration and apply the selection of the user to the current configurations.

    Keyword arguments:
    config -- The configuration from all configurations (by reference) which should be updated
    """
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
    """Load the text from a file in the assets folder next to the src folder."""
    base_path = Path().absolute().parent
    return (base_path / "assets" / name).read_text()


def get_selected(config: dict) -> object:
    """Return the selected option for a given configuration. The absence of such an option causes a ValueError.

    Keyword arguments:
    config -- The configuration from which to look up the selected option
    """
    for option in config["options"]:
        if "selected" in option and option["selected"]:
            return option
    raise ValueError(f'No option selected for config with key "{config["key"]}"')


def get_config_value(configs: list, key: str) -> object:
    """Return the selected value for a given configuration by its key.

    Keyword arguments:
    configs -- All configurations from which to return the value of the configuration with the specified key
    key -- The key of the configuration from which to return the value
    """
    for config in configs:
        if config["key"] == key:
            return get_selected(config)["value"]
    return ValueError(f'No config with key "{key}" found')


def load_generator(generator: str) -> object:
    """Load a generator module by name. (Process exists with error if the generator does not exist.)

    Keyword arguments:
    generator -- the filename without extension of the module in the generators folder
    """
    try:
        return importlib.import_module(f"generators.{generator}")
    except ModuleNotFoundError:
        print(f'Missing generator "{generator}"', file=sys.stderr)
        exit(1)


def generate(configs: list, filename: str) -> None:
    """Create the configuration by loading the template, applying the generators, and saving the resulting XML file.

    Keyword arguments:
    configs -- The list of configurations which to apply using the generators corresponding to the keys
    filename -- The output XML filename
    """
    xml = PartialXml(get_text_file("couette.xml.template"))
    print("Loaded configuration template")
    for config in configs:
        generator = load_generator(config["key"])
        # Make output path accessible to generators
        config_output_filename = [
            dict(key="output_filename", options=[dict(selected=True, value=filename)])
        ]
        generator.apply(
            xml, lambda k: get_config_value(configs + config_output_filename, k)
        )
    with open(filename, "w") as file:
        file.write(xml.get())
    print(f"\nWrote configuration to {filename}")


def validate(configs: list) -> str:
    """validate the configuration using the generators, returning a string of validation erros (empty if valid).

    Keyword arguments:
    configs -- The list of configurations which to validate using the generators corresponding to the keys
    """
    validation_errors = ""
    for config in configs:
        generator = load_generator(config["key"])
        try:
            validation_error = generator.validate(
                lambda k: get_config_value(configs, k)
            )
        except AttributeError:
            print(
                f"FIXME: Could not invoke validation for generator \"{config['key']}\""
            )
            continue
        if validation_error is not None:
            validation_errors += validation_error + "\n"
    return validation_errors.strip()


def main() -> None:
    os.chdir(Path(__file__).parent)
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
        title = "MaMiCo couette.xml Generator\n"
        main_menu = []
        for config in configs:
            main_menu.append(config["label"] + "\t" + get_selected(config)["label"])
        # 3. Validate
        validation_errors = validate(configs)
        is_valid = len(validation_errors.strip()) == 0
        if is_valid:
            main_menu[-1] += "\n"
            main_menu.append("Generate")
        else:
            title += "\n" + validation_errors + "\n"
        # 4. Show and the menu
        selected = term.select(
            main_menu,
            title=title,
            pre_selected=-1 if is_valid else None,
        )
        if selected == -1:
            exit(1)
        elif is_valid and selected == len(main_menu) - 1:
            # Generate
            generate(configs, sys.argv[1])
            ranks = get_config_value(configs, "mpi_ranks") * get_config_value(
                configs, "multi_md"
            )
            print(
                f'\nRun simulation using "mpirun -n {ranks} $MAMICO_DIR/build/couette"'
            )
            exit(0)
        else:
            # Update configuration
            select_option(configs[selected])


if __name__ == "__main__":
    main()
