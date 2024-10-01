#! /usr/bin/env python3

"""
Utility script to generate basic couette.xml configurations for MaMiCo.
Do NOT run this script directly! Use ../run instead.
"""

import argparse
import importlib
import json
import os
import sys
from pathlib import Path

import term
from utils import check_if_replacing, get_asset_text
from xml_templating import PartialXml


def select_option(configs: list, key: str, value: str) -> None:
    """Select an option based on a serialized key value pair without user input.

    Keyword arguments:
    configs -- All configurations in which to update the configuration with the specified key
    key -- The key by which to select the configuration
    value -- The serialized representation of the value which should be selected
    """
    for config in configs:
        if config["key"] == key:
            get_selected(config)["selected"] = False
            for option in config["options"]:
                if str(option["value"]) == value:
                    option["selected"] = True
                    return
            raise ValueError(
                f'No option with the value "{value}" exists for the configuration with the key "{key}"'
            )
    raise ValueError(f'No configuration with the key "{key}" exists')


def select_option_interactive(configs: list, label: str) -> None:
    """Show a terminal menu for a single option of the configuration and apply the selection of the user to the current configurations.

    Keyword arguments:
    configs -- All configurations in which to update the configuration with the specified label
    label -- The menu label by which to select the configuration
    """
    selected_config = None
    for config in configs:
        if config["label"] == label:
            if selected_config is not None:
                raise RuntimeError(
                    f'Ambiguous selection. Multiple configs match the label "{label}"'
                )
            else:
                selected_config = config
    del config  # Avoid bugs by re-using variable
    pre_selected = None
    option_labels = []
    for i, option in enumerate(selected_config["options"]):
        if "selected" in option and option["selected"]:
            pre_selected = i
            option["selected"] = False
        label = option["label"]
        if "description" in option:
            label += "\t" + option["description"]
        option_labels.append(label)
    selected = term.select(
        option_labels, title=selected_config["label"], pre_selected=pre_selected
    )
    selected_config["options"][selected]["selected"] = True


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


def generate(configs: list, output_dir: str, replace_existing: bool) -> None:
    """Create the configuration by loading the template, applying the generators, and saving the resulting XML file.
    This function terminates the process.

    Keyword arguments:
    configs -- The list of configurations which to apply using the generators corresponding to the keys
    output_dir -- The directory where couette.xml and other files will be created
    replace_existing -- If False and an output file already exists, an error is thrown
    """
    xml = PartialXml(get_asset_text("couette.xml.template"))
    print("Loaded configuration template")
    # Make output path accessible to generators and consider existing files
    config_output = [
        # Pre-selected configs with only one option
        dict(key="output_dir", options=[dict(selected=True, value=output_dir)]),
        dict(
            key="replace_existing",
            options=[dict(selected=True, value=replace_existing)],
        ),
    ]
    _get_config_value = lambda k: get_config_value(configs + config_output, k)
    for config in configs:
        generator = load_generator(config["key"])
        generator.apply(xml, _get_config_value)
    path = output_dir / "couette.xml"
    check_if_replacing(path, _get_config_value)
    with open(path, "w") as file:
        file.write(xml.get())
    print(f"\nWrote configuration to {path}")
    ranks = get_config_value(configs, "mpi_ranks") * get_config_value(
        configs, "multi_md"
    )
    print(f'\nRun simulation using "mpirun -n {ranks} $MAMICO_DIR/build/couette"')
    exit(0)


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
                f"Could not invoke validation for generator \"{config['key']}\"",
                file=sys.stderr,
            )
            continue
        if validation_error is not None:
            validation_errors += validation_error + "\n"
    return validation_errors.strip()


def parse_args(argv: dict = sys.argv[1:], configs: list = []) -> object:
    """Parse the arguments (from the command line) and return them as a dictionary.
    If the usage should be printed, the program is terminated instead.

    Keyword arguments:
    argv -- The command line arguments excluding the name of the script itself or the interpreter
    configs -- All configurations which may be overwritten using command line arguments
    """
    # Extract all possible key=value overrides for help text and validation
    all_override_kvs = set()
    all_overrides = ""
    for config in configs:
        key = config["key"]
        all_overrides += f"\n\t{key}\t=\t<"
        for i, option in enumerate(config["options"]):
            value = str(option["value"])
            all_override_kvs.add(f"{key}={value}")
            if i > 0:
                all_overrides += "|"
            all_overrides += value
        all_overrides += ">"
    # Parse arguments
    arg_parser = argparse.ArgumentParser(
        prog=Path(__file__).parent.parent.name + "/run",
        description="""A simple interactive command line utility to generate basic couette.xml configurations for MaMiCo.

If all options are provided through the command line, the script is executed non-interactively.
(This may be useful for generating configurations in an automated testing environment.)""",
        epilog="Available overrides:" + all_overrides,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    arg_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="The directory where couette.xml and other files will be created",
    )
    arg_parser.add_argument(
        "-r",
        "--replace-existing",
        action="store_true",
        help="Replace couette.xml and other files if they exist",
    )
    arg_parser.add_argument(
        "-O",
        "--override",
        type=str,
        default="",
        help="Comma separated list of key-value (key=value) config overrides. (Will be removed from interactive menu)",
    )
    # Parse key=value overrides
    args = arg_parser.parse_args(argv)
    override = {}
    for kv in args.override.split(","):
        if len(kv.strip()) == 0:
            continue
        parts = kv.split("=")
        if len(parts) != 2:
            arg_parser.error(
                "The value of --override must be a comma separated list of key-value pairs (key=value)."
            )
        k, v = (p.strip() for p in parts)
        override[k] = v
    args.override = override
    # Check if every key=value override is valid
    for key, value in args.override.items():
        if f"{key}={value}" not in all_override_kvs:
            arg_parser.error(f'Invalid override "{key}={value}"')
    if not Path(args.output).is_dir():
        arg_parser.error(
            "The value of --output must be a path to a directory, not a file."
        )
    return args


def main() -> None:
    os.chdir(Path(__file__).parent)
    # 1. Load configuration options
    configs = json.loads(get_asset_text("configuration_template.json"))
    args = parse_args(configs=configs)
    # Fill missing labels
    for config in configs:
        if "label" not in config:
            config["label"] = str(config["key"])
        for option in config["options"]:
            if "label" not in option:
                option["label"] = str(option["value"])

    while True:
        # 2. Create the application menu
        title = term.fmt_bold("MaMiCo couette.xml Generator\n")
        main_menu = []
        for config in configs:
            # Apply overrides from command line
            key = config["key"]
            if key in args.override:
                select_option(configs, key, args.override[key])
            else:
                main_menu.append(config["label"] + "\t" + get_selected(config)["label"])
        # 3. Validate
        validation_errors = validate(configs)
        is_valid = len(validation_errors.strip()) == 0
        is_interactive = (
            len(main_menu) > 0
        )  # Skip interactive mode if there are no options
        if not is_interactive:
            if is_valid:
                generate(configs, args.output)
            else:
                print(validation_errors, file=sys.stderr)
                exit(1)
        else:
            if is_valid:
                main_menu[-1] += "\n"
                main_menu.append("Generate config file(s)")
            else:
                validation_errors = "\n".join(
                    [term.fmt_red(s) for s in validation_errors.splitlines()]
                )
                title += (
                    term.fmt_red(term.fmt_bold("\nValidation errors:\n"))
                    + validation_errors
                    + "\n"
                )
            # 4. Show and the menu
            title += "\n(Press Ctrl+C to quit without generating anything)\n"
            selected = term.select(
                main_menu,
                title=title,
                pre_selected=-1 if is_valid else None,
            )
            if selected == -1:
                exit(1)
            elif is_valid and selected == len(main_menu) - 1:
                # Generate
                generate(configs, args.output, args.replace_existing)
            else:
                # Update configuration
                # Remove explanation behind label
                selected_label = main_menu[selected].split("\t")[0]
                select_option_interactive(configs, selected_label)


if __name__ == "__main__":
    main()
