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
