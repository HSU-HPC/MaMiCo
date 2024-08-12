"""Example for generators used to create the couette.xml file from the template.
This generator is applied for '"key": "example"' in the configuration dictionary.
In the file configuration_template.json, this is equivalent to:
```json
[
  ...
  {
    "key": "example",
    "options": [
      { "value": "world", "selected": true },
      { "value": "space" }
    ]
  }
  ...
]
```
"""


def apply(partial_xml, get_config_value) -> None:
    """Applies all substitutions for this generator to the current template.

    Keyword arguments:
    partial_xml -- The PartialXml object containing the current template to which the substitutions are applied
    get_config_value -- A function to get the configuration values by key
    """
    # Load the value of the generator from the configuration (Generators may load any other value for validation/substitution)
    key = __name__.split(".")[-1]
    value = get_config_value(key)
    placeholder = "examle-substitution"
    # Perform the substitution (The placeholder may be different from the generator name)
    # Here, "{example-substitution}" is replaced by "Hello world!" or " Hello space!"
    partial_xml.substitute(placeholder, f"Hello {value}!")
