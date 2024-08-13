# MaMiCo Config Generator
This simple utility can be used to generate basic `couette.xml` files for the Couette flow test scenario using MaMiCo.

Using [`assets/couette.xml.template`](./assets/couette.xml.template), this Python script substitutes all occurances of `{<variable-name>}` with concrete values.  
This is done by iterating the (meta) configuration (key-value pairs) and applying the corresponding [generators/](./src/generators/)`<key>` python module.  
A simple [example generator](./src/generators/example.py) is provided to explain how to implement/extend them.
The initial values are loaded from [`assets/configuration_template.json`](./assets/configuration_template.json).
This file has the following structure:
```json
[
  {
    "key": "the_generator_name",
    "label": "Optional label in menu (fallback to key)",
    "options": [
      {
        "selected": true,
        "value": 1,
        "label": "Optional label in menu (fallback to value)",
        "description": "Optional description of the option. Exactly one option must be pre-selected!"
      },
      {
        "value": 2,
        "label": "A second option"
      }
    ]
  },
  {
    "key": "another_generator",
    "options": [
      {
        "selected": true,
        "value": true
      }
    ]
  }
]
```

### Requirements
- Python 3.10+
- Bash

### Usage
1. Navigate to the target directory for the `couette.xml` file (most likely the MaMiCo build folder)
2. Execute `$MAMICO_DIR/examples/config-generator/run`  
Specific values can be fixed by providing them as a command line argument, e.g. `--override=domain_size=large,use_checkpoint=False,multi_md=2`.  
The corresponding entries will be removed from the interactive menu. If all configurations are overwritten, the script will run non-interactively. (May be used for automated testing.)
3. Customize the generator settings and generate the `couette.xml` file
4. Customize the generated configuration even further according to the [documentation](https://github.com/HSU-HPC/MaMiCo/wiki/couette.xml) (optional)
