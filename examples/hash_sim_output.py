#! /usr/bin/env python3

"""
Creates a list of hashes and output file suffixes ("..._<timestep>.<extension>").
Why? Allows for easy comparison of the simulation output between different implementations
(e.g. verify correctness of simulation after non-trivial refactoring)
"""

__author__ = "Ruben"

import os
import subprocess
import sys
from pathlib import Path

OUTPUT_FILE_EXTENSIONS = [
    ".csv",
    ".vtk"
]


if any(s in sys.argv for s in ["-h", "--help"]):
    the_script = Path(sys.argv[0]).name
    print(f"Usage: {the_script} path/to/sim/output/dir", file=sys.stderr)
    exit(1)

os.chdir(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent)

output_files = []
for file_extension in OUTPUT_FILE_EXTENSIONS:
    output_files += list(Path(".").glob("*" + file_extension))

for output_path in output_files:
    output_filename = str(output_path)
    hash = subprocess.check_output(
        ["md5sum", output_path], text=True).split(" ")[0]
    print(hash, "\t...", output_filename[output_filename.rindex("_"):], sep="")
