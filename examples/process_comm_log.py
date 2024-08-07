#! /usr/bin/env python3

import os
from pathlib import Path

# Merge files matching this pattern:
input_files_pattern = "comm.log*"

os.chdir(Path(__file__).parent)

input_files = list(Path().glob(input_files_pattern))

# Filter out the output of this script if it exists
input_files = [file for file in input_files if not file.name.endswith(input_files_pattern[:-1])]

if len(input_files) == 0:
    print(f"Files {input_files_pattern} do not exist")
    exit(1)

input_files.sort(key=lambda p: int(p.name[len(input_files_pattern)-1:]))

file_lines = [file.read_text().splitlines() for file in input_files]
i = 0
with open(input_files_pattern[:-1], "w") as file:
    while True:
        if len(file_lines[0]) <= i:
            break
        for j in range(len(file_lines)):
            print(file_lines[j][i], file=file)
        i += 1

for file in input_files:
    file.unlink()
