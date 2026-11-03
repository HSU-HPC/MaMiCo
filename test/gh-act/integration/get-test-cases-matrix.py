#! /usr/bin/env python3

"""
Parses the configurations in test-cases.txt (for the tool mamico-cfg).
Builds the integration tests job matrix (for GitHub Action as JSON).
"""

import json
from pathlib import Path


def clean_up_line(line):
    line = line.split("#")[0]  # Remove comment
    line = line.replace("\t", " ")
    line = line.strip()
    return line


def get_test_cases(filename):
    cells = []
    with open(filename, "r") as file:
        lines = [clean_up_line(s) for s in file.readlines()]
        lines = [s for s in lines if len(s) > 0]  # Remove empty lines
        cells = [s.split() for s in lines]
    for row in cells[1:]:
        yield {k: v for k, v in zip(cells[0], row)}


def get_gh_action_matrix():
    print("matrix=", end="")
    test_cases_path = Path(__file__).parent / "test-cases.txt"
    matrix = {"include": list(get_test_cases(test_cases_path))}
    print(json.dumps(matrix))


if __name__ == "__main__":
    get_gh_action_matrix()
