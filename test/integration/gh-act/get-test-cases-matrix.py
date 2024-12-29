#! /usr/bin/env python3

"""
Parses the configurations in test-cases.txt (for the tool mamico-cfg).
Builds the integration tests job matrix (for GitHub Action as JSON).
"""

import json
from pathlib import Path


def get_test_cases(filename):
    cells = []
    with open(filename, "r") as file:
        lines = [
            s.replace("\t", " ").replace("  ", " ").strip() for s in file.readlines()
        ]
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
