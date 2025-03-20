#! /usr/bin/env python3

"""
This script parses the current remote and local state of git and the HTML output of the test coverage
to suggest a list of files for which more tests need to be written.
"""

__author__ = "Ruben Horn"

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup, Tag


def iterate_coverage_indices(folder):
    if not folder.exists():
        return []
    index_files = []
    for f in folder.iterdir():
        if f.is_dir():
            index_files += iterate_coverage_indices(f)
        elif f.name.lower() == "index.html":
            index_files.append(f)
    return index_files


def get_test_coverage(index_files):
    data = {}

    def parse_file(file):
        source_folder = str(file.parent)[len(str(coverage_root)) + 1 :]
        soup = BeautifulSoup(file.read_text(), "html.parser")
        root_element = soup.body.center.table

        def parse_coverage_table_row(row, out_data):
            cells = [td.text.strip() for td in row]
            values = []
            for j, s in enumerate(s for s in cells if len(s) > 0):
                if 0 == j:
                    values.append(f"{source_folder}/{s}")
                elif "%" in s:
                    values.append(float(s[:-2]) / 100)
            if len(out_data.keys()) - 1 == len(values):
                values.append(np.nan)  # Has no function coverage
            elif len(out_data.keys()) != len(values):
                print(out_data.keys(), values)
                return False
            for key, value in zip(out_data.keys(), values):
                out_data[key].append(value)
            return True

        rows = list(e for e in root_element if isinstance(e, Tag))
        # First row is empty -> skip
        # Second row contains header (filename and coverage)
        # Third row contains more headers (type of coverage)
        header_row = rows[1]
        keys = [td.text.strip() for td in header_row]
        keys = [k for k in keys if len(k) > 0]
        if keys[0] != "Filename":
            print(f"Skipping {file}", file=sys.stderr)
            return
        if len(data) == 0:
            for key in keys:
                data[key] = []
        # Parse coverage for each file in this folder
        for i, row in enumerate(rows[3:]):
            if not parse_coverage_table_row(row, data):
                print(f"Error parsing row {4 + i} in {file}")

    for f in index_files:
        parse_file(f)

    return pd.DataFrame(data)


if __name__ == "__main__":
    base_dir = Path(__file__).parent.parent
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-d", "--build-directory", type=Path, default=base_dir / Path("build")
    )
    arg_parser.add_argument("-M", "--skip-make-coverage", action="store_true")
    args = arg_parser.parse_args()
    build_dir = args.build_directory.absolute()

    if not args.skip_make_coverage:
        os.chdir(build_dir)
        os.system("make coverage") # Calls this script again

    print("\n=== Analysing coverage of locally changed files ===\n", file=sys.stderr)
    os.chdir(base_dir)

    if not Path(".git").is_dir():
        print(f"{os.getcwd()} is not a git repository!", file=sys.stderr)
        exit(1)

    # Look up target branch for pull request in CI/CD pipeline or default to master when running locally
    target_branch = ""
    if "GITHUB_BASE_REF" in os.environ:
        target_branch = os.environ["GITHUB_BASE_REF"].strip()
    if len(target_branch) == 0:
        target_branch = "master"

    # Make sure to compare to the current target branch on the remote
    subprocess.call(
        ["git", "fetch"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    uncommited_files = [
        l.split(" ")[-1]
        for l in subprocess.check_output(
            ["git", "status", "--porcelain"],
            stderr=subprocess.DEVNULL,
        )
        .decode()
        .strip()
        .splitlines()
    ]
    current_branch = (
        subprocess.check_output(["git", "branch", "--show-current"]).decode().strip()
    )
    # Fall back on branch name from environment when running on a pull request in CI/CD pipeline
    if len(current_branch) == 0 and "GITHUB_HEAD_REF" in os.environ:
        current_branch = os.environ["GITHUB_HEAD_REF"]
    if len(current_branch) == 0:
        print("Could not determine current branch!", file=sys.stderr)
        exit(1)
    changed_files = []
    remote = "origin"
    for branch_prefix in ["", f"{remote}/"]:
        try:
            changed_files += (
                subprocess.check_output(
                    [
                        "git",
                        "diff",
                        "--name-only",
                        f"{remote}/{target_branch}...{branch_prefix}{current_branch}",
                    ],
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .strip()
                .splitlines()
            )
        except:
            pass  # No such branch
    touched_files = set(uncommited_files + changed_files)

    coverage_root = build_dir / "coverage"

    if not coverage_root.is_dir():
        print(
            f"No test coverage found in {str(coverage_root)}!", file=sys.stderr, end=""
        )
        coverage_root = Path("docs") / "html" / "coverage"  # Check alternative location
        print(f" (Switching to {str(coverage_root)})", file=sys.stderr)

    if not coverage_root.is_dir():
        print("Coverage has not been generated yet.", file=sys.stderr)
        exit(1)

    coverage_index_files = iterate_coverage_indices(coverage_root)
    df_coverage = get_test_coverage(coverage_index_files)

    # Remove completely tested (100% line coverage)
    mask_not_completed = ~(df_coverage[df_coverage.columns[1:]] >= 1).apply(all, axis=1)
    df_coverage = df_coverage[mask_not_completed]
    # Only keep those that were changed
    mask_touched = df_coverage.apply(lambda r: r.values[0] in touched_files, axis=1)
    df_coverage = df_coverage[mask_touched]

    if len(df_coverage) == 0:
        print("Sufficient test coverage for changed files.")
    else:
        print("You should write tests for:")
        print(
            df_coverage.to_string(
                index=False,
                formatters={c: "{:,.2%}".format for c in df_coverage.columns[1:]},
            )
        )
