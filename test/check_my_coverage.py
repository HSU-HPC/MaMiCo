#! /usr/bin/env python3

'''
This script parses the current remote and local state of git and the HTML output of the test coverage
to suggest a list of files that still need to be tested.
'''

__author__ = 'Ruben'

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup, Tag

os.chdir(Path(__file__).parent.parent)

# Look up target branch for pull request in CI/CD pipeline or default to master when running locally
target_branch = os.environ['GITHUB_BASE_REF'] if 'GITHUB_BASE_REF' in os.environ else 'master'


# Make sure to compare to the current target branch on the remote
subprocess.call(['git', 'fetch'], stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
uncommited_files = [l.split(' ')[-1] for l in subprocess.check_output(
    ['git', 'status', '--porcelain']).decode().strip().splitlines()]
current_branch = subprocess.check_output(
    ['git', 'branch', '--show-current']).decode().strip()
changed_files = subprocess.check_output(
    ['git', 'diff', '--name-only', f'origin/{target_branch}...{current_branch}']).decode().strip().splitlines()
touched_files = set(uncommited_files + changed_files)

coverage_root = Path('build') / 'coverage'

if not coverage_root.is_dir():
    print('Coverage has not been generated yet.', file=sys.stderr)
    exit(1)


def iterate_coverage_indices(folder):
    index_files = []
    for f in folder.iterdir():
        if f.is_dir():
            index_files += iterate_coverage_indices(f)
        elif f.name.lower() == 'index.html':
            index_files.append(f)
    return index_files


coverage_index_files = iterate_coverage_indices(coverage_root)


def get_test_coverage(index_files):
    data = {}

    def parse_file(file):
        soup = BeautifulSoup(f.read_text(), 'html.parser')
        root_element = soup.body.center.table
        for i, tr in enumerate(e for e in root_element if isinstance(e, Tag)):
            if 0 == i:
                continue  # First row is empty -> skip
            elif 1 == i:
                keys = [td.text.strip() for td in tr]
                keys = [k for k in keys if len(k) > 0]
                if keys[0] != 'Filename':
                    print(f'Skipping {f}', file=sys.stderr)
                    return
                if len(data) == 0:
                    # Second row contains header (filename and coverage)
                    for td in tr:
                        key = td.text.strip()
                        if len(key) == 0:
                            continue
                        data[key] = []
            else:
                folder = str(f.parent)[len(str(coverage_root))+1:]
                cells = [td.text.strip() for td in tr]
                values = []
                for j, s in enumerate(s for s in cells if len(s) > 0):
                    if 0 == j:
                        values.append(f'{folder}/{s}')
                    elif '%' in s:
                        values.append(float(s[:-2]) / 100)

                # assert (len(data.keys()) == len(values))
                if len(data.keys()) != len(values):
                    print(f'Error in {f}')
                    continue

                for key, value in zip(data.keys(), values):
                    data[key].append(value)

    for f in index_files:
        parse_file(f)

    return pd.DataFrame(data)


df_coverage = get_test_coverage(coverage_index_files)

# Remove completely tested
mask_not_completed = ~(
    df_coverage[df_coverage.columns[1:]] >= 1).apply(all, axis=1)
df_coverage = df_coverage[mask_not_completed]
# Only keep those that were changed
mask_touched = df_coverage.apply(
    lambda r: r.values[0] in touched_files, axis=1)
df_coverage = df_coverage[mask_touched]

if len(df_coverage) == 0:
    print('Sufficient test coverage for changed files.')
    exit(0)
else:
    print('You should write tests for:')
    print(df_coverage.to_string(index=False, formatters={
        c: '{:,.2%}'.format for c in df_coverage.columns[1:]}))
    exit(1)
