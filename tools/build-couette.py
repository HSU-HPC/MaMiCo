#! /usr/bin/env python3

import argparse
import getpass
import os
import shutil
import socket
import subprocess
import time
from pathlib import Path

"""
Script to build the MaMiCo couette flow example using a single command.
(Must be located in the same git repository.)
"""


class ChangeDir:
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir

    def __enter__(self):
        self.last_dir = os.getcwd()
        os.chdir(self.tmp_dir)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.last_dir)


def get_git_repo_root(cwd="."):
    if not Path(cwd).is_dir():
        cwd = Path(cwd).parent
    with ChangeDir(cwd):
        repo_root = subprocess.getoutput("git rev-parse --show-toplevel")
        if repo_root.startswith("fatal:"):
            return None
        return Path(repo_root)


# MaMiCo
MAMICO_REPO_DIR = get_git_repo_root(Path(__file__))
MAMICO_BUILD_DIR = MAMICO_REPO_DIR / "build"
MAMICO_BUILD_TYPE = "DebugOptimized"

# LAMMPS
LAMMPS_REPO_DIR = MAMICO_BUILD_DIR / "LAMMPS"
LAMMPS_REPO_URL = "https://github.com/lammps/lammps.git"
LAMMPS_REPO_BRANCH = "release"

# OpenFOAM
OPEN_FOAM_BASE_DIR = MAMICO_BUILD_DIR / "openfoam"
OPEN_FOAM_VERSION = "2206"
OPEN_FOAM_SRC_DIR = OPEN_FOAM_BASE_DIR / f"OpenFOAM-v{OPEN_FOAM_VERSION}"
OPEN_FOAM_URL = f"https://dl.openfoam.com/source/v{OPEN_FOAM_VERSION}/OpenFOAM-v{OPEN_FOAM_VERSION}.tgz"
OPEN_FOAM_THIRD_PARTY_REPO_DIR = OPEN_FOAM_SRC_DIR / "ThirdParty"
OPEN_FOAM_THIRD_PARTY_URL = "https://dl.openfoam.com/source/v2206/ThirdParty-v2206.tgz"


def shell(cmd):
    print(
        f"{getpass.getuser()}@{socket.gethostname()}:{os.getcwd()}$ {cmd}",
    )
    cmd = f"{cmd} 2>&1"
    return subprocess.call(["/bin/bash", "-c", "set -o pipefail; " + cmd])


def git_clone_shallow(repository_url, repository_dir, branch):
    cmd = f"git clone -b {branch} --depth 1 --single-branch {repository_url} {repository_dir}"
    return 0 != shell(cmd)


def build_ls1(mamico_repo_dir, with_mpi=False, jobs=8):
    print("Building ls1-MarDyn from source...")
    ls1_dir = mamico_repo_dir / "ls1"
    if ls1_dir.exists():
        # Avoid issues with initializing submodule
        shutil.rmtree(ls1_dir)
    had_error = 0 != shell("git submodule init && git submodule update")
    if had_error != 0:
        return had_error
    had_error = False
    build_dir = ls1_dir / "build"
    build_dir.mkdir(exist_ok=True)
    cmake_args = ""
    cmake_args += f"-S{ls1_dir} -B{build_dir}"
    cmake_args += " -DENABLE_ADIOS2=OFF"
    cmake_args += f" -DENABLE_MPI={'ON' if with_mpi else 'OFF'}"
    cmake_args += " -DOPENMP=OFF"
    cmake_args += " -DENABLE_AUTOPAS=OFF"
    cmake_args += " -DENABLE_UNIT_TESTS=OFF"
    cmake_args += " -DENABLE_ALLLBL=OFF"
    cmake_args += " -DMAMICO_COUPLING=ON"
    cmake_args += f" -DMAMICO_SRC_DIR={mamico_repo_dir}"
    had_error = 0 != shell(f"cmake {cmake_args}")
    if not had_error:
        with ChangeDir(build_dir):
            had_error = 0 != shell(f"make -j {jobs}")
    return had_error


def download_open_foam():
    had_error = False
    print(f"Downloading OpenFOAM-v{OPEN_FOAM_VERSION}...")
    OPEN_FOAM_BASE_DIR.mkdir(exist_ok=True)

    def download_and_extract(url, archive_filename):
        return 0 != shell(
            f"wget -q -nc -O {archive_filename} {url} && tar -xzf {archive_filename} && rm {archive_filename}"
        )

    with ChangeDir(OPEN_FOAM_BASE_DIR):
        had_error |= download_and_extract(
            OPEN_FOAM_URL, f"OpenFOAM-v{OPEN_FOAM_VERSION}.tgz"
        )
    if not had_error:
        with ChangeDir(OPEN_FOAM_SRC_DIR):
            print(f"Downloading OpenFOAM-v{OPEN_FOAM_VERSION} third party files...")
            had_error |= download_and_extract(
                OPEN_FOAM_THIRD_PARTY_URL, f"ThirdParty-v{OPEN_FOAM_VERSION}.tgz"
            )
            if not had_error:
                shell(f"mv ThirdParty-v{OPEN_FOAM_VERSION} ThirdParty")
    return had_error


def build_open_foam(jobs=8):
    print("Building OpenFOAM from source...")
    had_error = False
    if not OPEN_FOAM_BASE_DIR.exists():
        had_error |= download_open_foam()
        if had_error:
            shutil.rmtree(OPEN_FOAM_BASE_DIR)
            return had_error
    with ChangeDir(OPEN_FOAM_SRC_DIR):
        cmd_build = ""
        cmd_build += f"source {OPEN_FOAM_SRC_DIR / 'etc' / 'bashrc'};"
        cmd_build += f" foamSystemCheck; foam; ./Allwmake -s -l -q -j {jobs}"
        had_error |= 0 != shell(
            cmd_build + " || :"
        )  # Building OpenFOAM has some irrelevant errors (non-zero exit code)
    return had_error


def build_lammps():
    print("Building LAMMPS from source...")
    had_error = False
    if not LAMMPS_REPO_DIR.exists():
        had_error = git_clone_shallow(
            LAMMPS_REPO_URL, LAMMPS_REPO_DIR, LAMMPS_REPO_BRANCH
        )
    build_dir = LAMMPS_REPO_DIR / "build"
    build_dir.mkdir(exist_ok=True)
    with ChangeDir(build_dir):
        had_error_pull = 0 != shell("git pull")  # Track latest LAMMPS release
        had_error_cmake = 0 != shell(f"cmake -DBUILD_SHARED_LIBS=ON ../cmake/")
        had_error_cmake_build = 0 != shell("cmake --build .")
        had_error_install = 0 != shell("make install")
    had_error = (
        had_error_pull | had_error_cmake | had_error_cmake_build | had_error_install
    )
    return had_error


md_solvers = dict(
    md="SIMPLE_MD",
    ls1="LS1_MARDYN",
    lammps="LAMMPS_MD",
)


def build_mamico_couette_md(
    md_solver="md", with_openfoam=False, with_mpi=False, jobs=8
):
    run_info = f"Started {time.strftime('%Y-%m-%d %H:%M:%S %Z')}"
    print(run_info)
    had_error = False
    build_dir = MAMICO_REPO_DIR / "build"
    couette_bin_path = build_dir / "couette"
    build_info = f"Building {couette_bin_path} ({md_solver})"
    if with_openfoam:
        build_info += " with OpenFOAM"
    print(f"::: {build_info} :::")
    build_dir.mkdir(exist_ok=True)
    try:
        md_solver_cmake_flag = md_solvers[md_solver]
    except:
        raise ValueError("Unknown MD solver")
    cmake_args = ""
    cmake_args += f"-S{MAMICO_REPO_DIR} -B{build_dir}"
    cmake_args += f" -DCMAKE_BUILD_TYPE={MAMICO_BUILD_TYPE}"
    cmake_args += f" -DBUILD_WITH_MPI={'ON' if with_mpi else 'OFF'}"
    cmake_args += f" -DMD_SIM={md_solver_cmake_flag}"
    environement = dict()
    environement["PKG_CONFIG_PATH"] = "$PKG_CONFIG_PATH"
    cmake_args += " -DBUILD_WITH_LAMMPS=OFF"
    pre_cmake_cmd = ""
    if md_solver == "ls1":
        had_error |= build_ls1(MAMICO_REPO_DIR, with_mpi, jobs)
        cmake_args += f" -DLS1_SRC_DIR={MAMICO_REPO_DIR / 'ls1'}"
    elif md_solver == "lammps":
        had_error |= build_lammps()
        cmake_args = cmake_args.replace("LAMMPS=OFF", "LAMMPS=ON")
        for suffix in ["", "64"]:
            pkgconfig_path = Path.home() / ".local" / f"lib{suffix}" / "pkgconfig"
            environement["PKG_CONFIG_PATH"] += f":{pkgconfig_path}"
        shell(f"ln -sf {LAMMPS_REPO_DIR}/src {MAMICO_REPO_DIR}/lammps")
    if with_openfoam:
        had_error |= build_open_foam(jobs)
    cmake_args += f" -DBUILD_WITH_OPENFOAM={'ON' if with_openfoam else 'OFF'}"
    environement_prefix = " ".join(f"{k}={v}" for k, v in environement.items())
    cmd = f"{pre_cmake_cmd}{environement_prefix} cmake {cmake_args}"
    if with_openfoam:
        cmd = f"source {OPEN_FOAM_SRC_DIR / 'etc' / 'bashrc'}; {cmd}"
    had_error |= 0 != shell(cmd)
    if not had_error:
        with ChangeDir(build_dir):
            if couette_bin_path.exists():
                # Avoid confusing old build for different one
                couette_bin_path.unlink()
            had_error = 0 != shell(f"make -j {jobs} couette")
    print(f"::: Completed (Successful: {not had_error}) :::")
    return None if had_error else couette_bin_path


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.description = (
        "Script to build the MaMiCo couette flow example using a single command."
    )
    arg_parser.add_argument(
        "-s", "--md-solver", choices=md_solvers.keys(), required=True
    )
    arg_parser.add_argument("-F", "--with-foam", action="store_true")
    arg_parser.add_argument("-M", "--with-mpi", action="store_true")
    arg_parser.add_argument("-j", "--jobs", default=8)
    args = arg_parser.parse_args()

    exec_path = build_mamico_couette_md(
        md_solver=args.md_solver,
        with_openfoam=args.with_foam,
        with_mpi=args.with_mpi,
        jobs=args.jobs,
    )
    if exec_path is not None:
        print(exec_path)
    exit(1 if exec_path is None else 0)
