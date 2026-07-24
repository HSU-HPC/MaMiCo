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
MAMICO_BUILD_TYPE_DEFAULT = "DebugOptimized"
MAMICO_BUILD_TYPE_DEBUG = "DebugMD"

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
        flush=True
    )
    cmd = f"{cmd} 2>&1"
    return subprocess.call(["/bin/bash", "-c", "set -o pipefail; " + cmd])


def git_clone_shallow(repository_url, repository_dir, branch):
    cmd = f"git clone -b {branch} --depth 1 --single-branch {repository_url} {repository_dir}"
    return 0 != shell(cmd)


def build_ls1(mamico_repo_dir, with_mpi=False, jobs=8, force_gcc=False):
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
    if force_gcc:
        cmake_args += f" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
    cmake_args += f" -S{ls1_dir} -B{build_dir}"
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
                had_error = 0 != shell(f"mv ThirdParty-v{OPEN_FOAM_VERSION} ThirdParty")
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
    if not had_error:
        # Rename unused logging macro to avoid issue when compiling MaMiCo with both OpenFOAM and ls1
        had_error = 0 != shell(f"sed -i 's/^#define Log/#define FoamLog/g' {OPEN_FOAM_SRC_DIR}/src/OpenFOAM/lnInclude/messageStream.H")
    return had_error


def build_lammps(with_mpi, jobs=8):
    print("Building LAMMPS from source...")
    had_error = False
    if not LAMMPS_REPO_DIR.exists():
        had_error = git_clone_shallow(
            LAMMPS_REPO_URL, LAMMPS_REPO_DIR, LAMMPS_REPO_BRANCH
        )
    build_dir = LAMMPS_REPO_DIR / "build"
    install_dir = LAMMPS_REPO_DIR / "install"
    build_dir.mkdir(exist_ok=True)
    mamico_fix_dir = MAMICO_REPO_DIR / "coupling" / "interface" / "impl" / "LAMMPS" / "USER-MAMICO"
    with ChangeDir(build_dir):
        had_error_pull = 0 != shell("git pull")  # Track latest LAMMPS release
        shutil.copytree(mamico_fix_dir, LAMMPS_REPO_DIR / "src" / "USER-MAMICO", dirs_exist_ok=True)
        with open(LAMMPS_REPO_DIR / "cmake" / "CMakeLists.txt", 'r+') as file:
            cmakefile = file.readlines()
            inject = 0
            for line in cmakefile:
                if line.strip() == 'set(STANDARD_PACKAGES':
                    break
                inject += 1
            cmakefile.insert(inject + 1, 'USER-MAMICO' + os.linesep)
            file.seek(0)
            file.writelines(cmakefile)
        had_error_cmake = 0 != shell("cmake -DPKG_KOKKOS=ON -DPKG_USER-MAMICO=ON -DBUILD_SHARED_LIBS=ON"
                                     f" -DBUILD_MPI={'ON' if with_mpi else 'OFF'}"
                                     f" -DCMAKE_CXX_FLAGS=-I\\ {MAMICO_REPO_DIR} ../cmake/")
        had_error_cmake_build = 0 != shell(f"make -j {jobs}")
    had_error = (
        had_error_pull | had_error_cmake | had_error_cmake_build
    )
    return had_error


md_solvers = dict(
    md="SIMPLE_MD",
    ls1="LS1_MARDYN",
    lammps="LAMMPS_MD",
)


# FIXME: LAMMPS is not found during linking
def build_mamico_couette(
    md_solver="md", with_openfoam=False, with_mpi=False, with_tests=False, is_debug=False, jobs=8, clean=False, force_gcc=False
):
    run_info = f"Started {time.strftime('%Y-%m-%d %H:%M:%S %Z')}"
    print(run_info)
    had_error = False
    build_dir = MAMICO_REPO_DIR / "build"
    if clean:
        shutil.rmtree(build_dir, ignore_errors=True)
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
    if with_openfoam:
        had_error |= build_open_foam(jobs)
        if had_error:
            return None
    cmake_args = ""
    if force_gcc:
        cmake_args += f" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
    #cmake_args += f"-S{MAMICO_REPO_DIR} -B{build_dir}"
    cmake_args += f" -DCMAKE_BUILD_TYPE={MAMICO_BUILD_TYPE_DEBUG if is_debug else MAMICO_BUILD_TYPE_DEFAULT}"
    cmake_args += f" -DBUILD_WITH_MPI={'ON' if with_mpi else 'OFF'}"
    cmake_args += f" -DBUILD_TESTING={'ON' if with_tests else 'OFF'}"
    cmake_args += f" -DMD_SIM={md_solver_cmake_flag}"
    environement = dict(
        PKG_CONFIG_PATH = os.getenv("PKG_CONFIG_PATH", "")
    )
    cmake_args += " -DBUILD_WITH_LAMMPS=OFF"
    if md_solver == "ls1":
        had_error |= build_ls1(MAMICO_REPO_DIR, with_mpi, jobs, force_gcc)
        cmake_args += f" -DLS1_SRC_DIR={MAMICO_REPO_DIR / 'ls1'}"
    elif md_solver == "lammps":
        had_error |= build_lammps(with_mpi, jobs)
        cmake_args = cmake_args.replace("LAMMPS=OFF", "LAMMPS=ON")
        cmake_args += f" -DLAMMPS_DIR={LAMMPS_REPO_DIR}"
    cmake_args += f" -DBUILD_WITH_OPENFOAM={'ON' if with_openfoam else 'OFF'}"
    environement_prefix = " ".join(f"{k}={v}" for k, v in environement.items())
    cmd = f"{environement_prefix} cmake {cmake_args} .."
    if with_openfoam:
        cmd = f"source {OPEN_FOAM_SRC_DIR / 'etc' / 'bashrc'}; {cmd}"
    with ChangeDir(build_dir):
        had_error |= 0 != shell(cmd)
        if not had_error:
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
    arg_parser.add_argument("-T", "--with-tests", action="store_true")
    arg_parser.add_argument("-d", "--debug", action="store_true")
    arg_parser.add_argument("-j", "--jobs", default=8)
    arg_parser.add_argument("-g", "--force-gcc", action="store_true")
    arg_parser.add_argument("-c", "--clean", action="store_true")
    args = arg_parser.parse_args()

    exec_path = build_mamico_couette(
        md_solver=args.md_solver,
        with_openfoam=args.with_foam,
        with_mpi=args.with_mpi,
        with_tests=args.with_tests,
        is_debug=args.debug,
        jobs=args.jobs,
        clean=args.clean,
        force_gcc=args.force_gcc
    )
    if exec_path is not None:
        print(exec_path)
    exit(1 if exec_path is None else 0)
