#! /usr/bin/env python3

"""
Login-node half of the two-stage MaMiCo couette build.

Does everything that needs the internet: clones LAMMPS and runs a minimal
configure in a scratch directory so that the CMake FetchContent population of
Kokkos happens and lands in build/_deps/kokkos-src. Only network-derived
sources (LAMMPS/, _deps/, openfoam/) and the manifest are written to build/;
no CMake cache or object code is, so an existing build there is left alone.
build-couette-compute.py does the real build on the compute node, where
Kokkos_ARCH_NATIVE picks up the right microarchitecture.

An alternative to build-couette.py for clusters whose compute nodes have no
outbound network; build-couette.py itself is unchanged.
"""

import argparse
import getpass
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

LAMMPS_REPO_URL = "https://github.com/lammps/lammps.git"
LAMMPS_REPO_BRANCH = "release"

OPEN_FOAM_VERSION = "2206"
OPEN_FOAM_URL = f"https://dl.openfoam.com/source/v{OPEN_FOAM_VERSION}/OpenFOAM-v{OPEN_FOAM_VERSION}.tgz"
OPEN_FOAM_THIRD_PARTY_URL = f"https://dl.openfoam.com/source/v{OPEN_FOAM_VERSION}/ThirdParty-v{OPEN_FOAM_VERSION}.tgz"

MANIFEST_NAME = ".offline-prefetch.json"


class ChangeDir:
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir

    def __enter__(self):
        self.last_dir = os.getcwd()
        os.chdir(self.tmp_dir)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.last_dir)


def shell(cmd):
    print(f"{getpass.getuser()}@{socket.gethostname()}:{os.getcwd()}$ {cmd}", flush=True)
    return subprocess.call(["/bin/bash", "-c", "set -o pipefail; " + cmd + " 2>&1"])


def git_repo_root(start):
    start = Path(start)
    if not start.is_dir():
        start = start.parent
    if not start.is_dir():
        return None
    with ChangeDir(start):
        root = subprocess.getoutput("git rev-parse --show-toplevel")
    return None if root.startswith("fatal:") else Path(root)


def resolve_source(explicit):
    """-S > $MAMICO_DIR > repo holding this script > repo holding cwd."""
    candidates = [
        Path(explicit).resolve() if explicit else None,
        Path(os.environ["MAMICO_DIR"]).resolve() if os.environ.get("MAMICO_DIR") else None,
        git_repo_root(Path(__file__).resolve()),
        git_repo_root(Path.cwd()),
    ]
    for candidate in candidates:
        if candidate is None:
            continue
        if (candidate / "CMakeLists.txt").is_file() and (candidate / "coupling").is_dir():
            return candidate
    sys.exit(
        "Could not locate the MaMiCo source tree.\n"
        "Pass -S/--source <path-to-MaMiCo> or export MAMICO_DIR."
    )


def fetch_lammps(lammps_dir, pull):
    if not lammps_dir.exists():
        print(f"::: Cloning LAMMPS ({LAMMPS_REPO_BRANCH}) :::")
        cmd = shlex.join([
            "git", "clone", "-b", LAMMPS_REPO_BRANCH, "--depth", "1", "--single-branch",
            LAMMPS_REPO_URL, str(lammps_dir),
        ])
        if shell(cmd) != 0:
            return True
    elif pull:
        print("::: Updating existing LAMMPS checkout :::")
        with ChangeDir(lammps_dir):
            if shell("git pull") != 0:
                return True
    else:
        print("::: LAMMPS checkout present, skipping pull (--no-pull) :::")
    return not (lammps_dir / "src").is_dir()


def fetch_ls1(source_dir):
    print("::: Fetching ls1-MarDyn submodule :::")
    with ChangeDir(source_dir):
        return shell("git submodule init && git submodule update") != 0


def fetch_open_foam(build_dir):
    base = build_dir / "openfoam"
    src = base / f"OpenFOAM-v{OPEN_FOAM_VERSION}"
    if src.is_dir():
        print(f"::: OpenFOAM source already present at {src} :::")
        return False
    print(f"::: Downloading OpenFOAM-v{OPEN_FOAM_VERSION} :::")
    base.mkdir(parents=True, exist_ok=True)

    def grab(url, archive):
        return shell(f"wget -q -nc -O {archive} {url} && tar -xzf {archive} && rm {archive}") != 0
    with ChangeDir(base):
        if grab(OPEN_FOAM_URL, f"OpenFOAM-v{OPEN_FOAM_VERSION}.tgz"):
            return True
    with ChangeDir(src):
        if grab(OPEN_FOAM_THIRD_PARTY_URL, f"ThirdParty-v{OPEN_FOAM_VERSION}.tgz"):
            return True
        if shell(f"mv ThirdParty-v{OPEN_FOAM_VERSION} ThirdParty") != 0:
            return True
    return False


def prime_kokkos(source_dir, build_dir, force_gcc, kokkos_target):
    """Configure the simplest possible variant in a scratch directory.

    The point is FetchContent_MakeAvailable(Kokkos), which populates
    _deps/kokkos-src at configure time. Everything is populated under the
    scratch directory; only the fetched *-src trees are moved into build/_deps
    before it is removed. MPI, LAMMPS and OpenFOAM are all off so that this
    step cannot fail for reasons unrelated to fetching.
    """
    print("::: Priming Kokkos via a minimal SimpleMD configure :::")
    args = [
        "-DCMAKE_BUILD_TYPE=DebugOptimized",
        f"-DKOKKOS_TARGET={kokkos_target}",
        "-DBUILD_WITH_MPI=OFF",
        "-DBUILD_WITH_LAMMPS=OFF",
        "-DBUILD_WITH_OPENFOAM=OFF",
        "-DBUILD_TESTING=OFF",
        "-DMD_SIM=SIMPLE_MD",
    ]
    if force_gcc:
        args += ["-DCMAKE_CXX_COMPILER=g++", "-DCMAKE_C_COMPILER=gcc"]
    scratch = Path(tempfile.mkdtemp(prefix="mamico-prime-"))
    try:
        with ChangeDir(scratch):
            failed = shell(shlex.join(["cmake", *args, str(source_dir)])) != 0
        deps = build_dir / "_deps"
        for fetched in sorted((scratch / "_deps").glob("*-src")):
            if not (fetched / "CMakeLists.txt").is_file():
                continue
            deps.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(deps / fetched.name, ignore_errors=True)
            shutil.move(str(fetched), str(deps / fetched.name))
        return failed
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def write_manifest(build_dir, source_dir, lammps_dir):
    commit = "unknown"
    if (lammps_dir / ".git").exists():
        with ChangeDir(lammps_dir):
            commit = subprocess.getoutput("git rev-parse HEAD").strip()
    manifest = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "host": socket.gethostname(),
        "user": getpass.getuser(),
        "source_dir": str(source_dir),
        "lammps_commit": commit,
        "kokkos_src": str(build_dir / "_deps" / "kokkos-src"),
    }
    (build_dir / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def verify(build_dir, lammps_dir, want_foam_src):
    problems = []
    if not (lammps_dir / "src" / "lammps.cpp").is_file():
        problems.append(f"LAMMPS source missing under {lammps_dir}")
    if not (lammps_dir / "lib" / "kokkos" / "CMakeLists.txt").is_file():
        problems.append("LAMMPS bundled lib/kokkos missing, so PKG_KOKKOS would need the network")
    kokkos_src = build_dir / "_deps" / "kokkos-src"
    if not (kokkos_src / "CMakeLists.txt").is_file():
        problems.append(f"Kokkos was not populated at {kokkos_src}")
    if want_foam_src and not (build_dir / "openfoam" / f"OpenFOAM-v{OPEN_FOAM_VERSION}").is_dir():
        problems.append("OpenFOAM source tree missing")
    return problems


def main():
    parser = argparse.ArgumentParser(
        description="Login-node prefetch and priming build for the MaMiCo couette example (stage 1 of 2)."
    )
    parser.add_argument("-S", "--source", help="MaMiCo source dir (default: $MAMICO_DIR or the enclosing git repo)")
    parser.add_argument("-g", "--force-gcc", action="store_true")
    parser.add_argument("-k", "--kokkos-target", default="OPENMP", choices=["CUDA", "OPENMP", "SERIAL"])
    parser.add_argument("--no-pull", action="store_true", help="Do not 'git pull' an existing LAMMPS checkout")
    parser.add_argument("--with-ls1", action="store_true", help="Also init/update the ls1-MarDyn submodule")
    parser.add_argument(
        "--download-foam",
        action="store_true",
        help="Download OpenFOAM sources too. Not needed when OpenFOAM comes from spack or a module.",
    )
    args = parser.parse_args()

    source_dir = resolve_source(args.source)
    build_dir = source_dir / "build"
    lammps_dir = build_dir / "LAMMPS"
    print(f"::: MaMiCo source : {source_dir}")
    print(f"::: Build dir     : {build_dir}")
    build_dir.mkdir(parents=True, exist_ok=True)

    had_error = fetch_lammps(lammps_dir, pull=not args.no_pull)
    if not had_error and args.with_ls1:
        had_error = fetch_ls1(source_dir)
    if not had_error and args.download_foam:
        had_error = fetch_open_foam(build_dir)
    if not had_error:
        had_error = prime_kokkos(source_dir, build_dir, args.force_gcc, args.kokkos_target)

    problems = verify(build_dir, lammps_dir, args.download_foam)
    if problems:
        print("\n::: PREFETCH INCOMPLETE :::")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    if had_error:
        print("\n::: The priming configure failed, but every network artefact is present. :::")
        print("::: The compute-node stage can still proceed. :::")

    manifest = write_manifest(build_dir, source_dir, lammps_dir)
    compute = Path(__file__).resolve().parent / "build-couette-compute.py"
    print("\n::: Prefetch complete :::")
    print(f"  LAMMPS commit : {manifest['lammps_commit']}")
    print(f"  Kokkos source : {manifest['kokkos_src']}")
    print(f"  Manifest      : {build_dir / MANIFEST_NAME}")
    print("\nNext, on a compute node:")
    print(f"  python3 {shlex.quote(str(compute))} -S {shlex.quote(str(source_dir))} -s <md|ls1|lammps> -M -F")
    return 0


if __name__ == "__main__":
    sys.exit(main())
