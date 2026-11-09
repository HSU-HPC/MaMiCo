# Two-stage couette build for network-isolated compute nodes

An alternative to `MaMiCo/tools/build-couette.py`, which is left untouched.
Split so that everything needing the internet happens on a login node and the
real build happens on a compute node with no outbound network.

| | script | needs network | produces |
|---|---|---|---|
| Stage 1 | `build-couette-login.py` | yes | `build/LAMMPS` (git clone), `build/_deps/kokkos-src` (FetchContent tarball), `build/.offline-prefetch.json` |
| Stage 2 | `build-couette-compute.py` | no | `build/couette` |

Both locate the MaMiCo tree via `-S/--source`, else `$MAMICO_DIR`, else the git
repo enclosing the script, else the one enclosing the cwd. Paths containing
spaces are supported.

## Why two stages

The only network dependencies in a LAMMPS+OpenFOAM couette build are:

- the LAMMPS `git clone` / `git pull` in `build_lammps()`
- MaMiCo's `FetchContent_Declare(Kokkos URL .../4.6.02.tar.gz)` at
  `CMakeLists.txt:156`, populated at *configure* time into `build/_deps/kokkos-src`
- optionally the OpenFOAM v2206 tarballs, which you do not need when OpenFOAM
  comes from spack or a module

LAMMPS' own Kokkos is the bundled `lib/kokkos` (`DOWNLOAD_KOKKOS:BOOL=OFF`), so
once the clone exists LAMMPS builds fully offline.

## What stage 1 writes

Stage 1 only writes network-derived sources and the manifest into `build/`:
`LAMMPS/`, `_deps/` and `openfoam/` (with `--download-foam`), plus
`.offline-prefetch.json`. The Kokkos-populating CMake configure runs in a
temporary directory; only the fetched `_deps/*-src` trees are moved into
`build/_deps` before it is removed. No `CMakeCache.txt` or object code is
created in `build/`, so an existing build there is not reconfigured or
overwritten.

## What stage 2 does and does not keep

Stage 2 configures with `FETCHCONTENT_SOURCE_DIR_KOKKOS` pointing at the
prefetched source plus `FETCHCONTENT_FULLY_DISCONNECTED=ON`, so CMake uses the
tarball in place and never tries to download. It never runs `git pull`.

Kept: `LAMMPS/`, `_deps/*-src`, `_deps/*-subbuild`, `openfoam/`.

Deleted: `CMakeCache.txt`, `CMakeFiles/`, `_deps/*-build`, `LAMMPS/build`,
`LAMMPS/install`, and the other regenerable CMake artefacts.

`_deps/kokkos-build` and `LAMMPS/build` go by default because Kokkos is
configured with `Kokkos_ARCH_NATIVE ON` — object code compiled on a login node
can target the wrong microarchitecture. `--keep-kokkos-build` and
`--keep-lammps-build` override that.

Anything else in `build/` (CSVs, `log.lammps*`, `couette.xml`, `FoamSetup/`) is
**moved** into `build/_pruned_<timestamp>/`, never deleted.

Dropping `CMakeCache.txt` is also the fix for a cached `MPI_CXX_WORKS` failure:
`find_package(MPI)` caches its verdict, so a failed probe survives every
environment fix until the cache goes.

## Safety

Stage 2 refuses to prune when a process has the build directory as its cwd or
exe, or when anything in it changed within `--idle-minutes` (default 10). It
also asks for confirmation on a TTY. Use `--dry-run` first, `-y` in a batch
script, `--force` only deliberately.

## Usage

Login node:

```bash
module load stack/26.0.0 user_spack/26.0.0 gcc/14.3.0 cmake/3.31.11   # example module names
python3 build-couette-login.py -S ~/MaMiCo
```

Add `--download-foam` only if OpenFOAM is not provided by spack or a module,
and `--with-ls1` for the ls1-MarDyn submodule.

Compute node:

```bash
spack load openfoam        # or source its etc/bashrc
export I_MPI_CXX=g++ I_MPI_CC=gcc
python3 build-couette-compute.py -S ~/MaMiCo -s <md|ls1|lammps> -M -F -j 16 -y
```

`--mpi-cxx` / `--mpi-cc` default to `$I_MPI_ROOT/bin/{mpicxx,mpicc}`. Set
`--cxx`/`--cc` to the same compiler your MPI and OpenFOAM were built with.

`-F` requires OpenFOAM in the environment (`WM_PROJECT_DIR`), or a tree under
`build/openfoam/` whose `etc/bashrc` the script will source. It never downloads.
