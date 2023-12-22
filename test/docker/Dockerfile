# Used to generate ghcr.io/hsu-hpc/testmamico
FROM ubuntu:jammy

RUN apt -y update && apt -y upgrade
RUN apt install -y build-essential cmake git libcppunit-dev openmpi-bin openmpi-common libopenmpi-dev libeigen3-dev clang gcc ccache
RUN apt install -y pybind11-dev python3-pybind11 python-is-python3 python3 python3-dev python-dev-is-python3
ENV PATH "/usr/lib/ccache:$PATH"
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1
RUN apt install -y --no-install-recommends doxygen
RUN apt install -y graphviz lcov