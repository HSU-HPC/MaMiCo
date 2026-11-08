# Used to generate ghcr.io/hsu-hpc/testmamico
FROM ubuntu:26.04

# Install all dependencies in a single layer to keep image small
RUN apt -y update && \
    apt -y upgrade && \
    apt install --no-install-recommends \
        build-essential \
        cmake \
        gcc \
        git \
        libcppunit-dev \
        ccache \
        libeigen3-dev \
        clang-format-18 \
        `# MPI` \
        openmpi-bin \
        openmpi-common \
        libopenmpi-dev \
        `# Clang` \
        clang \
        g++-16 \
        libomp-dev \
        `# Python` \
        pybind11-dev \
        python3-pybind11 \
        python-is-python3 \
        python3 \
        python3-dev \
        python-dev-is-python3 \
        python3-pip \
        `# Coverage/documentation` \
        doxygen \
        graphviz \
        lcov \
        -y && \
    rm -rf /var/lib/apt/lists/* `# Remove package lists`

RUN ln -s /usr/bin/clang-format-18 /usr/local/bin/clang-format
ENV PATH="/usr/lib/ccache:$PATH"
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Python dependencies for reporting coverage and formatting
RUN pip install \
        pandas==3.0.5 \
        beautifulsoup4==4.15.0 \
        black==26.5.1 \
        cmakelang==0.6.13 \
        --break-system-packages `# Install at system-level`
