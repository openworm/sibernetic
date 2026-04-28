#!/usr/bin/env bash
#
# Setup script for Sibernetic on macOS or linux.
#
set -euo pipefail

# Install system dependencies for building Sibernetic and the PyTorch solver
if [[ "$(uname)" == "Darwin" ]]; then
    if ! command -v brew >/dev/null 2>&1; then
        echo "Homebrew not found. Please install it from https://brew.sh/" >&2
        exit 1
    fi
    brew update
    brew install python glew freeglut clinfo opencl-headers pipx || true

    # CLI tools (run as commands, not imported) — pipx-isolated venvs are fine
    pipx install ruff || true
    pipx install pytest || true

    # Runtime libraries embedded by the C++ binary must be importable from the
    # same Python that Sibernetic links against. pipx isolates each package in
    # its own venv, which makes them invisible to the embedded interpreter, so
    # use a project-local venv on top of Homebrew's python3 instead.
    if [[ ! -d .venv ]]; then
        python3 -m venv .venv
    fi
    .venv/bin/pip install --upgrade pip
    .venv/bin/pip install torch taichi numpy pyneuroml || \
        echo "Warning: failed to install one or more runtime Python packages"

    export PYTHONHEADERDIR="$(python3 -c 'import sysconfig; print(sysconfig.get_path("include"))')"
    export PYTHONLIBDIR=$(python3 -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')
    export PYTHONFRAMEWORKDIR=$(python3 -c 'import sysconfig; print(sysconfig.get_config_var("PYTHONFRAMEWORKPREFIX"))')

    make clean -f makefile.OSX
    make all -f makefile.OSX

    echo
    echo "Build complete. Run from the repo root, e.g.:"
    echo "  ./Release/Sibernetic -f worm backend=taichi"

else
    apt-get update
    apt-get install -y python3-dev ocl-icd-opencl-dev libglu1-mesa-dev freeglut3-dev libglew-dev clinfo pocl-opencl-icd

    # Install Python packages
    python3 -m pip install torch taichi ruff pytest pyneuroml || echo "Warning: failed to install pyneuroml"

fi

