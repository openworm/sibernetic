#!/usr/bin/env bash
#
# Setup script for Sibernetic on macOS or Linux.
#
# Builds the main C++ binary (OpenCL backend on Linux; on Apple Silicon the
# main binary builds without OpenCL, and the native Metal substrate is built
# separately via src/metal_diff/build.sh).
#
# Python is needed only for the embedded c302/NEURON bridge used by
# `sibernetic_c302.py`; if you don't need neural-coupled simulations you can
# ignore the Python plumbing — the binary itself still runs.
#
set -euo pipefail

SUPPORTED_PYTHONS=("python3.13" "python3.12" "python3.11" "python3.10")

find_supported_python() {
    for candidate in "${SUPPORTED_PYTHONS[@]}"; do
        if command -v "$candidate" >/dev/null 2>&1; then
            command -v "$candidate"
            return 0
        fi
    done
    return 1
}

if [[ "$(uname)" == "Darwin" ]]; then
    if ! command -v brew >/dev/null 2>&1; then
        echo "Homebrew not found. Please install it from https://brew.sh/" >&2
        exit 1
    fi
    brew update
    brew install python@3.13 glew freeglut clinfo opencl-headers pipx || true

    pipx install ruff || true
    pipx install pytest || true

    BUILD_PY="$(find_supported_python || true)"
    if [[ -z "$BUILD_PY" ]]; then
        BUILD_PY="$(brew --prefix python@3.13)/bin/python3.13"
    fi
    if [[ ! -x "$BUILD_PY" ]]; then
        echo "Error: could not locate a working Python 3.10-3.13 interpreter." >&2
        exit 1
    fi
    echo "Using $BUILD_PY ($("$BUILD_PY" --version)) for the venv and the C++ build."

    # The embedded Python (used by c302/NEURON only — not for any backend
    # solver) needs to share an ABI with the venv. Recreate the venv if the
    # interpreter version drifted.
    if [[ -d .venv ]]; then
        existing_minor="$(.venv/bin/python -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")' 2>/dev/null || echo "?")"
        required_minor="$("$BUILD_PY" -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')"
        if [[ "$existing_minor" != "$required_minor" ]]; then
            echo "Existing .venv uses Python $existing_minor but build requires $required_minor — recreating."
            rm -rf .venv
        fi
    fi
    if [[ ! -d .venv ]]; then
        "$BUILD_PY" -m venv .venv
    fi
    .venv/bin/pip install --upgrade pip

    # numpy is needed for the position-buffer renderer (render_movie.py) and
    # by the c302/NEURON bridge.
    .venv/bin/pip install numpy

    # Optional: pyneuroml only matters if you're running neural-coupled
    # simulations via sibernetic_c302.py.
    .venv/bin/pip install pyneuroml || \
        echo "Warning: pyneuroml not installed; sibernetic_c302.py will use the subprocess fallback."

    export PYTHONHEADERDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_path("include"))')"
    export PYTHONLIBDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')"
    export PYTHONFRAMEWORKDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_config_var("PYTHONFRAMEWORKPREFIX"))')"
    export NUMPYHEADERDIR="$(.venv/bin/python -c 'import numpy; print(numpy.get_include())')"

    PYTHON_CONFIG="$(dirname "$BUILD_PY")/$(basename "$BUILD_PY")-config"
    if [[ -x "$PYTHON_CONFIG" ]]; then
        export PYTHON_LIB_FLAGS="$("$PYTHON_CONFIG" --embed --ldflags)"
    fi

    make clean -f makefile.OSX
    make all -f makefile.OSX

    echo
    echo "Main Sibernetic binary built. On Apple Silicon, OpenCL is disabled"
    echo "in this binary; for GPU simulation, build the native Metal substrate:"
    echo "  cd src/metal_diff && ./build.sh"
    echo "  ./src/metal_diff/sib_metal --help"

else
    apt-get update
    apt-get install -y python3-dev python3-venv ocl-icd-opencl-dev libglu1-mesa-dev freeglut3-dev libglew-dev clinfo pocl-opencl-icd

    python3 -m pip install numpy
    python3 -m pip install ruff pytest pyneuroml || \
        echo "Warning: optional packages (ruff/pytest/pyneuroml) failed to install."
fi
