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
set -euox pipefail

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

    ###########################################################################################################
    # NOTE: This is not yet tested on most Sibernetic branches, Work in progress, targeted at new backends! 
    ###########################################################################################################

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
    echo "Installing Sibernetic on Linux."

    # Use sudo for privileged commands if it's available and we're not already
    # root; otherwise run them directly as the current user (e.g. in a container).
    RUN_SUDO=""
    if [[ "$(id -u)" -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
        RUN_SUDO="sudo"
    fi

    $RUN_SUDO apt update
    $RUN_SUDO apt install -y --only-upgrade python3-pip

    lscpu

    echo "Installing OpenCL Drivers and Java (for jNeuroML)"
        
    $RUN_SUDO apt-get update --fix-missing
    $RUN_SUDO apt-get install -y \
          build-essential \
          git \
          python3-dev \
          default-jre-headless \
          libopengl0 \
          freeglut3-dev \
          libglu1-mesa-dev \
          clinfo \
          ocl-icd-libopencl1 \
          ocl-icd-opencl-dev \
          opencl-clhpp-headers \
          pocl-opencl-icd \
          mesa-opencl-icd

    echo "OpenCL Driver Installation Complete"

    echo "CLINFO:"
    clinfo

    $RUN_SUDO pip install --ignore-installed -r requirements.txt

    export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}."

    make clean
    make -j4


    echo "Sibernetic installation complete. You can now run the main binary:"
    echo "  ./Release/Sibernetic --help"
fi
