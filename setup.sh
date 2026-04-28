#!/usr/bin/env bash
#
# Setup script for Sibernetic on macOS or linux.
#
# Taichi 1.7.4 ships pre-built wheels only for CPython 3.7-3.13. Python 3.14+
# has no wheel and a from-source build needs LLVM 15 + matching toolchain, so
# this script pins Sibernetic to a Python in the 3.10-3.13 range. If you need
# a newer Python, see the "Python 3.14 and beyond" section in the README.
#
set -euo pipefail

# Taichi 1.7.4 wheel coverage — keep newest first so the script prefers
# the most recent supported interpreter.
SUPPORTED_PYTHONS=("python3.13" "python3.12" "python3.11" "python3.10")

# Find a usable interpreter on PATH; returns empty string if none match.
find_supported_python() {
    for candidate in "${SUPPORTED_PYTHONS[@]}"; do
        if command -v "$candidate" >/dev/null 2>&1; then
            command -v "$candidate"
            return 0
        fi
    done
    return 1
}

verify_runtime_imports() {
    local py="$1"
    if ! "$py" -c "import taichi, torch, numpy" 2>/dev/null; then
        echo "Error: taichi/torch/numpy are installed but not importable from $py" >&2
        echo "Diagnostics:" >&2
        "$py" -c "import sys; print('  python:', sys.executable, sys.version)" >&2
        "$py" -m pip show taichi torch numpy 2>&1 | sed 's/^/  /' >&2
        return 1
    fi
}

if [[ "$(uname)" == "Darwin" ]]; then
    if ! command -v brew >/dev/null 2>&1; then
        echo "Homebrew not found. Please install it from https://brew.sh/" >&2
        exit 1
    fi
    brew update
    # Pin to python@3.13 explicitly so `brew install python` flipping to 3.14+
    # in a future Homebrew release doesn't silently break Taichi.
    brew install python@3.13 glew freeglut clinfo opencl-headers pipx || true

    # CLI tools (run as commands, not imported) — pipx-isolated venvs are fine
    pipx install ruff || true
    pipx install pytest || true

    # Pick a supported Python; install python@3.13 if none is on PATH.
    BUILD_PY="$(find_supported_python || true)"
    if [[ -z "$BUILD_PY" ]]; then
        echo "No python 3.10-3.13 on PATH after brew install python@3.13."
        echo "Falling back to the Homebrew-managed python@3.13."
        BUILD_PY="$(brew --prefix python@3.13)/bin/python3.13"
    fi
    if [[ ! -x "$BUILD_PY" ]]; then
        echo "Error: could not locate a working Python 3.10-3.13 interpreter." >&2
        echo "Install one (brew install python@3.13) and re-run setup.sh." >&2
        exit 1
    fi
    echo "Using $BUILD_PY ($("$BUILD_PY" --version)) for the venv and the C++ build."

    # Runtime libraries are embedded by the C++ binary, so they must live in a
    # venv tied to the SAME Python the binary links against. Recreate the venv
    # if it was created with a different Python version (binary extensions
    # like numpy/torch/taichi are ABI-incompatible across minor versions).
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

    # Hard requirements — fail loudly if any one of these can't be installed.
    .venv/bin/pip install torch taichi numpy

    # Optional: pyneuroml is only needed for neural-coupled simulations and
    # the test suite has a subprocess fallback.
    .venv/bin/pip install pyneuroml || \
        echo "Warning: pyneuroml not installed; sibernetic_c302.py will use the subprocess fallback."

    verify_runtime_imports .venv/bin/python

    # Use the same Python for the C++ build so the embedded interpreter and
    # the venv's binary extensions agree on Python version + ABI.
    export PYTHONHEADERDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_path("include"))')"
    export PYTHONLIBDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')"
    export PYTHONFRAMEWORKDIR="$("$BUILD_PY" -c 'import sysconfig; print(sysconfig.get_config_var("PYTHONFRAMEWORKPREFIX"))')"

    make clean -f makefile.OSX
    make all -f makefile.OSX

    echo
    echo "Build complete. Run from the repo root, e.g.:"
    echo "  ./Release/Sibernetic -f worm backend=taichi"

else
    apt-get update
    apt-get install -y python3-dev python3-venv ocl-icd-opencl-dev libglu1-mesa-dev freeglut3-dev libglew-dev clinfo pocl-opencl-icd

    # Warn (don't abort) if the system Python is outside Taichi's wheel range.
    py_minor="$(python3 -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')"
    case "$py_minor" in
        3.10|3.11|3.12|3.13) ;;
        *) echo "Warning: Python $py_minor is outside Taichi 1.7.4's supported range (3.10-3.13). Install one of those versions if taichi fails to install." >&2 ;;
    esac

    # Hard requirements
    python3 -m pip install torch taichi numpy

    # Optional + dev tools
    python3 -m pip install ruff pytest pyneuroml || \
        echo "Warning: optional packages (ruff/pytest/pyneuroml) failed to install."

    verify_runtime_imports python3
fi
