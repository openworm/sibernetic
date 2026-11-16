#!/bin/bash
set -ex

if command -v ruff >/dev/null 2>&1; then
    ruff format *.py wcon/*.py tests/*.py
    ruff check *.py wcon/*.py tests/*.py
fi

make clean
make -j4

# Run tests using native ./Release/Sibernetic executable with OpenCL
time ./run_opencl_tests.sh $@

# Run tests using sibernetic_c302.py
time ./run_sibernetic_c302_tests.sh $@ 


