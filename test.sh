#!/bin/bash
set -ex

ruff format *.py wcon/*.py tests/*.py
ruff check *.py wcon/*.py tests/*.py

make clean
make

# Run tests using native ./Release/Sibernetic executable with OpenCL
time ./run_opencl_tests.sh $@

# Run tests using sibernetic_c302.py
time ./run_sibernetic_c302_tests.sh $@ 


