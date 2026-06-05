#!/bin/bash
set -ex

#./setup.sh  # Run this outside test.sh

ruff format *.py wcon/*.py tests/*.py
ruff check *.py wcon/*.py tests/*.py

make clean
make

time ./run_opencl_tests.sh $@
time ./run_all_tests.sh $@

