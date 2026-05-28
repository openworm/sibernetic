#!/bin/bash
set -ex

./setup.sh

ruff format *.py wcon/*.py tests/*.py
ruff check *.py wcon/*.py tests/*.py

make clean
make

time ./run_all_tests.sh $@

