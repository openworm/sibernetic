#!/bin/bash
set -ex

ruff format *.py
ruff check *.py

mv *dat /tmp || true

omv test -V .test.cube.omt
omv test -V .test.demo1.omt

if [[ ($# -eq 1) && ($1 == '-all') ]]; then

echo "Running longer tests"

fi


