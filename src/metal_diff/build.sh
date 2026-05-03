#!/bin/bash
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"
clang++ -std=c++17 -O2 \
    -framework Metal -framework Foundation \
    sib_metal.mm -o sib_metal
echo "Built $DIR/sib_metal"
