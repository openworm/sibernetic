#!/usr/bin/env bash
# Regression runner for the CUDA differentiable substrate.
#
# Runs every test_*.py in this directory sequentially. Each test is
# self-contained: builds its own fixture, invokes sib_cuda (or sib_cuda.exe)
# directly, validates outputs in-process. Exit 0 iff all PASS.
#
# Mirrors the simple shell-loop style of openworm/sibernetic/run_all_tests.sh
# (top-level) so the CI conventions stay consistent.
#
# Prerequisite: sib_cuda built. Run ./build.sh first (or cmd /c build.bat on
# Windows; build.sh is unverified on Linux as of this writing).
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/../.."   # cd to repo root so relative paths in tests work

PY="${PYTHON:-python3}"
if ! command -v "$PY" >/dev/null 2>&1; then
    PY=python
fi

n_pass=0
n_fail=0
failed=()

for t in "$HERE"/test_*.py; do
    name="$(basename "$t" .py)"
    printf '=== %s ===\n' "$name"
    if "$PY" "$t"; then
        printf '[PASS] %s\n\n' "$name"
        n_pass=$((n_pass + 1))
    else
        printf '[FAIL] %s (exit $?)\n\n' "$name"
        n_fail=$((n_fail + 1))
        failed+=("$name")
    fi
done

printf '=== Summary ===\n'
printf '  passed: %d\n' "$n_pass"
printf '  failed: %d\n' "$n_fail"
if [ "$n_fail" -gt 0 ]; then
    printf '  failed tests:\n'
    for f in "${failed[@]}"; do
        printf '    - %s\n' "$f"
    done
    exit 1
fi
printf '[OVERALL PASS] all %d CUDA tests\n' "$n_pass"
exit 0
