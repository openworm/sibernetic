# Repository Guide

This repository contains the **Sibernetic** simulator along with Python
bindings and tests.  The code base mixes C++ (in `src/` and `inc/`), OpenCL
kernels and a number of helper Python scripts.

## Layout
- `src/` – main C++ sources and OpenCL kernels (`sphFluid.cl`).
- `inc/` – C++ headers with physics constants, solver classes and helpers.
- `configuration/` – example configuration files for the simulator.
- `buffers/` – output data (created at runtime).
- `tests/` – Python tests; `run_all_tests.sh` wraps `sibernetic_c302.py`
  for automated runs.
- Python utilities such as `main_sim.py`, `sibernetic_c302.py` and
  `plot_positions.py` can drive the simulator or analyse its output.

## Building and Testing
Before building you should install dependencies via `./setup.sh`.
To compile the C++ code use `make`.  A convenience script `test.sh`
runs code formatting, static checks via **ruff**, builds the simulator
and executes the test suite:

```bash
./test.sh
```

`test.sh` in turn calls `run_all_tests.sh` which performs multiple
runs of `sibernetic_c302.py` and verifies the produced output files.

For Apple Silicon GPU simulation, build the native Metal substrate
separately via `cd src/metal_diff && ./build.sh`; the resulting
`sib_metal` binary exposes the differentiable XPBD pipeline (see
README.md "Differentiable physics with native Metal").

For NVIDIA GPU simulation, build the native CUDA substrate separately
via `./src/cuda_diff/build.sh` (Linux) or `cmd /c src\cuda_diff\build.bat`
(Windows); the resulting `sib_cuda` binary exposes the same
differentiable XPBD CLI as `sib_metal` (`xpbd_full_fwd` /
`xpbd_full_bwd` / `xpbd_step`) so the shared Python harnesses
(`sgd_true.py`, `dump_*_trajectory.py`) work against either substrate
by swapping the binary path. Per-phase status, the FD-validated test
suite, and op-name maps are in `src/cuda_diff/README.md`.

## Contributing Notes
- Keep C++ headers and sources under `inc/` and `src/` respectively.
- OpenCL kernels reside in `src/*.cl`.
- When adding Python scripts ensure they pass `ruff format` and
  `ruff check`.
- Test additions should be placed under `tests/` and runnable via
  `test.sh`.

