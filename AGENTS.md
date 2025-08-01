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
To compile the C++ code use `make`.  A convenience script `test.sh`
runs code formatting, static checks via **ruff**, builds the simulator
and executes the test suite:

```bash
./test.sh
```

`test.sh` in turn calls `run_all_tests.sh` which performs multiple
runs of `sibernetic_c302.py` and verifies the produced output files.

## Contributing Notes
- Keep C++ headers and sources under `inc/` and `src/` respectively.
- OpenCL kernels reside in `src/*.cl`.
- When adding Python scripts ensure they pass `ruff format` and
  `ruff check`.
- Test additions should be placed under `tests/` and runnable via
  `test.sh`.

