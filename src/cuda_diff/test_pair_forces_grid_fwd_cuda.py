"""Smoke test for the CUDA pair_forces_grid_fwd kernel.

Mirrors the FORWARD portion of src/metal_diff/test_pair_forces_bwd.py
(setup: 2 active + 4 boundary particles, h=3.34, sim_scale=7.4e-6).
Windows-portable: uses tempfile.gettempdir() instead of /tmp, and the
grid_origin path is passed explicitly on the CLI (not hardcoded).

Pass criterion: the kernel runs end-to-end without crashing, the
returned ext_accel is non-NaN, and its magnitude lands in a sensible
range (between 1e-8 and 1e8). FD validation against the backward
kernel is sub-task C, not this test. Local-only.
"""
import os
import platform
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()


def write_bin(path, arr, dtype=np.float32):
    arr.astype(dtype).tofile(path)


def make_grid(pos_static, h):
    """Mini grid build for a few static particles (mirrors
    metal_diff/test_pair_forces_bwd.py::make_grid)."""
    box_min = pos_static.min(axis=0) - 0.1
    box_max = pos_static.max(axis=0) + 0.1
    grid_dim = np.ceil((box_max - box_min) / h).astype(np.int32)
    n_cells = int(grid_dim.prod())
    cells = np.floor((pos_static - box_min) / h).astype(np.int32)
    cell_ids = (cells[:, 0]
                + cells[:, 1] * grid_dim[0]
                + cells[:, 2] * grid_dim[0] * grid_dim[1])
    perm = np.argsort(cell_ids, kind='stable')
    sorted_pos = pos_static[perm].astype(np.float32)
    sorted_ids = cell_ids[perm]
    counts = np.bincount(sorted_ids, minlength=n_cells)
    cell_start = np.zeros(n_cells + 1, dtype=np.int32)
    cell_start[1:] = np.cumsum(counts)
    return sorted_pos, cell_start, grid_dim, box_min.astype(np.float32)


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    print("Test: pair_forces_grid_fwd smoke test\n")

    # Scenario mirrors metal_diff/test_pair_forces_bwd.py lines 91-146.
    sim_scale = 7.4e-6
    h = 3.34
    mass = 2.0e-12
    visc_pair_coef = 1e-4

    n_active = 2
    pos_a = np.array([[1.0, 1.0, 1.0],
                      [1.0, 1.0, 2.5]], dtype=np.float32)
    vel_a = np.array([[0.0, -0.5, 0.0],
                      [0.0, +0.3, 0.1]], dtype=np.float32)

    pos_static = np.array([[0.0, 0.0, 0.0],
                           [2.5, 0.0, 0.0],
                           [0.0, 0.0, 3.0],
                           [2.5, 0.0, 3.0]], dtype=np.float32)
    n_static = len(pos_static)

    sorted_static, cell_start, grid_dim, grid_origin = make_grid(
        pos_static, h)

    # Density: hand-picked reasonable values (same as Metal test).
    density = np.array([4.05e-13, 4.05e-13], dtype=np.float32)

    paths = {
        'pos_a':         os.path.join(TMP, 'cuda_pair_pos_a.bin'),
        'vel_a':         os.path.join(TMP, 'cuda_pair_vel_a.bin'),
        'sorted_static': os.path.join(TMP, 'cuda_pair_sorted_static.bin'),
        'cell_start':    os.path.join(TMP, 'cuda_pair_cell_start.bin'),
        'density':       os.path.join(TMP, 'cuda_pair_density.bin'),
        'grid_origin':   os.path.join(TMP, 'cuda_pair_grid_origin.bin'),
        'ext_accel':     os.path.join(TMP, 'cuda_pair_ext_accel.bin'),
    }
    write_bin(paths['pos_a'], pos_a)
    write_bin(paths['vel_a'], vel_a)
    write_bin(paths['sorted_static'], sorted_static)
    cell_start.astype(np.int32).tofile(paths['cell_start'])
    write_bin(paths['density'], density)
    write_bin(paths['grid_origin'], grid_origin)

    cmd = [BINARY, "pair_forces_grid_fwd",
           str(n_active), str(n_static),
           str(h), str(mass), str(sim_scale), str(visc_pair_coef),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(int(grid_dim[0])), str(int(grid_dim[1])), str(int(grid_dim[2])),
           paths['grid_origin'], paths['ext_accel']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr, file=sys.stderr)
        print("STDOUT:", r.stdout, file=sys.stderr)
        return 1

    ext = np.fromfile(paths['ext_accel'], dtype=np.float32).reshape(-1, 3)
    print(f"ext_accel (shape={ext.shape}):")
    print(ext)
    mag = float(np.linalg.norm(ext))
    print(f"  ||ext_accel|| = {mag:.3e}")

    if np.isnan(ext).any():
        print("[FAIL] NaN in ext_accel")
        return 1
    if np.isinf(ext).any():
        print("[FAIL] Inf in ext_accel")
        return 1
    # Sanity range: a real number, neither vanishing nor exploding.
    # Without the Metal binary to compare to we can't reproduce the
    # exact magnitude, so we just check it's in a sensible band.
    if not (1e-8 <= mag <= 1e8):
        print(f"[FAIL] ||ext_accel|| = {mag:.3e} outside [1e-8, 1e8]")
        return 1

    print("\n[PASS] pair_forces_grid_fwd ran cleanly, output non-NaN and in range")
    return 0


if __name__ == "__main__":
    sys.exit(main())
