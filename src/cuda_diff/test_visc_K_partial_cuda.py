"""FD validation of CUDA visc_K_partial kernel.

visc_K_partial computes the per-particle dot of (d visc_accel_i /
d visc_pair_coef) with upstream grad_ext_accel[i]. Summing those scalars
gives the total dL/d(visc_pair_coef) for use in xpbd_full_bwd.

Procedure (mirrors test_pair_forces_grid_cuda.py scaffolding):
  1. Run pair_forces_grid_fwd at baseline visc_pair_coef = 1e-4.
  2. Pick ones-vector pullback: grad_ext = ones_like(ext_accel), so L = sum(ext).
  3. Run visc_K_partial, read back per-particle array, sum to get analytic.
  4. Numerical FD: two forward calls at visc_pair_coef +/- eps,
     numerical = (sum(ext+) - sum(ext-)) / (2*eps).
  5. Threshold rel_err < 1e-2 (same as pair_forces_grid bwd test).

Windows-portable: tempfile.gettempdir() instead of /tmp; grid_origin
passed explicitly on the CLI. Local-only.
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
    test_pair_forces_grid_cuda.py::make_grid)."""
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


def run_fwd(args, paths, visc_pair_coef):
    cmd = [BINARY, "pair_forces_grid_fwd",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']),
           str(args['sim_scale']), str(visc_pair_coef),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(int(args['grid_dim'][0])),
           str(int(args['grid_dim'][1])),
           str(int(args['grid_dim'][2])),
           paths['grid_origin'], paths['ext_accel']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("FWD STDERR:", r.stderr, file=sys.stderr)
        print("FWD STDOUT:", r.stdout, file=sys.stderr)
        sys.exit(1)
    return np.fromfile(paths['ext_accel'], dtype=np.float32).reshape(-1, 3)


def run_visc_K_partial(args, paths, grad_ext):
    write_bin(paths['grad_ext'], grad_ext)
    cmd = [BINARY, "visc_K_partial",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']), str(args['sim_scale']),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(int(args['grid_dim'][0])),
           str(int(args['grid_dim'][1])),
           str(int(args['grid_dim'][2])),
           paths['grid_origin'], paths['grad_ext'],
           paths['per_particle']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("VISC_K STDERR:", r.stderr, file=sys.stderr)
        print("VISC_K STDOUT:", r.stdout, file=sys.stderr)
        sys.exit(1)
    return np.fromfile(paths['per_particle'], dtype=np.float32)


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    print("Test: visc_K_partial (FD validation)\n")

    sim_scale = 7.4e-6
    h = 3.34
    mass = 2.0e-12
    visc_pair_coef_0 = 1e-4

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

    density = np.array([4.05e-13, 4.05e-13], dtype=np.float32)

    paths = {
        'pos_a':         os.path.join(TMP, 'cuda_viscK_pos_a.bin'),
        'vel_a':         os.path.join(TMP, 'cuda_viscK_vel_a.bin'),
        'sorted_static': os.path.join(TMP, 'cuda_viscK_sorted_static.bin'),
        'cell_start':    os.path.join(TMP, 'cuda_viscK_cell_start.bin'),
        'density':       os.path.join(TMP, 'cuda_viscK_density.bin'),
        'grid_origin':   os.path.join(TMP, 'cuda_viscK_grid_origin.bin'),
        'ext_accel':     os.path.join(TMP, 'cuda_viscK_ext_accel.bin'),
        'grad_ext':      os.path.join(TMP, 'cuda_viscK_grad_ext.bin'),
        'per_particle':  os.path.join(TMP, 'cuda_viscK_per_particle.bin'),
    }
    write_bin(paths['pos_a'], pos_a)
    write_bin(paths['vel_a'], vel_a)
    write_bin(paths['sorted_static'], sorted_static)
    cell_start.astype(np.int32).tofile(paths['cell_start'])
    write_bin(paths['density'], density)
    write_bin(paths['grid_origin'], grid_origin)

    args = {
        'n_active': n_active, 'n_static': n_static,
        'h': h, 'mass': mass, 'sim_scale': sim_scale,
        'grid_dim': grid_dim,
    }

    # ── Forward sanity at baseline visc_pair_coef ──
    ext0 = run_fwd(args, paths, visc_pair_coef_0)
    print(f"Baseline ext_accel:\n{ext0}")
    print(f"  ||ext_accel|| = {np.linalg.norm(ext0):.3e}")
    if np.isnan(ext0).any() or np.isinf(ext0).any():
        print("[FAIL] NaN/Inf in forward output")
        return 1

    # ── Analytic: visc_K_partial with grad_ext = ones, then sum ──
    grad_ext = np.ones_like(ext0)
    per_particle = run_visc_K_partial(args, paths, grad_ext)
    print(f"\nPer-particle:\n{per_particle}")
    analytic = float(per_particle.sum())
    print(f"Analytic dL/d(visc_pair_coef) = {analytic:.6e}")

    # ── Numerical FD ──
    eps = 1e-7  # baseline is 1e-4, relative perturbation = 1e-3
    ext_p = run_fwd(args, paths, visc_pair_coef_0 + eps)
    ext_m = run_fwd(args, paths, visc_pair_coef_0 - eps)
    numerical = float((ext_p.sum() - ext_m.sum()) / (2.0 * eps))
    print(f"Numerical dL/d(visc_pair_coef) = {numerical:.6e}")

    rel_err = abs(analytic - numerical) / max(abs(numerical), 1e-30)
    print(f"\nrel_err = {rel_err:.3e}")

    if rel_err < 1e-2:
        print("[PASS] visc_K_partial matches finite-diff")
        return 0
    print("[FAIL] visc_K_partial does not match finite-diff")
    return 1


if __name__ == "__main__":
    sys.exit(main())
