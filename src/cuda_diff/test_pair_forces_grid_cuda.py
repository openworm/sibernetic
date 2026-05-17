"""FD validation of CUDA pair_forces_grid forward + backward kernels.

Mirrors src/metal_diff/test_pair_forces_bwd.py end-to-end:
- Forward call as the sanity check.
- Backward call with grad_ext = ones (L = sum(ext_accel) so dL/d(ext) == 1).
- Numerical: perturb each (i, d) component of pos and vel by +/-eps and FD.
- Compare kernel grads vs numerical grads. Loose threshold (1e-2) because
  pair forces have steep gradients near r=0 and r=h, and fp32 mantissa
  only has ~6 digits.

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


def run_fwd(args, paths):
    cmd = [BINARY, "pair_forces_grid_fwd",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']),
           str(args['sim_scale']), str(args['visc_pair_coef']),
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


def run_bwd(args, paths, grad_ext):
    write_bin(paths['grad_ext'], grad_ext)
    cmd = [BINARY, "pair_forces_grid_bwd",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']),
           str(args['sim_scale']), str(args['visc_pair_coef']),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(int(args['grid_dim'][0])),
           str(int(args['grid_dim'][1])),
           str(int(args['grid_dim'][2])),
           paths['grid_origin'], paths['grad_ext'],
           paths['grad_pos'], paths['grad_vel']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("BWD STDERR:", r.stderr, file=sys.stderr)
        print("BWD STDOUT:", r.stdout, file=sys.stderr)
        sys.exit(1)
    gp = np.fromfile(paths['grad_pos'], dtype=np.float32).reshape(-1, 3)
    gv = np.fromfile(paths['grad_vel'], dtype=np.float32).reshape(-1, 3)
    return gp, gv


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    print("Test: pair_forces_grid forward + backward (FD validation)\n")

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

    density = np.array([4.05e-13, 4.05e-13], dtype=np.float32)

    paths = {
        'pos_a':         os.path.join(TMP, 'cuda_pair_pos_a.bin'),
        'vel_a':         os.path.join(TMP, 'cuda_pair_vel_a.bin'),
        'sorted_static': os.path.join(TMP, 'cuda_pair_sorted_static.bin'),
        'cell_start':    os.path.join(TMP, 'cuda_pair_cell_start.bin'),
        'density':       os.path.join(TMP, 'cuda_pair_density.bin'),
        'grid_origin':   os.path.join(TMP, 'cuda_pair_grid_origin.bin'),
        'ext_accel':     os.path.join(TMP, 'cuda_pair_ext_accel.bin'),
        'grad_ext':      os.path.join(TMP, 'cuda_pair_grad_ext.bin'),
        'grad_pos':      os.path.join(TMP, 'cuda_pair_grad_pos.bin'),
        'grad_vel':      os.path.join(TMP, 'cuda_pair_grad_vel.bin'),
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
        'visc_pair_coef': visc_pair_coef,
        'grid_dim': grid_dim,
    }

    # ── Forward sanity ──
    ext0 = run_fwd(args, paths)
    print(f"Initial ext_accel:\n{ext0}")
    print(f"  ||ext_accel|| = {np.linalg.norm(ext0):.3e}")
    if np.isnan(ext0).any() or np.isinf(ext0).any():
        print("[FAIL] NaN/Inf in forward output")
        return 1

    # ── Backward with grad_ext = ones (so L = sum(ext_accel)) ──
    grad_ext = np.ones_like(ext0)
    gp_kernel, gv_kernel = run_bwd(args, paths, grad_ext)

    # ── Numerical FD ──
    eps = 1e-3
    gp_num = np.zeros_like(pos_a)
    gv_num = np.zeros_like(vel_a)

    for i in range(n_active):
        for d in range(3):
            # pos perturbation
            pos_pert = pos_a.copy()
            pos_pert[i, d] += eps
            write_bin(paths['pos_a'], pos_pert)
            ext_p = run_fwd(args, paths)
            pos_pert[i, d] -= 2 * eps
            write_bin(paths['pos_a'], pos_pert)
            ext_m = run_fwd(args, paths)
            gp_num[i, d] = (ext_p - ext_m).sum() / (2 * eps)
            write_bin(paths['pos_a'], pos_a)  # restore

            # vel perturbation
            vel_pert = vel_a.copy()
            vel_pert[i, d] += eps
            write_bin(paths['vel_a'], vel_pert)
            ext_p = run_fwd(args, paths)
            vel_pert[i, d] -= 2 * eps
            write_bin(paths['vel_a'], vel_pert)
            ext_m = run_fwd(args, paths)
            gv_num[i, d] = (ext_p - ext_m).sum() / (2 * eps)
            write_bin(paths['vel_a'], vel_a)  # restore

    print(f"\nKernel grad_pos:\n{gp_kernel}")
    print(f"Numer  grad_pos:\n{gp_num}")
    print(f"\nKernel grad_vel:\n{gv_kernel}")
    print(f"Numer  grad_vel:\n{gv_num}")

    def rel(k, n):
        denom = max(np.abs(n).max(), 1.0)
        return float((np.abs(k - n) / denom).max())

    err_pos = rel(gp_kernel, gp_num)
    err_vel = rel(gv_kernel, gv_num)
    print(f"\nMax rel err pos: {err_pos:.3e}")
    print(f"Max rel err vel: {err_vel:.3e}")

    pos_pass = err_pos < 1e-2
    vel_pass = err_vel < 1e-2
    print(f"\nPos backward:  {'[PASS]' if pos_pass else '[FAIL]'}")
    print(f"Vel backward:  {'[PASS]' if vel_pass else '[FAIL]'}")

    if pos_pass and vel_pass:
        print("\n[OVERALL PASS] pair_forces_grid backward matches finite-diff")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
