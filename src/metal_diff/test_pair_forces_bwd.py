"""Validate pair_forces_grid_backward against finite differences.

Setup:
- 2 active particles + 4 boundary particles in a tiny box.
- Forward: pair_forces_grid → ext_accel (3-vec per active).
- Loss: L = sum(ext_accel) — i.e. ones-vector pullback so dL/d(ext_accel) ≡ 1.
- Backward: pair_forces_grid_backward → grad_pos, grad_vel.
- Numerical: perturb each component of pos and vel by ±ε and finite-diff L.

Pass criterion: max abs (kernel_grad - num_grad) / max(|num_grad|, 1) < 1e-2.
The threshold is loose because:
- Pair forces have steep ∇W terms near r=0 / r=h boundaries.
- fp32 host formula for surf_amp has ~6 digits of mantissa.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def write_bin(path, arr, dtype=np.float32):
    arr.astype(dtype).tofile(path)


def make_grid(pos_static, h):
    """Mini grid build for a few static particles."""
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
    # NB: C op reads grid_origin from /tmp/pair_grid_origin.bin (fixed
    # path) — we wrote it there during test setup.
    cmd = [BIN, "pair_forces_fwd",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']),
           str(args['sim_scale']), str(args['visc_pair_coef']),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(args['grid_dim'][0]), str(args['grid_dim'][1]),
           str(args['grid_dim'][2])]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr, file=sys.stderr)
        sys.exit(1)
    return np.fromfile("/tmp/pair_ext_accel.bin",
                       dtype=np.float32).reshape(-1, 3)


def run_bwd(args, paths, grad_ext):
    write_bin("/tmp/pair_grad_ext.bin", grad_ext)
    # 16th arg slot is unused (placeholder for grid_origin path the C
    # op no longer reads) — pass any string. 17th is grad_ext path.
    cmd = [BIN, "pair_forces_bwd",
           str(args['n_active']), str(args['n_static']),
           str(args['h']), str(args['mass']),
           str(args['sim_scale']), str(args['visc_pair_coef']),
           paths['pos_a'], paths['vel_a'], paths['sorted_static'],
           paths['cell_start'], paths['density'],
           str(args['grid_dim'][0]), str(args['grid_dim'][1]),
           str(args['grid_dim'][2]),
           "/tmp/pair_grad_ext.bin"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr, file=sys.stderr)
        sys.exit(1)
    gp = np.fromfile("/tmp/pair_grad_pos.bin",
                     dtype=np.float32).reshape(-1, 3)
    gv = np.fromfile("/tmp/pair_grad_vel.bin",
                     dtype=np.float32).reshape(-1, 3)
    return gp, gv


def main():
    print("Test: pair_forces_grid_backward via finite differences\n")

    np.random.seed(42)
    sim_scale = 7.4e-6
    h = 3.34
    mass = 2.0e-12
    visc_pair_coef = 1e-4

    # 2 active particles, neighbors of each other and of boundary.
    n_active = 2
    pos_a = np.array([[1.0, 1.0, 1.0],
                      [1.0, 1.0, 2.5]], dtype=np.float32)
    vel_a = np.array([[0.0, -0.5, 0.0],
                      [0.0, +0.3, 0.1]], dtype=np.float32)

    # 4 boundary particles surrounding them.
    pos_static = np.array([[0.0, 0.0, 0.0],
                           [2.5, 0.0, 0.0],
                           [0.0, 0.0, 3.0],
                           [2.5, 0.0, 3.0]], dtype=np.float32)
    n_static = len(pos_static)

    sorted_static, cell_start, grid_dim, grid_origin = make_grid(
        pos_static, h)

    # Density: just hand-pick reasonable values
    density = np.array([4.05e-13, 4.05e-13], dtype=np.float32)

    # Write all files
    paths = {
        'pos_a':         '/tmp/pair_pos_a.bin',
        'vel_a':         '/tmp/pair_vel_a.bin',
        'sorted_static': '/tmp/pair_sorted_static.bin',
        'cell_start':    '/tmp/pair_cell_start.bin',
        'density':       '/tmp/pair_density.bin',
        'grid_origin':   '/tmp/pair_grid_origin.bin',
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

    # Forward sanity
    ext0 = run_fwd(args, paths)
    print(f"Initial ext_accel: {ext0}")
    print(f"  magnitude: {np.linalg.norm(ext0):.3e}")
    assert not np.isnan(ext0).any(), "NaN in forward"

    # Backward with grad_ext = ones
    grad_ext = np.ones_like(ext0)
    gp_kernel, gv_kernel = run_bwd(args, paths, grad_ext)

    # Numerical: L = sum(ext_accel)
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

    # Relative error
    def rel(k, n):
        denom = max(np.abs(n).max(), 1e-30)
        return (np.abs(k - n) / denom).max()

    err_pos = rel(gp_kernel, gp_num)
    err_vel = rel(gv_kernel, gv_num)
    print(f"\nMax rel err pos: {err_pos:.3e}")
    print(f"Max rel err vel: {err_vel:.3e}")

    pos_pass = err_pos < 1e-2
    vel_pass = err_vel < 1e-2
    print(f"\nPos backward:  {'[PASS]' if pos_pass else '[FAIL]'}")
    print(f"Vel backward:  {'[PASS]' if vel_pass else '[FAIL]'}")

    if pos_pass and vel_pass:
        print("\n[OVERALL PASS] pair_forces_grid_backward matches finite-diff")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
