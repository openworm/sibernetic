"""Validate spring_bonds_force_backward against finite differences.

Setup: 4 active particles connected by 3 bonds (a chain).
Forward: spring_bonds_force → ext_accel
Loss: L = sum(ext_accel)
Backward: spring_bonds_force_backward → grad_pos
Validate: max abs (kernel - num) / max(|num|, 1) < 1e-2
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def write_bin(path, arr, dtype):
    arr.astype(dtype).tofile(path)


def run_fwd(n_active, n_bonds, spring_K, paths):
    # NB: C op writes output to /tmp/spring_ext_accel.bin (fixed path).
    cmd = [BIN, "spring_bonds_fwd",
           str(n_active), str(n_bonds), str(spring_K),
           paths['pos'], paths['bond_ij'], paths['bond_rest']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr, file=sys.stderr)
        sys.exit(1)
    return np.fromfile("/tmp/spring_ext_accel.bin",
                       dtype=np.float32).reshape(-1, 3)


def run_bwd(n_active, n_bonds, spring_K, paths, grad_ext):
    write_bin("/tmp/spring_grad_ext.bin", grad_ext, np.float32)
    cmd = [BIN, "spring_bonds_bwd",
           str(n_active), str(n_bonds), str(spring_K),
           paths['pos'], paths['bond_ij'], paths['bond_rest'],
           "/tmp/spring_grad_ext.bin"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr, file=sys.stderr)
        sys.exit(1)
    return np.fromfile("/tmp/spring_grad_pos.bin",
                       dtype=np.float32).reshape(-1, 3)


def main():
    print("Test: spring_bonds_force_backward via finite differences\n")

    np.random.seed(7)
    n_active = 4
    n_bonds = 3
    spring_K = 100.0

    # Chain: 0-1, 1-2, 2-3, with rest length 1.0 each.
    pos = np.array([[0.0, 0.0, 0.0],
                    [1.05, 0.0, 0.0],   # slightly stretched
                    [2.10, 0.0, 0.0],
                    [3.15, 0.0, 0.0]], dtype=np.float32)

    bond_ij = np.array([(0, 1), (1, 2), (2, 3)], dtype=[('i', 'i4'), ('j', 'i4')])
    bond_rest = np.array([1.0, 1.0, 1.0], dtype=np.float32)

    paths = {'pos': '/tmp/spring_pos.bin',
             'bond_ij': '/tmp/spring_bij.bin',
             'bond_rest': '/tmp/spring_brl.bin'}
    write_bin(paths['pos'], pos, np.float32)
    bond_ij.tofile(paths['bond_ij'])
    write_bin(paths['bond_rest'], bond_rest, np.float32)

    ext0 = run_fwd(n_active, n_bonds, spring_K, paths)
    print(f"Initial ext_accel:\n{ext0}")
    print(f"  magnitude: {np.linalg.norm(ext0):.3e}")
    assert not np.isnan(ext0).any(), "NaN in forward"

    # Asymmetric grad_ext_accel: per-particle weight breaks Newton's-3rd-law
    # cancellation that makes loss=sum(ext_accel) trivially zero.
    grad_ext = np.zeros_like(ext0)
    for i in range(n_active):
        grad_ext[i] = (i + 1, 0.5 * (i + 1), 0.25 * (i + 1))
    gp_kernel = run_bwd(n_active, n_bonds, spring_K, paths, grad_ext)

    # Numerical: L = sum_i <grad_ext[i], ext_accel[i]>
    eps = 1e-3
    gp_num = np.zeros_like(pos)
    for i in range(n_active):
        for d in range(3):
            pos_p = pos.copy(); pos_p[i, d] += eps
            write_bin(paths['pos'], pos_p, np.float32)
            ext_p = run_fwd(n_active, n_bonds, spring_K, paths)

            pos_p[i, d] -= 2 * eps
            write_bin(paths['pos'], pos_p, np.float32)
            ext_m = run_fwd(n_active, n_bonds, spring_K, paths)

            L_p = (grad_ext * ext_p).sum()
            L_m = (grad_ext * ext_m).sum()
            gp_num[i, d] = (L_p - L_m) / (2 * eps)
            write_bin(paths['pos'], pos, np.float32)  # restore

    print(f"\nKernel grad_pos:\n{gp_kernel}")
    print(f"Numer  grad_pos:\n{gp_num}")

    denom = max(np.abs(gp_num).max(), 1e-30)
    err = (np.abs(gp_kernel - gp_num) / denom).max()
    print(f"\nMax rel err: {err:.3e}")

    if err < 1e-2:
        print("\n[PASS] spring_bonds_force_backward matches finite-diff")
        return 0
    print("\n[FAIL]")
    return 1


if __name__ == "__main__":
    sys.exit(main())
