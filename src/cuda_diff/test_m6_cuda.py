"""Validate the CUDA M6 kernels against numpy references.

Run after build:
    src/cuda_diff/build.bat       (Windows)
    src/cuda_diff/build.sh        (Linux)
    python src/cuda_diff/test_m6_cuda.py

Mirrors src/metal_diff/test_dist.py — same numpy refs and tolerances —
but Windows-portable (uses tempfile and the .exe binary). Local-only.
"""
import math
import os
import platform
import subprocess
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()


# ─────────────────────────────  M6.0  ─────────────────────────────
def run_dist_as(active, static_p):
    n_active, n_static = len(active), len(static_p)
    a = os.path.join(TMP, "sib_cuda_a.bin")
    s = os.path.join(TMP, "sib_cuda_s.bin")
    o = os.path.join(TMP, "sib_cuda_o.bin")
    active.astype(np.float32).tofile(a)
    static_p.astype(np.float32).tofile(s)
    t0 = time.perf_counter()
    subprocess.run([BINARY, "dist_active_static",
                    str(n_active), str(n_static), a, s, o], check=True)
    wall = time.perf_counter() - t0
    return np.fromfile(o, dtype=np.float32).reshape(n_active, n_static), wall


def numpy_dist_as(active, static_p):
    diff = active[:, None, :] - static_p[None, :, :]
    return np.sum(diff * diff, axis=-1).astype(np.float32)


def test_dist_as(name, n_active, n_static, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * 5.0
    static_p = np.random.randn(n_static, 3).astype(np.float32) * 5.0
    cuda, wall = run_dist_as(active, static_p)
    expected = numpy_dist_as(active, static_p)
    err_max = float(np.abs(cuda - expected).max())
    print(f"=== {name} (n_active={n_active}, n_static={n_static}) ===")
    print(f"  Max error:  {err_max:.6e}")
    print(f"  CUDA wall:  {wall*1000:.1f} ms")
    ok = err_max < tolerance
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.3  ─────────────────────────────
def run_dist_aa(active):
    n_active = len(active)
    a = os.path.join(TMP, "sib_cuda_a.bin")
    o = os.path.join(TMP, "sib_cuda_aa.bin")
    active.astype(np.float32).tofile(a)
    t0 = time.perf_counter()
    subprocess.run([BINARY, "dist_active_active", str(n_active), a, o],
                   check=True)
    wall = time.perf_counter() - t0
    return np.fromfile(o, dtype=np.float32).reshape(n_active, n_active), wall


def numpy_dist_aa(active):
    diff = active[:, None, :] - active[None, :, :]
    return np.sum(diff * diff, axis=-1).astype(np.float32)


def test_dist_aa(name, n_active, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * 5.0
    cuda, wall = run_dist_aa(active)
    expected = numpy_dist_aa(active)
    err_max = float(np.abs(cuda - expected).max())
    diag_max = float(np.abs(np.diag(cuda)).max())
    print(f"=== {name} (n_active={n_active}) ===")
    print(f"  Diagonal max:  {diag_max:.3e}  (should be ~0)")
    print(f"  Max error:     {err_max:.3e}")
    print(f"  CUDA wall:     {wall*1000:.1f} ms")
    ok = err_max < tolerance
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.2  ─────────────────────────────
def run_rowsum(W, mass):
    n_rows, n_cols = W.shape
    p_W = os.path.join(TMP, "sib_cuda_W.bin")
    p_o = os.path.join(TMP, "sib_cuda_density.bin")
    W.astype(np.float32).tofile(p_W)
    t0 = time.perf_counter()
    subprocess.run([BINARY, "rowsum_density",
                    str(n_rows), str(n_cols), str(mass), p_W, p_o], check=True)
    wall = time.perf_counter() - t0
    return np.fromfile(p_o, dtype=np.float32), wall


def numpy_rowsum(W, mass):
    return (mass * W.sum(axis=1)).astype(np.float32)


def test_rowsum(name, W, mass, tolerance=1e-2):
    cuda, wall = run_rowsum(W, mass)
    expected = numpy_rowsum(W, mass)
    rel = np.abs(cuda - expected) / (np.abs(expected) + 1e-9)
    err_max = float(rel.max())
    print(f"=== {name} (n_rows={W.shape[0]}, n_cols={W.shape[1]}) ===")
    print(f"  Mean density:   {expected.mean():.4e}")
    print(f"  Max rel error:  {err_max:.6e}")
    print(f"  CUDA wall:      {wall*1000:.1f} ms")
    ok = err_max < tolerance
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.4  ─────────────────────────────
def numpy_density_grad(active, static_p, h, mass, rho_rest):
    spiky_const = -45.0 / (math.pi * h ** 6)
    diff_aa = active[:, None, :] - active[None, :, :]
    r_aa = np.linalg.norm(diff_aa, axis=-1)
    mask_aa = (r_aa < h) & (r_aa > 1e-9)
    coef_aa = np.where(mask_aa,
                       spiky_const * (h - r_aa) ** 2 / (r_aa + 1e-7), 0.0)
    grad_aa = (mass * coef_aa[..., None] * diff_aa).sum(axis=1)
    diff_as = active[:, None, :] - static_p[None, :, :]
    r_as = np.linalg.norm(diff_as, axis=-1)
    mask_as = r_as < h
    coef_as = np.where(mask_as,
                       spiky_const * (h - r_as) ** 2 / (r_as + 1e-7), 0.0)
    grad_as = (mass * coef_as[..., None] * diff_as).sum(axis=1)
    return ((grad_aa + grad_as) / rho_rest).astype(np.float32)


def run_density_grad(active, static_p, r2_aa, r2_as, h, mass, rho_rest):
    n_active, n_static = len(active), len(static_p)
    p_a = os.path.join(TMP, "sib_cuda_grad_a.bin")
    p_s = os.path.join(TMP, "sib_cuda_grad_s.bin")
    p_aa = os.path.join(TMP, "sib_cuda_grad_r2aa.bin")
    p_as = os.path.join(TMP, "sib_cuda_grad_r2as.bin")
    p_out = os.path.join(TMP, "sib_cuda_grad_out.bin")
    p_den = os.path.join(TMP, "sib_cuda_grad_denom.bin")
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    r2_aa.astype(np.float32).tofile(p_aa)
    r2_as.astype(np.float32).tofile(p_as)
    t0 = time.perf_counter()
    subprocess.run([BINARY, "density_constraint_grad",
                    str(n_active), str(n_static), str(h), str(mass),
                    str(rho_rest), p_a, p_s, p_aa, p_as, p_out, p_den],
                   check=True)
    wall = time.perf_counter() - t0
    grad = np.fromfile(p_out, dtype=np.float32).reshape(n_active, 3)
    denom = np.fromfile(p_den, dtype=np.float32)
    return grad, denom, wall


def test_density_grad(name, n_active, n_static, h, mass, rho_rest,
                      tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(n_static, 3).astype(np.float32) * (h * 0.7)
    diff_aa = active[:, None, :] - active[None, :, :]
    r2_aa = np.sum(diff_aa * diff_aa, axis=-1).astype(np.float32)
    diff_as = active[:, None, :] - static_p[None, :, :]
    r2_as = np.sum(diff_as * diff_as, axis=-1).astype(np.float32)

    cuda_grad, _denom, wall = run_density_grad(
        active, static_p, r2_aa, r2_as, h, mass, rho_rest)
    expected = numpy_density_grad(active, static_p, h, mass, rho_rest)

    err_abs = np.abs(cuda_grad - expected)
    err_max = float(err_abs.max())
    grad_mag = float(np.linalg.norm(expected, axis=1).mean())
    err_rel = err_max / (grad_mag + 1e-9)
    print(f"=== {name} (n_active={n_active}, n_static={n_static}) ===")
    print(f"  Mean |grad C|:   {grad_mag:.4e}")
    print(f"  Max abs error:   {err_max:.4e}")
    print(f"  Max rel error:   {err_rel:.4e}")
    print(f"  CUDA wall:       {wall*1000:.1f} ms")
    ok = err_rel < tolerance
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    all_pass = True
    h = 5.0

    # M6.0 active×static distance
    all_pass &= test_dist_as("M6.0 small", 100, 200)
    all_pass &= test_dist_as("M6.0 demo1-realistic", 343, 17498)

    # M6.3 active×active distance
    all_pass &= test_dist_aa("M6.3 small", 50)
    all_pass &= test_dist_aa("M6.3 demo1-realistic", 343)

    # M6.2 rowsum (over wpoly6 output for realism)
    from test_wpoly6_cuda import run_wpoly6, numpy_wpoly6
    np.random.seed(0)
    active = np.random.randn(343, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(17498, 3).astype(np.float32) * (h * 0.7)
    diff = active[:, None, :] - static_p[None, :, :]
    r2 = np.sum(diff * diff, axis=-1).astype(np.float32)
    W = numpy_wpoly6(r2, h)
    all_pass &= test_rowsum("M6.2 demo1-realistic", W, mass=1.0)

    # M6.4 density-constraint gradient
    all_pass &= test_density_grad("M6.4 small (20+30)",  20,  30, h, 1.0, 10.0)
    all_pass &= test_density_grad("M6.4 demo1-realistic (343+17498)",
                                  343, 17498, h, 1.0, 10.0)

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
