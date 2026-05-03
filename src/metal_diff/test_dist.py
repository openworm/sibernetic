"""Validate hand-written Metal kernels against numpy references.

Run:
    .venv/bin/python src/metal_diff/test_dist.py

Expects sib_metal binary built via src/metal_diff/build.sh.
"""
import math
import os
import subprocess
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def run_metal(active: np.ndarray, static_p: np.ndarray):
    n_active = len(active)
    n_static = len(static_p)
    active_path = "/tmp/sib_metal_active.bin"
    static_path = "/tmp/sib_metal_static.bin"
    out_path = "/tmp/sib_metal_out.bin"

    active.astype(np.float32).tofile(active_path)
    static_p.astype(np.float32).tofile(static_path)

    t0 = time.perf_counter()
    subprocess.run(
        [BINARY, "dist_active_static",
         str(n_active), str(n_static),
         active_path, static_path, out_path],
        check=True,
    )
    wall = time.perf_counter() - t0

    result = np.fromfile(out_path, dtype=np.float32).reshape(n_active, n_static)
    return result, wall


def numpy_reference(active: np.ndarray, static_p: np.ndarray):
    diff = active[:, None, :] - static_p[None, :, :]
    return np.sum(diff * diff, axis=-1)


def test_case(name, n_active, n_static, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * 5.0
    static_p = np.random.randn(n_static, 3).astype(np.float32) * 5.0

    metal_result, metal_wall = run_metal(active, static_p)
    t0 = time.perf_counter()
    expected = numpy_reference(active, static_p).astype(np.float32)
    numpy_wall = time.perf_counter() - t0

    err_max = float(np.abs(metal_result - expected).max())
    err_mean = float(np.abs(metal_result - expected).mean())
    matrix_mb = n_active * n_static * 4 / 1024 / 1024

    print(f"=== {name} (n_active={n_active}, n_static={n_static}) ===")
    print(f"  Matrix size:        {matrix_mb:.1f} MB")
    print(f"  Max error:          {err_max:.6e}")
    print(f"  Mean error:         {err_mean:.6e}")
    print(f"  Metal wall (incl IO): {metal_wall*1000:.1f} ms")
    print(f"  numpy wall:           {numpy_wall*1000:.1f} ms")
    print(f"  Speedup vs numpy:     {numpy_wall/metal_wall:.2f}x")
    passed = err_max < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed


def run_wpoly6(r2: np.ndarray, h: float):
    inout_path = "/tmp/sib_metal_wpoly6.bin"
    r2.astype(np.float32).tofile(inout_path)
    t0 = time.perf_counter()
    subprocess.run(
        [BINARY, "wpoly6_inplace", str(r2.size), str(h), inout_path],
        check=True,
    )
    wall = time.perf_counter() - t0
    out = np.fromfile(inout_path, dtype=np.float32).reshape(r2.shape)
    return out, wall


def numpy_wpoly6(r2: np.ndarray, h: float):
    h2 = h * h
    poly6_const = 315.0 / (64.0 * math.pi * h ** 9)
    diff = h2 - r2
    out = poly6_const * diff ** 3
    out = np.where(r2 >= h2, 0.0, out)
    return out.astype(np.float32)


def run_density(W: np.ndarray, mass: float):
    n_rows, n_cols = W.shape
    w_path = "/tmp/sib_metal_W.bin"
    out_path = "/tmp/sib_metal_density.bin"
    W.astype(np.float32).tofile(w_path)
    t0 = time.perf_counter()
    subprocess.run(
        [BINARY, "rowsum_density",
         str(n_rows), str(n_cols), str(mass), w_path, out_path],
        check=True,
    )
    wall = time.perf_counter() - t0
    out = np.fromfile(out_path, dtype=np.float32)
    return out, wall


def numpy_density(W: np.ndarray, mass: float):
    return (mass * W.sum(axis=1)).astype(np.float32)


def test_wpoly6(name, n_active, n_static, h, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(n_static, 3).astype(np.float32) * (h * 0.7)
    diff = active[:, None, :] - static_p[None, :, :]
    r2 = np.sum(diff * diff, axis=-1).astype(np.float32)

    metal_W, metal_wall = run_wpoly6(r2, h)
    expected_W = numpy_wpoly6(r2, h)

    err_max = float(np.abs(metal_W - expected_W).max())
    err_mean = float(np.abs(metal_W - expected_W).mean())
    nonzero_frac = float((metal_W > 0).sum()) / metal_W.size

    print(f"=== {name} (n={r2.size}, h={h}) ===")
    print(f"  Nonzero W fraction: {nonzero_frac:.1%}")
    print(f"  Max error:          {err_max:.6e}")
    print(f"  Mean error:         {err_mean:.6e}")
    print(f"  Metal wall (incl IO): {metal_wall*1000:.1f} ms")
    passed = err_max < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed, metal_W


def test_density(name, W: np.ndarray, mass: float, tolerance=1e-2):
    metal_rho, metal_wall = run_density(W, mass)
    expected_rho = numpy_density(W, mass)

    # density values can be large; use relative-error tolerance for a fair
    # PASS threshold (rather than absolute).
    rel = np.abs(metal_rho - expected_rho) / (np.abs(expected_rho) + 1e-9)
    err_max = float(rel.max())
    err_mean = float(rel.mean())

    print(f"=== {name} (n_rows={W.shape[0]}, n_cols={W.shape[1]}) ===")
    print(f"  Mean density: {expected_rho.mean():.4e}")
    print(f"  Max rel error:  {err_max:.6e}")
    print(f"  Mean rel error: {err_mean:.6e}")
    print(f"  Metal wall (incl IO): {metal_wall*1000:.1f} ms")
    passed = err_max < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed


def run_dist_aa(active: np.ndarray):
    n_active = len(active)
    active_path = "/tmp/sib_metal_active.bin"
    out_path = "/tmp/sib_metal_dist_aa.bin"
    active.astype(np.float32).tofile(active_path)
    t0 = time.perf_counter()
    subprocess.run(
        [BINARY, "dist_active_active",
         str(n_active), active_path, out_path],
        check=True,
    )
    wall = time.perf_counter() - t0
    out = np.fromfile(out_path, dtype=np.float32).reshape(n_active, n_active)
    return out, wall


def numpy_dist_aa(active: np.ndarray):
    diff = active[:, None, :] - active[None, :, :]
    return np.sum(diff * diff, axis=-1).astype(np.float32)


def test_dist_aa(name, n_active, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * 5.0
    metal_d, metal_wall = run_dist_aa(active)
    expected = numpy_dist_aa(active)

    err_max = float(np.abs(metal_d - expected).max())
    err_mean = float(np.abs(metal_d - expected).mean())
    diag_max = float(np.abs(np.diag(metal_d)).max())
    print(f"=== {name} (n_active={n_active}) ===")
    print(f"  Matrix size:   {n_active * n_active * 4 / 1024:.1f} KB")
    print(f"  Diagonal max:  {diag_max:.3e}  (should be ~0; self-distance)")
    print(f"  Max error:     {err_max:.3e}")
    print(f"  Metal wall:    {metal_wall*1000:.1f} ms")
    passed = err_max < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed


def numpy_density_grad(active, static_p, h, mass, rho_rest):
    """Reference implementation of XPBD density constraint gradient.

    ∇C_i = (1/ρ_rest) Σ_j m_j ∇W_spiky(p_i - p_j)
    where the sum is over active neighbors (j != i) and static neighbors.
    """
    spiky_const = -45.0 / (math.pi * h ** 6)
    n_active = len(active)

    # active-active contributions
    diff_aa = active[:, None, :] - active[None, :, :]
    r_aa = np.linalg.norm(diff_aa, axis=-1)
    mask_aa = (r_aa < h) & (r_aa > 1e-9)
    coef_aa = np.where(mask_aa,
                       spiky_const * (h - r_aa) ** 2 / (r_aa + 1e-7),
                       0.0)
    grad_aa = (mass * coef_aa[..., None] * diff_aa).sum(axis=1)

    # active-static contributions
    diff_as = active[:, None, :] - static_p[None, :, :]
    r_as = np.linalg.norm(diff_as, axis=-1)
    mask_as = r_as < h
    coef_as = np.where(mask_as,
                       spiky_const * (h - r_as) ** 2 / (r_as + 1e-7),
                       0.0)
    grad_as = (mass * coef_as[..., None] * diff_as).sum(axis=1)

    return ((grad_aa + grad_as) / rho_rest).astype(np.float32)


def run_density_grad(active, static_p, r2_aa, r2_as, h, mass, rho_rest):
    n_active = len(active)
    n_static = len(static_p)
    p_a = "/tmp/sm_grad_a.bin"
    p_s = "/tmp/sm_grad_s.bin"
    p_aa = "/tmp/sm_grad_r2aa.bin"
    p_as = "/tmp/sm_grad_r2as.bin"
    p_out = "/tmp/sm_grad_out.bin"
    p_denom = "/tmp/sm_grad_denom.bin"
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    r2_aa.astype(np.float32).tofile(p_aa)
    r2_as.astype(np.float32).tofile(p_as)
    t0 = time.perf_counter()
    subprocess.run(
        [BINARY, "density_constraint_grad",
         str(n_active), str(n_static), str(h), str(mass), str(rho_rest),
         p_a, p_s, p_aa, p_as, p_out, p_denom],
        check=True,
    )
    wall = time.perf_counter() - t0
    out = np.fromfile(p_out, dtype=np.float32).reshape(n_active, 3)
    denom = np.fromfile(p_denom, dtype=np.float32)
    return out, denom, wall


def test_density_grad(name, n_active, n_static, h, mass, rho_rest,
                      tolerance=1e-3):
    np.random.seed(0)
    # Particles within smoothing radius of each other so kernel actually fires.
    active = np.random.randn(n_active, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(n_static, 3).astype(np.float32) * (h * 0.7)

    diff_aa = active[:, None, :] - active[None, :, :]
    r2_aa = np.sum(diff_aa * diff_aa, axis=-1).astype(np.float32)
    diff_as = active[:, None, :] - static_p[None, :, :]
    r2_as = np.sum(diff_as * diff_as, axis=-1).astype(np.float32)

    metal_grad, metal_denom, metal_wall = run_density_grad(
        active, static_p, r2_aa, r2_as, h, mass, rho_rest
    )
    expected = numpy_density_grad(active, static_p, h, mass, rho_rest)

    # Component-wise relative error.
    err_abs = np.abs(metal_grad - expected)
    err_max = float(err_abs.max())
    grad_mag = float(np.linalg.norm(expected, axis=1).mean())
    err_rel = err_max / (grad_mag + 1e-9)

    print(f"=== {name} (n_active={n_active}, n_static={n_static}) ===")
    print(f"  Mean |∇C|:      {grad_mag:.4e}")
    print(f"  Max abs error:  {err_max:.4e}")
    print(f"  Max rel error:  {err_rel:.4e}")
    print(f"  Metal wall:     {metal_wall*1000:.1f} ms")
    passed = err_rel < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    all_pass = True
    # M6.0
    all_pass &= test_case("M6.0 small correctness", 100, 200)
    all_pass &= test_case("M6.0 demo1-realistic", 343, 17498)

    # M6.1: Wpoly6 elementwise
    h = 5.0
    p1, W_small = test_wpoly6("M6.1 small (50×80)", 50, 80, h)
    p2, W_big = test_wpoly6("M6.1 demo1-realistic (343×17498)", 343, 17498, h)
    all_pass &= p1 and p2

    # M6.2: density via row-sum (uses W_big from M6.1)
    p3 = test_density("M6.2 demo1-realistic", W_big, mass=1.0)
    all_pass &= p3

    # M6.3: active×active distance
    all_pass &= test_dist_aa("M6.3 small", 50)
    all_pass &= test_dist_aa("M6.3 demo1-realistic", 343)

    # M6.4: density constraint gradient
    all_pass &= test_density_grad("M6.4 small (20+30)", 20, 30, h=5.0,
                                  mass=1.0, rho_rest=10.0)
    all_pass &= test_density_grad("M6.4 demo1-realistic (343+17498)",
                                  343, 17498, h=5.0,
                                  mass=1.0, rho_rest=10.0)

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
