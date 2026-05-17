"""FD-validate the M9 density-chain backward CUDA kernels.

Tests the 5 M6.x backwards:
  M6.0_bwd  dist_active_static_bwd       (∂r²/∂active)
  M6.1_bwd  wpoly6_inplace_bwd           (∂W/∂r², in-place)
  M6.2_bwd  rowsum_density_bwd           (∂density/∂W)
  M6.3_bwd  dist_active_active_bwd       (∂r²/∂active, both ends differentiable)
  M6.4_bwd  density_constraint_grad_bwd  (∂(grad_C, denom_helper)/∂active)

Strategy: for each kernel, draw a small random input, pick random upstream
gradients ω (and ψ for M6.4), evaluate L(x) = ω·output(x) numerically,
compute FD gradient via central differences, and compare to the CUDA
backward kernel's output.

Sizes kept small enough for an O(n) FD sweep over input components to run
in well under a second.

Local-only; no Metal needed.
"""
import math
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


def _run(args):
    subprocess.run([BINARY] + [str(a) for a in args], check=True)


# ─────────────────────────────  M6.0_bwd  ─────────────────────────────
def numpy_dist_as(active, static_p):
    """r²[i,j] = ||a_i - s_j||²"""
    diff = active[:, None, :] - static_p[None, :, :]
    return np.sum(diff * diff, axis=-1)


def cuda_dist_as_bwd(active, static_p, grad_r2):
    n_active, n_static = len(active), len(static_p)
    p_a = os.path.join(TMP, "m6bwd_a.bin")
    p_s = os.path.join(TMP, "m6bwd_s.bin")
    p_g = os.path.join(TMP, "m6bwd_gr2.bin")
    p_o = os.path.join(TMP, "m6bwd_gact.bin")
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    grad_r2.astype(np.float32).tofile(p_g)
    np.zeros((n_active, 3), np.float32).tofile(p_o)
    _run(["dist_active_static_bwd", n_active, n_static, p_a, p_s, p_g, p_o])
    return np.fromfile(p_o, dtype=np.float32).reshape(n_active, 3)


def test_m6_0_bwd(eps=1e-3, tol=1e-3):
    np.random.seed(1)
    n_active, n_static = 6, 9
    active   = np.random.randn(n_active, 3).astype(np.float64) * 0.8
    static_p = np.random.randn(n_static, 3).astype(np.float64) * 0.8
    grad_r2  = np.random.randn(n_active, n_static).astype(np.float64) * 0.5

    cuda_grad = cuda_dist_as_bwd(active.astype(np.float32),
                                 static_p.astype(np.float32),
                                 grad_r2.astype(np.float32))

    # FD: L(active) = sum(grad_r2 * r²(active, static))
    def loss(a):
        return float((grad_r2 * numpy_dist_as(a, static_p)).sum())
    fd = np.zeros_like(active)
    for i in range(n_active):
        for k in range(3):
            a_p = active.copy(); a_p[i, k] += eps
            a_m = active.copy(); a_m[i, k] -= eps
            fd[i, k] = (loss(a_p) - loss(a_m)) / (2 * eps)

    err_max = float(np.abs(cuda_grad.astype(np.float64) - fd).max())
    mag = float(np.abs(fd).max()) + 1e-12
    print(f"=== M6.0_bwd dist_active_static_bwd ===")
    print(f"  max |analytic|:  {np.abs(cuda_grad).max():.4e}")
    print(f"  max abs err:     {err_max:.4e}")
    print(f"  max rel err:     {err_max / mag:.4e}")
    ok = err_max / mag < tol
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.3_bwd  ─────────────────────────────
def numpy_dist_aa(active):
    diff = active[:, None, :] - active[None, :, :]
    return np.sum(diff * diff, axis=-1)


def cuda_dist_aa_bwd(active, grad_r2):
    n_active = len(active)
    p_a = os.path.join(TMP, "m6bwd_aa_a.bin")
    p_g = os.path.join(TMP, "m6bwd_aa_gr2.bin")
    p_o = os.path.join(TMP, "m6bwd_aa_gact.bin")
    active.astype(np.float32).tofile(p_a)
    grad_r2.astype(np.float32).tofile(p_g)
    np.zeros((n_active, 3), np.float32).tofile(p_o)
    _run(["dist_active_active_bwd", n_active, p_a, p_g, p_o])
    return np.fromfile(p_o, dtype=np.float32).reshape(n_active, 3)


def test_m6_3_bwd(eps=1e-3, tol=1e-3):
    np.random.seed(2)
    n_active = 7
    active  = np.random.randn(n_active, 3).astype(np.float64) * 0.7
    # Asymmetric upstream so we exercise both (i,j) and (j,i) contributions.
    grad_r2 = np.random.randn(n_active, n_active).astype(np.float64) * 0.5
    np.fill_diagonal(grad_r2, 0.0)  # diagonal r² is identically 0 anyway

    cuda_grad = cuda_dist_aa_bwd(active.astype(np.float32),
                                 grad_r2.astype(np.float32))

    def loss(a):
        return float((grad_r2 * numpy_dist_aa(a)).sum())
    fd = np.zeros_like(active)
    for i in range(n_active):
        for k in range(3):
            a_p = active.copy(); a_p[i, k] += eps
            a_m = active.copy(); a_m[i, k] -= eps
            fd[i, k] = (loss(a_p) - loss(a_m)) / (2 * eps)

    err_max = float(np.abs(cuda_grad.astype(np.float64) - fd).max())
    mag = float(np.abs(fd).max()) + 1e-12
    print(f"=== M6.3_bwd dist_active_active_bwd ===")
    print(f"  max |analytic|:  {np.abs(cuda_grad).max():.4e}")
    print(f"  max abs err:     {err_max:.4e}")
    print(f"  max rel err:     {err_max / mag:.4e}")
    ok = err_max / mag < tol
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.1_bwd  ─────────────────────────────
def numpy_wpoly6(r2, h):
    h2 = h * h
    poly6_const = 315.0 / (64.0 * math.pi * h ** 9)
    out = np.where(r2 < h2, poly6_const * (h2 - r2) ** 3, 0.0)
    return out


def cuda_wpoly6_bwd_inplace(r2_saved, grad_W, h):
    """In-place: input grad_W → output grad_r2 in same buffer."""
    n = r2_saved.size
    p_r2 = os.path.join(TMP, "m6bwd_w_r2.bin")
    p_b  = os.path.join(TMP, "m6bwd_w_buf.bin")
    r2_saved.astype(np.float32).tofile(p_r2)
    grad_W.astype(np.float32).tofile(p_b)
    _run(["wpoly6_inplace_bwd", n, h, p_r2, p_b])
    return np.fromfile(p_b, dtype=np.float32).reshape(r2_saved.shape)


def test_m6_1_bwd(eps=1e-3, tol=1e-3):
    np.random.seed(3)
    n_rows, n_cols = 5, 7
    h = 1.5
    # Mix of in-range and out-of-range r² values so we hit both branches.
    r2 = np.random.uniform(0.0, h * h * 1.5, size=(n_rows, n_cols)).astype(np.float64)
    grad_W = np.random.randn(n_rows, n_cols).astype(np.float64) * 0.5

    cuda_grad = cuda_wpoly6_bwd_inplace(r2.astype(np.float32),
                                        grad_W.astype(np.float32), h)

    def loss(r2v):
        return float((grad_W * numpy_wpoly6(r2v, h)).sum())
    fd = np.zeros_like(r2)
    for i in range(n_rows):
        for j in range(n_cols):
            r_p = r2.copy(); r_p[i, j] += eps
            r_m = r2.copy(); r_m[i, j] -= eps
            fd[i, j] = (loss(r_p) - loss(r_m)) / (2 * eps)

    # Cells where r² is near h² flip the mask under perturbation — skip those
    # in the comparison since the function is not smooth there. The CUDA
    # kernel returns 0 for out-of-range cells (correct one-sided gradient);
    # FD will return ~half the analytic value for cells straddling h².
    h2 = h * h
    safe = (r2 < h2 - 2 * eps) | (r2 > h2 + 2 * eps)
    err_max = float(np.abs(cuda_grad.astype(np.float64) - fd)[safe].max())
    mag = float(np.abs(fd[safe]).max()) + 1e-12
    print(f"=== M6.1_bwd wpoly6_inplace_bwd ===")
    print(f"  smooth cells:    {int(safe.sum())} / {safe.size}")
    print(f"  max |analytic|:  {np.abs(cuda_grad).max():.4e}")
    print(f"  max abs err:     {err_max:.4e}")
    print(f"  max rel err:     {err_max / mag:.4e}")
    ok = err_max / mag < tol
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.2_bwd  ─────────────────────────────
def cuda_rowsum_bwd(grad_density, n_cols, mass):
    n_rows = grad_density.size
    p_gd = os.path.join(TMP, "m6bwd_rs_gd.bin")
    p_gw = os.path.join(TMP, "m6bwd_rs_gw.bin")
    grad_density.astype(np.float32).tofile(p_gd)
    _run(["rowsum_density_bwd", n_rows, n_cols, mass, p_gd, p_gw])
    return np.fromfile(p_gw, dtype=np.float32).reshape(n_rows, n_cols)


def test_m6_2_bwd(eps=1e-3, tol=1e-4):
    np.random.seed(4)
    n_rows, n_cols = 6, 11
    mass = 1.3
    W = np.random.randn(n_rows, n_cols).astype(np.float64) * 0.4
    grad_density = np.random.randn(n_rows).astype(np.float64) * 0.5

    cuda_grad = cuda_rowsum_bwd(grad_density.astype(np.float32), n_cols, mass)

    # density[i] = mass · sum_j W[i,j]  →  ∂L/∂W[i,j] = mass · grad_density[i].
    expected = (mass * grad_density[:, None] * np.ones((n_rows, n_cols))).astype(np.float64)

    def loss(W_arr):
        return float((grad_density * (mass * W_arr.sum(axis=1))).sum())
    fd = np.zeros_like(W)
    for i in range(n_rows):
        for j in range(n_cols):
            Wp = W.copy(); Wp[i, j] += eps
            Wm = W.copy(); Wm[i, j] -= eps
            fd[i, j] = (loss(Wp) - loss(Wm)) / (2 * eps)

    err_cuda = float(np.abs(cuda_grad.astype(np.float64) - fd).max())
    err_expected = float(np.abs(expected - fd).max())
    print(f"=== M6.2_bwd rowsum_density_bwd ===")
    print(f"  CUDA vs FD max:     {err_cuda:.4e}")
    print(f"  analytic vs FD max: {err_expected:.4e}")
    ok = (err_cuda / (np.abs(fd).max() + 1e-12) < tol)
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


# ─────────────────────────────  M6.4_bwd  ─────────────────────────────
def numpy_density_grad_full(active, static_p, h, mass, rho_rest):
    """Returns (grad_C, denom_helper). Mirrors CUDA forward density_constraint_grad."""
    spiky_const = -45.0 / (math.pi * h ** 6)
    n_active = len(active)
    grad_C = np.zeros((n_active, 3))
    denom_h = np.zeros(n_active)
    h2 = h * h
    for i in range(n_active):
        v_aa = active[i] - active
        r2_aa = (v_aa * v_aa).sum(axis=1)
        for j in range(n_active):
            if j == i: continue
            if r2_aa[j] >= h2: continue
            r = math.sqrt(r2_aa[j])
            coef = spiky_const * (h - r) ** 2 / (r + 1e-7)
            gW = coef * v_aa[j]
            grad_C[i] += mass * gW
            denom_h[i] += np.dot(gW, gW)
        v_as = active[i] - static_p
        r2_as = (v_as * v_as).sum(axis=1)
        for k in range(len(static_p)):
            if r2_as[k] >= h2: continue
            r = math.sqrt(r2_as[k])
            coef = spiky_const * (h - r) ** 2 / (r + 1e-7)
            gW = coef * v_as[k]
            grad_C[i] += mass * gW
            denom_h[i] += np.dot(gW, gW)
    return grad_C / rho_rest, denom_h


def cuda_density_grad_bwd(active, static_p, h, mass, rho_rest,
                          grad_grad_C, grad_denom_h):
    n_active, n_static = len(active), len(static_p)
    diff_aa = active[:, None, :] - active[None, :, :]
    r2_aa = (diff_aa * diff_aa).sum(axis=-1).astype(np.float32)
    diff_as = active[:, None, :] - static_p[None, :, :]
    r2_as = (diff_as * diff_as).sum(axis=-1).astype(np.float32)

    p_a   = os.path.join(TMP, "m6bwd_dg_a.bin")
    p_s   = os.path.join(TMP, "m6bwd_dg_s.bin")
    p_raa = os.path.join(TMP, "m6bwd_dg_raa.bin")
    p_ras = os.path.join(TMP, "m6bwd_dg_ras.bin")
    p_gg  = os.path.join(TMP, "m6bwd_dg_gg.bin")
    p_gp  = os.path.join(TMP, "m6bwd_dg_gp.bin")
    p_o   = os.path.join(TMP, "m6bwd_dg_gact.bin")
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    r2_aa.tofile(p_raa)
    r2_as.tofile(p_ras)
    grad_grad_C.astype(np.float32).tofile(p_gg)
    grad_denom_h.astype(np.float32).tofile(p_gp)
    np.zeros((n_active, 3), np.float32).tofile(p_o)
    _run(["density_constraint_grad_bwd",
          n_active, n_static, h, mass, rho_rest,
          p_a, p_s, p_raa, p_ras, p_gg, p_gp, p_o])
    return np.fromfile(p_o, dtype=np.float32).reshape(n_active, 3)


def test_m6_4_bwd(eps=1e-3, tol=2e-3):
    np.random.seed(5)
    n_active, n_static = 5, 8
    h = 2.0
    mass, rho_rest = 1.0, 10.0
    # Particles spread enough that pair distances are O(h) — well above the
    # R_MIN cutoff (1e-4) so the backward's r-safety guard isn't engaged.
    active   = np.random.randn(n_active, 3).astype(np.float64) * (h * 0.5)
    static_p = np.random.randn(n_static, 3).astype(np.float64) * (h * 0.5)
    omega = np.random.randn(n_active, 3).astype(np.float64) * 0.3
    psi   = np.random.randn(n_active).astype(np.float64) * 0.2

    cuda_grad = cuda_density_grad_bwd(
        active.astype(np.float32), static_p.astype(np.float32),
        h, mass, rho_rest,
        omega.astype(np.float32), psi.astype(np.float32))

    def loss(a):
        gC, dh = numpy_density_grad_full(a, static_p, h, mass, rho_rest)
        return float((omega * gC).sum() + (psi * dh).sum())
    fd = np.zeros_like(active)
    for i in range(n_active):
        for k in range(3):
            a_p = active.copy(); a_p[i, k] += eps
            a_m = active.copy(); a_m[i, k] -= eps
            fd[i, k] = (loss(a_p) - loss(a_m)) / (2 * eps)

    err_max = float(np.abs(cuda_grad.astype(np.float64) - fd).max())
    mag = float(np.abs(fd).max()) + 1e-12
    print(f"=== M6.4_bwd density_constraint_grad_bwd ===")
    print(f"  max |analytic|:  {np.abs(cuda_grad).max():.4e}")
    print(f"  max |FD|:        {np.abs(fd).max():.4e}")
    print(f"  max abs err:     {err_max:.4e}")
    print(f"  max rel err:     {err_max / mag:.4e}")
    ok = err_max / mag < tol
    print(f"  {'[PASS]' if ok else '[FAIL]'}\n")
    return ok


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1
    all_pass = True
    all_pass &= test_m6_0_bwd()
    all_pass &= test_m6_3_bwd()
    all_pass &= test_m6_1_bwd()
    all_pass &= test_m6_2_bwd()
    all_pass &= test_m6_4_bwd()
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
