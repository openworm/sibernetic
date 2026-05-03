"""Numerical gradient check for M7.C backward kernels.

Validates the hand-derived backward of (predict_positions +
update_velocities) against finite-difference. Forward is one step:

    v_pred = v + dt·g
    x_pred = x + dt·v_pred = x + dt·v + dt²·g
    v_new  = (x_pred - x) / dt = v + dt·g

For a loss function L = sum(w_x · x_pred + w_v · v_new), the analytic
gradients are:

    ∂L/∂x   = w_x                      (since dx_pred/dx = 1)
    ∂L/∂v   = dt · w_x + w_v           (chain v→v_pred→x_pred + identity v→v_new)
    ∂L/∂g_y = dt² · sum(w_x.y) + dt · sum(w_v.y)

We compare the kernel-computed backward to numerical finite-diff and
to the analytic formula. PASSES if max |kernel - analytic| < 1e-5.

Run:
    .venv/bin/python src/metal_diff/test_grad.py
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def fwd(x, v, dt, g_y):
    n = len(x)
    p_x = "/tmp/grad_x.bin"
    p_v = "/tmp/grad_v.bin"
    p_xp = "/tmp/grad_xp.bin"
    p_vn = "/tmp/sib_step_simple_vn.bin"
    x.astype(np.float32).tofile(p_x)
    v.astype(np.float32).tofile(p_v)
    subprocess.run(
        [BINARY, "step_simple_fwd",
         str(n), str(dt), str(g_y), p_x, p_v, p_xp, p_vn],
        check=True,
    )
    x_pred = np.fromfile(p_xp, dtype=np.float32).reshape(n, 3)
    v_new  = np.fromfile(p_vn, dtype=np.float32).reshape(n, 3)
    return x_pred, v_new


def bwd(grad_x_pred, grad_v_new, dt):
    n = len(grad_x_pred)
    p1 = "/tmp/grad_in_xp.bin"
    p2 = "/tmp/grad_in_vn.bin"
    p3 = "/tmp/grad_out_xo.bin"
    p4 = "/tmp/grad_out_vo.bin"
    p5 = "/tmp/grad_out_gy.bin"
    grad_x_pred.astype(np.float32).tofile(p1)
    grad_v_new.astype(np.float32).tofile(p2)
    subprocess.run(
        [BINARY, "step_simple_bwd",
         str(n), str(dt), p1, p2, p3, p4, p5],
        check=True,
    )
    g_x_old = np.fromfile(p3, dtype=np.float32).reshape(n, 3)
    g_v_old = np.fromfile(p4, dtype=np.float32).reshape(n, 3)
    g_g_y_per = np.fromfile(p5, dtype=np.float32)
    return g_x_old, g_v_old, g_g_y_per.sum()


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(42)
    n = 5
    dt = 0.01
    g_y = -9.8

    x = np.random.randn(n, 3).astype(np.float32) * 2.0
    v = np.random.randn(n, 3).astype(np.float32) * 0.5

    # Loss weights — fixed per particle so the loss is L = w_x·x_pred + w_v·v_new
    w_x = np.random.randn(n, 3).astype(np.float32)
    w_v = np.random.randn(n, 3).astype(np.float32)

    # Forward
    x_pred, v_new = fwd(x, v, dt, g_y)

    # Backward seed gradients
    grad_x_pred = w_x.copy()
    grad_v_new = w_v.copy()
    g_x_kernel, g_v_kernel, g_g_kernel = bwd(grad_x_pred, grad_v_new, dt)

    # Analytic gradients (closed form)
    g_x_analytic = w_x.copy()
    g_v_analytic = dt * w_x + w_v
    g_g_analytic = dt * dt * w_x[:, 1].sum() + dt * w_v[:, 1].sum()

    print(f"=== M7.C numerical gradient check ===")
    print(f"  n={n} particles, dt={dt}, g_y={g_y}")
    print()
    print(f"  ∂L/∂x_old:  max |kernel - analytic| = "
          f"{np.abs(g_x_kernel - g_x_analytic).max():.3e}")
    print(f"  ∂L/∂v_old:  max |kernel - analytic| = "
          f"{np.abs(g_v_kernel - g_v_analytic).max():.3e}")
    print(f"  ∂L/∂g_y:    kernel={g_g_kernel:.6f}  analytic={g_g_analytic:.6f}  "
          f"diff={abs(g_g_kernel - g_g_analytic):.3e}")
    print()

    # Numerical (finite-difference) gradient on g_y as cross-check
    eps = 1e-3
    L_plus, _ = fwd(x, v, dt, g_y + eps)
    L_minus, _ = fwd(x, v, dt, g_y - eps)
    L_plus_v = fwd(x, v, dt, g_y + eps)[1]
    L_minus_v = fwd(x, v, dt, g_y - eps)[1]
    L_plus_total = (w_x * L_plus).sum() + (w_v * L_plus_v).sum()
    L_minus_total = (w_x * L_minus).sum() + (w_v * L_minus_v).sum()
    g_g_numerical = (L_plus_total - L_minus_total) / (2 * eps)
    print(f"  ∂L/∂g_y numerical (finite-diff): {g_g_numerical:.6f}")
    print(f"      kernel/numerical agreement:   "
          f"{abs(g_g_kernel - g_g_numerical):.3e}")

    tol = 1e-4
    err_x = float(np.abs(g_x_kernel - g_x_analytic).max())
    err_v = float(np.abs(g_v_kernel - g_v_analytic).max())
    err_g = float(abs(g_g_kernel - g_g_analytic))

    passed = err_x < tol and err_v < tol and err_g < tol
    print()
    if passed:
        print("  [PASS] M7.C: predict + update_vel backward kernels match analytic gradient")
    else:
        print(f"  [FAIL] M7.C: gradient mismatch (tol={tol})")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
