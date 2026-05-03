"""Multi-step xpbd_full_fwd / _bwd validation.

Forward: K steps of (predict + density solve + update_vel).
Backward: chain through all M9 backwards + predict_bwd + update_vel_bwd
          per step, walking K steps in reverse.

Validates ∂L/∂rho_rest via finite-difference over the full multi-step
trajectory.
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def _w(p, a):
    a.astype(np.float32).tofile(p)


def fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        pos0, vel0, static_p):
    p_x = "/tmp/xfu_x.bin"; p_v = "/tmp/xfu_v.bin"; p_s = "/tmp/xfu_s.bin"
    p_state = "/tmp/xfu_state.bin"
    _w(p_x, pos0); _w(p_v, vel0); _w(p_s, static_p)
    subprocess.run(
        [BINARY, "xpbd_full_fwd",
         str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
         str(dt), str(g_y), str(alpha),
         p_x, p_v, p_s, p_state],
        check=True,
    )
    # State file: K × per_step + (K+1) × n*3 + n*3
    per_step = n_a * (3 + 3 + 1 + 3 + 1)
    raw = np.fromfile(p_state, dtype=np.float32)
    state = raw[:K * per_step]
    traj = raw[K * per_step:K * per_step + (K + 1) * n_a * 3].reshape(K + 1, n_a, 3)
    vel_final = raw[K * per_step + (K + 1) * n_a * 3:].reshape(n_a, 3)
    return traj, vel_final, p_state


def bwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        state_path, static_p, grad_x_final):
    p_s = "/tmp/xfu_s.bin"; p_gxf = "/tmp/xfu_gxf.bin"
    _w(p_s, static_p); _w(p_gxf, grad_x_final)
    p_gxi = "/tmp/xfu_gxi.bin"; p_gvi = "/tmp/xfu_gvi.bin"
    p_gr = "/tmp/xfu_gr.bin"
    subprocess.run(
        [BINARY, "xpbd_full_bwd",
         str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
         str(dt), str(g_y), str(alpha),
         state_path, p_s, p_gxf, p_gxi, p_gvi, p_gr],
        check=True,
    )
    g_x = np.fromfile(p_gxi, dtype=np.float32).reshape(n_a, 3)
    g_v = np.fromfile(p_gvi, dtype=np.float32).reshape(n_a, 3)
    g_r = float(np.fromfile(p_gr, dtype=np.float32)[0])
    return g_x, g_v, g_r


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(17)
    n_a, n_s = 6, 12
    K = 5
    h = 1.5
    mass = 1.0
    # Pick rho_rest below typical initial density so the constraint fires.
    rho_rest = 0.3
    dt = 0.01
    g_y = -1.0
    alpha = 0.5

    # Cluster particles tightly so they overlap kernel support.
    pos0 = np.random.randn(n_a, 3).astype(np.float32) * 0.3
    vel0 = np.zeros((n_a, 3), dtype=np.float32)
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 0.5

    print(f"=== Multi-step xpbd_full validation ===")
    print(f"  n_active={n_a}, n_static={n_s}, K={K}, dt={dt}, ρ={rho_rest}")

    # Generate "observed" final positions at TRUE rho_rest
    true_rho = rho_rest
    traj_obs, _, _ = fwd(n_a, n_s, K, h, mass, true_rho, dt, g_y, alpha,
                         pos0, vel0, static_p)
    pos_obs = traj_obs[-1].copy()
    print(f"  Final pos (obs): {pos_obs[0]}")

    # Compute gradient at perturbed rho via kernel
    bad_rho = 0.5
    traj_bad, _, state_path = fwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                    pos0, vel0, static_p)
    grad_x_final = (traj_bad[-1] - pos_obs).astype(np.float32)
    L = 0.5 * (grad_x_final ** 2).sum()
    g_x, g_v, g_r_kernel = bwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                state_path, static_p, grad_x_final)

    # Numerical
    eps = 1e-3
    traj_p, _, _ = fwd(n_a, n_s, K, h, mass, bad_rho + eps, dt, g_y, alpha,
                       pos0, vel0, static_p)
    traj_m, _, _ = fwd(n_a, n_s, K, h, mass, bad_rho - eps, dt, g_y, alpha,
                       pos0, vel0, static_p)
    L_p = 0.5 * ((traj_p[-1] - pos_obs) ** 2).sum()
    L_m = 0.5 * ((traj_m[-1] - pos_obs) ** 2).sum()
    g_r_num = (L_p - L_m) / (2 * eps)

    print(f"  At ρ={bad_rho}: L = {L:.4e}")
    print(f"  ∂L/∂ρ kernel    = {g_r_kernel:+.6e}")
    print(f"  ∂L/∂ρ numerical = {g_r_num:+.6e}")
    rel = abs(g_r_kernel - g_r_num) / max(abs(g_r_num), 1e-9)
    print(f"  relative error  = {rel:.3e}")
    grad_pass = rel < 0.1
    print(f"  {'[PASS]' if grad_pass else '[FAIL]'} multi-step ∂L/∂ρ matches finite-diff")

    if not grad_pass:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
