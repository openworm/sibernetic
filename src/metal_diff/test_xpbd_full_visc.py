"""Validate xpbd_full_fwd / _bwd with pair_forces enabled.

Same setup as test_xpbd_full.py but with visc_pair_coef > 0 — exercises
the new pair_forces forward/backward dispatches inside the multi-step
trainable chain.

For pair_forces backward, we treat density as stop-gradient (the chain
doesn't backprop through density-compute via the viscosity divisor).
The dL/d(rho_rest) gradient should still be correct because rho_rest
appears in the density constraint chain (which is the dominant path
for rho gradient flow).
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def _w(p, a):
    a.astype(np.float32).tofile(p)


def fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        pos0, vel0, static_p, sim_scale, visc_pair_coef):
    p_x, p_v, p_s = "/tmp/xfv_x.bin", "/tmp/xfv_v.bin", "/tmp/xfv_s.bin"
    p_state = "/tmp/xfv_state.bin"
    _w(p_x, pos0); _w(p_v, vel0); _w(p_s, static_p)
    subprocess.run(
        [BIN, "xpbd_full_fwd",
         str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
         str(dt), str(g_y), str(alpha),
         p_x, p_v, p_s, p_state,
         str(sim_scale), str(visc_pair_coef)],
        check=True,
    )
    extra = 1 if visc_pair_coef != 0 else 0
    per_step = n_a * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(p_state, dtype=np.float32)
    state = raw[:K * per_step]
    traj = raw[K * per_step:K * per_step + (K + 1) * n_a * 3].reshape(K + 1, n_a, 3)
    return traj, p_state


def bwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        state_path, static_p, grad_x_final, sim_scale, visc_pair_coef):
    p_s, p_gxf = "/tmp/xfv_s.bin", "/tmp/xfv_gxf.bin"
    _w(p_s, static_p); _w(p_gxf, grad_x_final)
    p_gxi, p_gvi, p_gr = "/tmp/xfv_gxi.bin", "/tmp/xfv_gvi.bin", "/tmp/xfv_gr.bin"
    subprocess.run(
        [BIN, "xpbd_full_bwd",
         str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
         str(dt), str(g_y), str(alpha),
         state_path, p_s, p_gxf, p_gxi, p_gvi, p_gr,
         str(sim_scale), str(visc_pair_coef)],
        check=True,
    )
    g_x = np.fromfile(p_gxi, dtype=np.float32).reshape(n_a, 3)
    g_v = np.fromfile(p_gvi, dtype=np.float32).reshape(n_a, 3)
    g_r = float(np.fromfile(p_gr, dtype=np.float32)[0])
    return g_x, g_v, g_r


def main():
    if not os.path.exists(BIN):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(17)
    # Same toy config as test_xpbd_full but with sim_scale=1 (toy
    # convention; pair_forces still fires via visc_pair_coef > 0).
    n_a, n_s = 6, 12
    K = 5
    h = 1.5
    mass = 1.0
    rho_rest = 0.3
    dt = 0.01
    g_y = -1.0
    alpha = 0.5
    sim_scale = 1.0
    visc_pair_coef = 1e-2  # weak but non-zero — pair forces actually fire

    pos0 = np.random.randn(n_a, 3).astype(np.float32) * 0.3
    vel0 = np.zeros((n_a, 3), dtype=np.float32)
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 0.5

    print(f"=== xpbd_full WITH pair_forces (visc_pair_coef={visc_pair_coef}) ===")
    print(f"  n_active={n_a}, n_static={n_s}, K={K}, dt={dt}, ρ={rho_rest}")

    # Forward at true ρ
    traj_obs, _ = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                      pos0, vel0, static_p, sim_scale, visc_pair_coef)
    pos_obs = traj_obs[-1].copy()

    # Forward at perturbed ρ
    bad_rho = 0.5
    traj_bad, state_path = fwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                pos0, vel0, static_p, sim_scale, visc_pair_coef)
    grad_x_final = (traj_bad[-1] - pos_obs).astype(np.float32)
    L = 0.5 * (grad_x_final ** 2).sum()
    g_x, g_v, g_r_kernel = bwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                state_path, static_p, grad_x_final,
                                sim_scale, visc_pair_coef)

    eps = 1e-3
    traj_p, _ = fwd(n_a, n_s, K, h, mass, bad_rho + eps, dt, g_y, alpha,
                    pos0, vel0, static_p, sim_scale, visc_pair_coef)
    traj_m, _ = fwd(n_a, n_s, K, h, mass, bad_rho - eps, dt, g_y, alpha,
                    pos0, vel0, static_p, sim_scale, visc_pair_coef)
    L_p = 0.5 * ((traj_p[-1] - pos_obs) ** 2).sum()
    L_m = 0.5 * ((traj_m[-1] - pos_obs) ** 2).sum()
    g_r_num = (L_p - L_m) / (2 * eps)

    print(f"  At ρ={bad_rho}: L = {L:.4e}")
    print(f"  ∂L/∂ρ kernel    = {g_r_kernel:+.6e}")
    print(f"  ∂L/∂ρ numerical = {g_r_num:+.6e}")
    rel = abs(g_r_kernel - g_r_num) / max(abs(g_r_num), 1e-9)
    print(f"  relative error  = {rel:.3e}")
    grad_pass = rel < 0.1
    print(f"  {'[PASS]' if grad_pass else '[FAIL]'} multi-step ∂L/∂ρ matches finite-diff (with pair_forces)")

    return 0 if grad_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
