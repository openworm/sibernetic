"""Validate xpbd_full_fwd/bwd with springs enabled.

Same toy as test_xpbd_full but with `spring_K > 0` active. Exercises
the new spring_bonds_force forward + backward dispatches inside the
multi-step trainable chain.

Spring backward chain feeds dL/dpos contributions into bGx_old_new
via the same ext_accel pathway as pair_forces. dL/drho_rest should
still match finite-diff (springs don't depend on rho_rest).
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def _w(p, a, dtype=np.float32):
    a.astype(dtype).tofile(p)


def fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        pos0, vel0, static_p, sim_scale, visc_pair_coef,
        spring_K, bonds_path):
    p_x, p_v, p_s = "/tmp/xfs_x.bin", "/tmp/xfs_v.bin", "/tmp/xfs_s.bin"
    p_state = "/tmp/xfs_state.bin"
    _w(p_x, pos0); _w(p_v, vel0); _w(p_s, static_p)
    cmd = [BIN, "xpbd_full_fwd",
           str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
           str(dt), str(g_y), str(alpha),
           p_x, p_v, p_s, p_state,
           str(sim_scale), str(visc_pair_coef),
           str(spring_K), bonds_path]
    subprocess.run(cmd, check=True, capture_output=True)
    extra = 1 if visc_pair_coef != 0 else 0
    per_step = n_a * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(p_state, dtype=np.float32)
    state = raw[:K * per_step]
    traj = raw[K * per_step:K * per_step + (K + 1) * n_a * 3].reshape(K + 1, n_a, 3)
    return traj, p_state


def bwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        state_path, static_p, grad_x_final, sim_scale, visc_pair_coef,
        spring_K, bonds_path):
    p_s, p_gxf = "/tmp/xfs_s.bin", "/tmp/xfs_gxf.bin"
    _w(p_s, static_p); _w(p_gxf, grad_x_final)
    p_gxi, p_gvi, p_gr = "/tmp/xfs_gxi.bin", "/tmp/xfs_gvi.bin", "/tmp/xfs_gr.bin"
    cmd = [BIN, "xpbd_full_bwd",
           str(n_a), str(n_s), str(K), str(h), str(mass), str(rho_rest),
           str(dt), str(g_y), str(alpha),
           state_path, p_s, p_gxf, p_gxi, p_gvi, p_gr,
           str(sim_scale), str(visc_pair_coef),
           str(spring_K), bonds_path]
    subprocess.run(cmd, check=True, capture_output=True)
    g_x = np.fromfile(p_gxi, dtype=np.float32).reshape(n_a, 3)
    g_v = np.fromfile(p_gvi, dtype=np.float32).reshape(n_a, 3)
    g_r = float(np.fromfile(p_gr, dtype=np.float32)[0])
    return g_x, g_v, g_r


def main():
    if not os.path.exists(BIN):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(17)
    n_a, n_s = 6, 12
    K = 5
    h = 1.5
    mass = 1.0
    rho_rest = 0.3
    dt = 0.01
    g_y = -1.0
    alpha = 0.5
    sim_scale = 1.0
    visc_pair_coef = 1e-2
    spring_K = 50.0

    pos0 = np.random.randn(n_a, 3).astype(np.float32) * 0.3
    vel0 = np.zeros((n_a, 3), dtype=np.float32)
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 0.5

    # Make a chain of 3 bonds: 0-1, 1-2, 2-3 with rest=0.5.
    bonds = np.zeros(3, dtype=[('i', 'i4'), ('j', 'i4'),
                                ('rest', 'f4'), ('pad', 'f4')])
    for k, (i, j) in enumerate([(0, 1), (1, 2), (2, 3)]):
        bonds[k] = (i, j, 0.5, 0.0)
    bonds_path = "/tmp/xfs_bonds.bin"
    bonds.tofile(bonds_path)

    print(f"=== xpbd_full WITH springs (K={spring_K}, visc={visc_pair_coef}) ===")
    print(f"  n_active={n_a}, n_static={n_s}, K={K}, dt={dt}")

    traj_obs, _ = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                      pos0, vel0, static_p, sim_scale, visc_pair_coef,
                      spring_K, bonds_path)
    pos_obs = traj_obs[-1].copy()

    bad_rho = 0.5
    traj_bad, state_path = fwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                pos0, vel0, static_p, sim_scale, visc_pair_coef,
                                spring_K, bonds_path)
    grad_x_final = (traj_bad[-1] - pos_obs).astype(np.float32)
    L = 0.5 * (grad_x_final ** 2).sum()
    g_x, g_v, g_r_kernel = bwd(n_a, n_s, K, h, mass, bad_rho, dt, g_y, alpha,
                                state_path, static_p, grad_x_final,
                                sim_scale, visc_pair_coef,
                                spring_K, bonds_path)

    eps = 1e-3
    traj_p, _ = fwd(n_a, n_s, K, h, mass, bad_rho + eps, dt, g_y, alpha,
                    pos0, vel0, static_p, sim_scale, visc_pair_coef,
                    spring_K, bonds_path)
    traj_m, _ = fwd(n_a, n_s, K, h, mass, bad_rho - eps, dt, g_y, alpha,
                    pos0, vel0, static_p, sim_scale, visc_pair_coef,
                    spring_K, bonds_path)
    L_p = 0.5 * ((traj_p[-1] - pos_obs) ** 2).sum()
    L_m = 0.5 * ((traj_m[-1] - pos_obs) ** 2).sum()
    g_r_num = (L_p - L_m) / (2 * eps)

    print(f"  At ρ={bad_rho}: L = {L:.4e}")
    print(f"  ∂L/∂ρ kernel    = {g_r_kernel:+.6e}")
    print(f"  ∂L/∂ρ numerical = {g_r_num:+.6e}")
    rel = abs(g_r_kernel - g_r_num) / max(abs(g_r_num), 1e-9)
    print(f"  relative error  = {rel:.3e}")
    grad_pass = rel < 0.1
    print(f"  {'[PASS]' if grad_pass else '[FAIL]'} multi-step ∂L/∂ρ matches finite-diff (with springs)")

    return 0 if grad_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
