"""Validate analytic ∂L/∂(spring_K) and ∂L/∂(visc_pair_coef) via finite-diff.

Same toy as test_xpbd_full_spring but checks the new param-gradient
outputs from xpbd_full_bwd against numerical FD.
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
    _w("/tmp/tpg_x.bin", pos0); _w("/tmp/tpg_v.bin", vel0); _w("/tmp/tpg_s.bin", static_p)
    cmd = [BIN, "xpbd_full_fwd",
           str(n_a), str(n_s), str(K), str(h), str(mass), f"{rho_rest:.15e}",
           str(dt), str(g_y), str(alpha),
           "/tmp/tpg_x.bin", "/tmp/tpg_v.bin", "/tmp/tpg_s.bin", "/tmp/tpg_state.bin",
           f"{sim_scale:.15e}", f"{visc_pair_coef:.15e}",
           f"{spring_K:.15e}", bonds_path]
    subprocess.run(cmd, check=True, capture_output=True)
    extra = 1 if visc_pair_coef != 0 else 0
    per_step = n_a * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile("/tmp/tpg_state.bin", dtype=np.float32)
    state = raw[:K * per_step]
    traj = raw[K * per_step:K * per_step + (K + 1) * n_a * 3].reshape(K + 1, n_a, 3)
    return traj


def bwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
        static_p, grad_x_final, sim_scale, visc_pair_coef,
        spring_K, bonds_path):
    _w("/tmp/tpg_s.bin", static_p); _w("/tmp/tpg_gxf.bin", grad_x_final)
    cmd = [BIN, "xpbd_full_bwd",
           str(n_a), str(n_s), str(K), str(h), str(mass), f"{rho_rest:.15e}",
           str(dt), str(g_y), str(alpha),
           "/tmp/tpg_state.bin", "/tmp/tpg_s.bin", "/tmp/tpg_gxf.bin",
           "/tmp/tpg_gxi.bin", "/tmp/tpg_gvi.bin", "/tmp/tpg_gr.bin",
           f"{sim_scale:.15e}", f"{visc_pair_coef:.15e}",
           f"{spring_K:.15e}", bonds_path,
           "/tmp/tpg_gK.bin", "/tmp/tpg_gv.bin"]
    subprocess.run(cmd, check=True, capture_output=True)
    g_K = float(np.fromfile("/tmp/tpg_gK.bin", dtype=np.float32)[0])
    g_v = float(np.fromfile("/tmp/tpg_gv.bin", dtype=np.float32)[0])
    return g_K, g_v


def main():
    if not os.path.exists(BIN):
        print(f"sib_metal binary missing")
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

    bonds = np.zeros(3, dtype=[('i', 'i4'), ('j', 'i4'),
                                ('rest', 'f4'), ('pad', 'f4')])
    for k, (i, j) in enumerate([(0, 1), (1, 2), (2, 3)]):
        bonds[k] = (i, j, 0.5, 0.0)
    bonds_path = "/tmp/tpg_bonds.bin"
    bonds.tofile(bonds_path)

    print(f"=== Param-grad validation (springs + visc) ===")
    print(f"  K={spring_K} v_pair={visc_pair_coef} ρ={rho_rest}")

    # Reference trajectory at "true" params.
    traj_obs = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                   pos0, vel0, static_p, sim_scale, visc_pair_coef,
                   spring_K, bonds_path)
    pos_obs = traj_obs[-1].copy()

    # Perturbed params, compute gradient.
    bad_K = 60.0
    bad_v = 1.5e-2
    traj_bad = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                   pos0, vel0, static_p, sim_scale, bad_v,
                   bad_K, bonds_path)
    grad_x_final = (traj_bad[-1] - pos_obs).astype(np.float32)
    L = 0.5 * (grad_x_final ** 2).sum()
    g_K_kernel, g_v_kernel = bwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                                  static_p, grad_x_final, sim_scale, bad_v,
                                  bad_K, bonds_path)
    print(f"  L = {L:.4e}")
    print(f"  ∂L/∂(spring_K) kernel    = {g_K_kernel:+.6e}")
    print(f"  ∂L/∂(visc_K)   kernel    = {g_v_kernel:+.6e}")

    # FD on spring_K
    eps_K = 0.5
    traj_p = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                 pos0, vel0, static_p, sim_scale, bad_v,
                 bad_K + eps_K, bonds_path)
    traj_m = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                 pos0, vel0, static_p, sim_scale, bad_v,
                 bad_K - eps_K, bonds_path)
    L_p = 0.5 * ((traj_p[-1] - pos_obs) ** 2).sum()
    L_m = 0.5 * ((traj_m[-1] - pos_obs) ** 2).sum()
    g_K_num = (L_p - L_m) / (2 * eps_K)

    # FD on visc_K — use larger eps to escape fp32 FD precision floor.
    # Loss is ~2e-4 and ∂L/∂visc_K ~ 5e-7, so L_p - L_m at small eps is
    # in the 1e-10 range, dominated by fp32 cancellation noise.
    eps_v = 5e-3
    traj_p = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                 pos0, vel0, static_p, sim_scale, bad_v + eps_v,
                 bad_K, bonds_path)
    traj_m = fwd(n_a, n_s, K, h, mass, rho_rest, dt, g_y, alpha,
                 pos0, vel0, static_p, sim_scale, bad_v - eps_v,
                 bad_K, bonds_path)
    L_p = 0.5 * ((traj_p[-1] - pos_obs) ** 2).sum()
    L_m = 0.5 * ((traj_m[-1] - pos_obs) ** 2).sum()
    g_v_num = (L_p - L_m) / (2 * eps_v)

    print(f"  ∂L/∂(spring_K) numerical = {g_K_num:+.6e}")
    print(f"  ∂L/∂(visc_K)   numerical = {g_v_num:+.6e}")

    rel_K = abs(g_K_kernel - g_K_num) / max(abs(g_K_num), 1e-12)
    rel_v = abs(g_v_kernel - g_v_num) / max(abs(g_v_num), 1e-12)
    print(f"  rel err spring_K: {rel_K:.3e}")
    print(f"  rel err visc_K:   {rel_v:.3e}")

    pass_K = rel_K < 0.1
    pass_v = rel_v < 0.1
    print(f"  {'[PASS]' if pass_K else '[FAIL]'} spring_K analytic grad")
    print(f"  {'[PASS]' if pass_v else '[FAIL]'} visc_K analytic grad")

    return 0 if (pass_K and pass_v) else 1


if __name__ == "__main__":
    sys.exit(main())
