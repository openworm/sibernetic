"""TRUE analytic-gradient SGD using xpbd_full_fwd / xpbd_full_bwd directly.

This is what the differentiable substrate was designed for. No FD,
no random search. Each iteration:
  1. xpbd_full_fwd at current (rho, K, v) → trajectory + state file
  2. Read final positions from trajectory frame
  3. Compute scalar loss L = w_dm·(Δm_y - target)² + w_ext·(ext_y - target)²
  4. Compute ∂L/∂x_final per particle (analytic chain through Δm_y, ext_y)
  5. xpbd_full_bwd → analytic ∂L/∂(rho_rest, spring_K, visc_pair_coef)
  6. Adam update in log-space, repeat

Gradients are EXACT (modulo the ~17% stop-gradient through density on
visc_K — see test_param_grads.py). No FD noise; chaos doesn't hurt
the gradient itself, only how the optimizer responds to it.

Cost: K=50000 sim steps → ~90s fwd + ~8min bwd = ~10min per iter.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")

TARGET_DELTA_Y = -36.75
TARGET_EXTENT_Y_RATIO = 1.085


def load_demo1_inputs():
    """Set up demo1 binary buffers; return info dict."""
    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers
    return load_to_metal_buffers('demo1', out_dir='/tmp/demo1_metal', h=3.34)


def run_fwd(info, rho, K, v, n_steps, dt=2e-5, sim_scale=7.4e-6,
            mass=2e-12, h=3.34, alpha_dens=1e-3, gravity_y=-9.8, floor_y=0.0):
    """xpbd_full_fwd → state file. Returns final positions [n_active, 3]."""
    state_path = "/tmp/sgd_true_state.bin"
    cmd = [BIN, "xpbd_full_fwd",
           str(info['n_active']), str(info['n_static']), str(n_steps),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity_y), str(alpha_dens),
           info['paths']['pos_active'],
           info['paths']['vel_active'],
           info['paths']['pos_static'],
           state_path,
           f"{sim_scale:.15e}", f"{v:.15e}",
           f"{K:.15e}", info['paths']['bonds'],
           f"{floor_y:.6f}"]
    subprocess.run(cmd, check=True, capture_output=True)
    # State file layout: K * per_step + (K+1) * n*3 + n*3
    n_active = info['n_active']
    extra = (1 if v != 0 else 0) + 1  # +1 for floor (always set here)
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(state_path, dtype=np.float32)
    traj_off = n_steps * per_step
    traj = raw[traj_off:traj_off + (n_steps + 1) * n_active * 3].reshape(n_steps + 1, n_active, 3)
    return traj, state_path


def compute_metrics(traj_final, n_elastic):
    """Compute Δmean_y and extent_y_ratio for the elastic cube."""
    pos_init = np.fromfile('/tmp/demo1_metal/demo1_pos_active.bin',
                            dtype=np.float32).reshape(-1, 3)
    init_mean_y = float(pos_init[:, 1].mean())
    elastic_init = pos_init[:n_elastic]
    init_ext_y = float(elastic_init[:, 1].max() - elastic_init[:, 1].min())
    final_mean_y = float(traj_final[:, 1].mean())
    elastic_fin = traj_final[:n_elastic]
    final_ext_y = float(elastic_fin[:, 1].max() - elastic_fin[:, 1].min())
    return final_mean_y - init_mean_y, final_ext_y / init_ext_y


def loss_and_grad_x(traj_final, n_elastic, w_dm=1.0, w_ext=1.0):
    """Compute L (percentage-error) and ∂L/∂x_final per particle.

    L = w_dm · err_dm² + w_ext · err_ext²
    where err_dm = (Δm_y - target_dm) / |target_dm|
          err_ext = (ext_y - target_ext) / |target_ext|

    ∂L/∂x_final[i, 1] (only y-component nonzero):
    - From Δm_y: 2·w_dm·err_dm/|target_dm| · ∂Δm_y/∂x[i, 1]
                where ∂Δm_y/∂x[i, 1] = 1/n_active (mean over all active)
    - From ext_y: 2·w_ext·err_ext/|target_ext| · ∂ext_y/∂x[i, 1] / init_ext
                where ∂ext_y_raw/∂x[i, 1] = +1 if i = argmax_y(elastic),
                                            -1 if i = argmin_y(elastic),
                                             0 else (only elastic, only y-axis)

    Returns (L, dL/dx [n_active, 3], dm, ext_ratio).
    """
    pos_init = np.fromfile('/tmp/demo1_metal/demo1_pos_active.bin',
                            dtype=np.float32).reshape(-1, 3)
    init_mean_y = float(pos_init[:, 1].mean())
    elastic_init = pos_init[:n_elastic]
    init_ext_y = float(elastic_init[:, 1].max() - elastic_init[:, 1].min())

    n_active = traj_final.shape[0]
    final_mean_y = float(traj_final[:, 1].mean())
    elastic_fin = traj_final[:n_elastic]
    argmax = int(elastic_fin[:, 1].argmax())  # within elastic
    argmin = int(elastic_fin[:, 1].argmin())
    final_ext_y_raw = float(elastic_fin[argmax, 1] - elastic_fin[argmin, 1])
    final_ext_y_ratio = final_ext_y_raw / init_ext_y

    dm = final_mean_y - init_mean_y
    err_dm = (dm - TARGET_DELTA_Y) / abs(TARGET_DELTA_Y)
    err_ext = (final_ext_y_ratio - TARGET_EXTENT_Y_RATIO) / abs(TARGET_EXTENT_Y_RATIO)
    L = w_dm * err_dm * err_dm + w_ext * err_ext * err_ext

    # Gradient
    grad = np.zeros_like(traj_final)
    # ∂L/∂(Δm_y) = 2·w_dm·err_dm/|target_dm|
    dL_ddm = 2.0 * w_dm * err_dm / abs(TARGET_DELTA_Y)
    # Δm_y = mean(pos[:,1]) - init_mean_y; ∂Δm_y/∂x[i,1] = 1/n_active
    grad[:, 1] += dL_ddm / n_active
    # ∂L/∂(ext_y_ratio) = 2·w_ext·err_ext/|target_ext|
    dL_dext = 2.0 * w_ext * err_ext / abs(TARGET_EXTENT_Y_RATIO)
    # ext_y_ratio = (max - min) / init_ext; ∂/∂x[argmax,1] = +1/init_ext, ∂/∂x[argmin,1] = -1/init_ext
    grad[argmax, 1] += dL_dext / init_ext_y
    grad[argmin, 1] -= dL_dext / init_ext_y

    return L, grad, dm, final_ext_y_ratio


def run_bwd(info, rho, K, v, n_steps, grad_x_final, dt=2e-5, sim_scale=7.4e-6,
            mass=2e-12, h=3.34, alpha_dens=1e-3, gravity_y=-9.8, floor_y=0.0):
    """xpbd_full_bwd → (∂L/∂rho, ∂L/∂spring_K, ∂L/∂visc_K)."""
    np.asarray(grad_x_final, dtype=np.float32).tofile('/tmp/sgd_true_gxf.bin')
    cmd = [BIN, "xpbd_full_bwd",
           str(info['n_active']), str(info['n_static']), str(n_steps),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity_y), str(alpha_dens),
           '/tmp/sgd_true_state.bin',
           info['paths']['pos_static'],
           '/tmp/sgd_true_gxf.bin',
           '/tmp/sgd_true_gx0.bin',
           '/tmp/sgd_true_gv0.bin',
           '/tmp/sgd_true_grho.bin',
           f"{sim_scale:.15e}", f"{v:.15e}",
           f"{K:.15e}", info['paths']['bonds'],
           '/tmp/sgd_true_gK.bin',
           '/tmp/sgd_true_gvK.bin',
           f"{floor_y:.6f}"]
    subprocess.run(cmd, check=True, capture_output=True)
    g_rho = float(np.fromfile('/tmp/sgd_true_grho.bin', dtype=np.float32)[0])
    g_K = float(np.fromfile('/tmp/sgd_true_gK.bin', dtype=np.float32)[0])
    g_v = float(np.fromfile('/tmp/sgd_true_gvK.bin', dtype=np.float32)[0])
    return g_rho, g_K, g_v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-steps', type=int, default=15)
    ap.add_argument('--lr', type=float, default=0.05)
    ap.add_argument('--n-sim-steps', type=int, default=50000)
    ap.add_argument('--rho-init', type=float, default=2.5e-13)
    ap.add_argument('--K-init', type=float, default=1000.0)
    ap.add_argument('--v-init', type=float, default=5e-5)
    args = ap.parse_args()

    info = load_demo1_inputs()
    n_elastic = info['n_elastic']

    rho, K, v = args.rho_init, args.K_init, args.v_init
    print(f"# True analytic-gradient SGD (xpbd_full_fwd + xpbd_full_bwd)")
    print(f"# n_active={info['n_active']} n_static={info['n_static']}")
    print(f"# n_elastic={n_elastic} (cube particles for extent metric)")
    print(f"# Targets: Δm_y={TARGET_DELTA_Y}  ext_y_ratio={TARGET_EXTENT_Y_RATIO}")
    print(f"# Initial: rho={rho:.4e} K={K:.2f} v={v:.4e}")
    print(f"# K_steps={args.n_sim_steps} lr={args.lr}")
    print()

    m = np.zeros(3); var = np.zeros(3)
    beta1, beta2, eps = 0.9, 0.999, 1e-8
    history = []
    best_L = float('inf'); best_params = (rho, K, v)

    for it in range(args.max_steps):
        t0 = time.perf_counter()
        traj, _ = run_fwd(info, rho, K, v, args.n_sim_steps)
        traj_final = traj[-1]
        L, grad_x, dm, ext_ratio = loss_and_grad_x(traj_final, n_elastic)
        err_dm = abs(dm - TARGET_DELTA_Y) / abs(TARGET_DELTA_Y)
        err_ext = abs(ext_ratio - TARGET_EXTENT_Y_RATIO) / abs(TARGET_EXTENT_Y_RATIO)
        elapsed_fwd = time.perf_counter() - t0

        if L < best_L: best_L = L; best_params = (rho, K, v); marker = " ★"
        else: marker = ""

        print(f"[step {it}] L={L:.4e}  Δm err {err_dm*100:5.2f}%  "
              f"ext err {err_ext*100:5.2f}%  "
              f"rho={rho:.4e} K={K:.2f} v={v:.4e}  fwd={elapsed_fwd:.1f}s{marker}")
        history.append({'step': it, 'L': float(L),
                        'delta_y': float(dm), 'extent_y_ratio': float(ext_ratio),
                        'err_delta_pct': float(err_dm*100),
                        'err_extent_pct': float(err_ext*100),
                        'rho': float(rho), 'K': float(K), 'v': float(v)})

        if err_dm < 0.01 and err_ext < 0.01:
            print(f"\n[CONVERGED] both errors < 1% at step {it}")
            break

        t0 = time.perf_counter()
        g_rho, g_K, g_v = run_bwd(info, rho, K, v, args.n_sim_steps, grad_x)
        elapsed_bwd = time.perf_counter() - t0

        # Convert to log-space gradients via chain rule:
        # ∂L/∂(log p) = p · ∂L/∂p
        g = np.array([g_rho * rho, g_K * K, g_v * v])
        print(f"   grad: rho_log={g[0]:+.3e} K_log={g[1]:+.3e} v_log={g[2]:+.3e}  bwd={elapsed_bwd:.1f}s")

        m = beta1 * m + (1 - beta1) * g
        var = beta2 * var + (1 - beta2) * (g * g)
        m_hat = m / (1 - beta1 ** (it + 1))
        v_hat = var / (1 - beta2 ** (it + 1))
        delta = -args.lr * m_hat / (np.sqrt(v_hat) + eps)
        rho = float(np.exp(np.log(rho) + delta[0]))
        K   = float(np.exp(np.log(K)   + delta[1]))
        v   = float(np.exp(np.log(v)   + delta[2]))

    print(f"\n=== Best: L={best_L:.4e} ===")
    rho_b, K_b, v_b = best_params
    print(f"  rho={rho_b:.6e} K={K_b:.3f} v={v_b:.6e}")

    out = "/tmp/sgd_true_history.json"
    with open(out, 'w') as f: json.dump(history, f, indent=2)
    print(f"History: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
