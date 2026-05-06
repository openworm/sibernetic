"""Finite-difference SGD for Metal substrate parameters vs OpenCL targets.

The analytic-gradient SGD (sgd_true.py) returns NaN gradients on
demo1 — the backward pass blows up through the chaotic post-impact
dynamics. FD avoids this by using only forward passes.

Each iteration:
  1. fwd at (K, v)                → L_0
  2. fwd at (K·(1+ε), v)          → L_K+
  3. fwd at (K, v·(1+ε))          → L_v+
  4. ∂L/∂(log K) ≈ (L_K+ - L_0) / log(1+ε)
  5. ∂L/∂(log v) ≈ (L_v+ - L_0) / log(1+ε)
  6. Adam step in log-space

3 forwards × ~1.5s = ~5s per iteration. 20 iterations = ~100s.

Loss: w_dm·(Δm_y - target_dm)² + w_ext·(ext_ratio - target_ext)²
"""
from __future__ import annotations
import argparse, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def load_demo1():
    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers
    return load_to_metal_buffers('demo1', out_dir='/tmp/demo1_metal_fd', h=3.34)


def run_fwd(info, rho, K, v, n_steps, alpha_dens=1e-3, dt=2e-5,
            sim_scale=7.4e-6, mass=2e-12, h=3.34, gravity_y=-9.8,
            floor_y=0.0):
    """Forward pass via xpbd_full_fwd, returns final positions."""
    state_path = "/tmp/sgd_fd_state.bin"
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
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"xpbd_full_fwd failed: {r.stderr[-500:]}")
    n_active = info['n_active']
    extra = (1 if v != 0 else 0) + 1
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(state_path, dtype=np.float32)
    traj_off = n_steps * per_step
    traj = raw[traj_off:traj_off + (n_steps + 1) * n_active * 3].reshape(n_steps + 1, n_active, 3)
    return traj[-1]


def compute_loss(traj_final, n_elastic, target_dm, target_ext, w_dm=1.0, w_ext=1.0):
    """Match OpenCL Δm_y and ext_y_ratio."""
    pos_init = np.fromfile('/tmp/demo1_metal_fd/demo1_pos_active.bin',
                            dtype=np.float32).reshape(-1, 3)
    init_my = float(pos_init[:, 1].mean())
    init_ext = float(pos_init[:n_elastic, 1].max() - pos_init[:n_elastic, 1].min())
    if not np.all(np.isfinite(traj_final)):
        return float('inf'), float('nan'), float('nan')
    final_my = float(traj_final[:, 1].mean())
    el_fin = traj_final[:n_elastic]
    final_ext = float(el_fin[:, 1].max() - el_fin[:, 1].min())
    dm = final_my - init_my
    rat = final_ext / init_ext
    err_dm = (dm - target_dm) / abs(target_dm)
    err_ex = (rat - target_ext) / abs(target_ext)
    L = w_dm * err_dm * err_dm + w_ext * err_ex * err_ex
    return L, dm, rat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rho', type=float, default=1000.0)
    ap.add_argument('--K-init', type=float, default=5500.0)
    ap.add_argument('--v-init', type=float, default=5e-5)
    ap.add_argument('--target-dm', type=float, default=-37.74)
    ap.add_argument('--target-ext-ratio', type=float, default=1.121)
    ap.add_argument('--n-sim-steps', type=int, default=1250)
    ap.add_argument('--max-steps', type=int, default=20)
    ap.add_argument('--lr', type=float, default=0.1)
    ap.add_argument('--eps', type=float, default=0.03,
                    help="relative perturbation for FD")
    ap.add_argument('--clip-norm', type=float, default=2.0)
    ap.add_argument('--w-dm', type=float, default=1.0)
    ap.add_argument('--w-ext', type=float, default=1.0)
    args = ap.parse_args()

    info = load_demo1()
    n_elastic = info['n_elastic']
    rho = args.rho
    K = args.K_init
    v = args.v_init
    print(f"# FD-SGD: matching demo1 OpenCL @ {args.n_sim_steps} steps")
    print(f"#   targets: Δm_y={args.target_dm}  ext_y_ratio={args.target_ext_ratio}")
    print(f"#   init:    rho={rho} K={K} v={v}")
    print(f"#   lr={args.lr}  eps={args.eps}  max_steps={args.max_steps}")
    print()

    m = np.zeros(2); var = np.zeros(2)
    beta1, beta2, eps_adam = 0.9, 0.999, 1e-8
    best_L = float('inf'); best = (K, v)
    history = []
    eps = args.eps

    for it in range(args.max_steps):
        t0 = time.perf_counter()
        # Center forward
        try:
            xf = run_fwd(info, rho, K, v, args.n_sim_steps)
        except RuntimeError as e:
            print(f"[step {it}] center fwd failed: {e}"); break
        L0, dm, rat = compute_loss(xf, n_elastic, args.target_dm, args.target_ext_ratio,
                                     args.w_dm, args.w_ext)
        if not np.isfinite(L0):
            print(f"[step {it}] center loss not finite — backing off and retrying"); break

        # K perturb
        K_p = K * (1 + eps)
        try:
            xf_K = run_fwd(info, rho, K_p, v, args.n_sim_steps)
            L_K, _, _ = compute_loss(xf_K, n_elastic, args.target_dm, args.target_ext_ratio,
                                       args.w_dm, args.w_ext)
        except RuntimeError:
            L_K = float('inf')

        # v perturb
        v_p = v * (1 + eps)
        try:
            xf_v = run_fwd(info, rho, K, v_p, args.n_sim_steps)
            L_v, _, _ = compute_loss(xf_v, n_elastic, args.target_dm, args.target_ext_ratio,
                                       args.w_dm, args.w_ext)
        except RuntimeError:
            L_v = float('inf')

        # Log-space gradient
        log_eps = np.log(1 + eps)
        g_logK = (L_K - L0) / log_eps if np.isfinite(L_K) else 0.0
        g_logv = (L_v - L0) / log_eps if np.isfinite(L_v) else 0.0
        g = np.array([g_logK, g_logv])

        elapsed = time.perf_counter() - t0
        marker = ""
        if L0 < best_L:
            best_L = L0; best = (K, v); marker = " ★"

        err_dm = abs(dm - args.target_dm) / abs(args.target_dm)
        err_ex = abs(rat - args.target_ext_ratio) / abs(args.target_ext_ratio)
        print(f"[step {it}] L={L0:.4e}  Δm={dm:+.2f}({err_dm*100:5.1f}%)  ext={rat:.3f}({err_ex*100:5.1f}%)  "
              f"K={K:>8.1f} v={v:.3e}  g_K={g[0]:+.3e} g_v={g[1]:+.3e}  t={elapsed:.1f}s{marker}")
        history.append({'it': it, 'L': L0, 'dm': dm, 'rat': rat,
                        'K': K, 'v': v, 'g_K': float(g[0]), 'g_v': float(g[1])})

        if err_dm < 0.005 and err_ex < 0.005:
            print(f"\nCONVERGED at step {it}"); break

        # Filter NaN/Inf
        g = np.where(np.isfinite(g), g, 0.0)
        # Clip
        n = float(np.linalg.norm(g))
        if n > args.clip_norm:
            g = g * (args.clip_norm / n)
        # Adam in log-space
        m = beta1 * m + (1 - beta1) * g
        var = beta2 * var + (1 - beta2) * (g * g)
        m_hat = m / (1 - beta1 ** (it + 1))
        v_hat = var / (1 - beta2 ** (it + 1))
        delta = -args.lr * m_hat / (np.sqrt(v_hat) + eps_adam)
        K = float(np.exp(np.log(K) + delta[0]))
        v = float(np.exp(np.log(v) + delta[1]))

    print(f"\n=== Best: L={best_L:.4e}  K={best[0]:.1f}  v={best[1]:.3e} ===")
    import json
    with open('/tmp/sgd_fd_history.json', 'w') as f:
        json.dump(history, f, indent=2)


if __name__ == "__main__":
    main()
