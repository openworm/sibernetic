"""SGD on demo2 with M10 membranes — tune spring_K to minimize liquid overshoot.

Uses xpbd_full_fwd / xpbd_full_bwd directly with membrane args. The
membrane backward chain is FD-validated (test_xpbd_full_membrane.py
passes K=1/2/3/5). This script extends sgd_true.py with:
  - demo2 instead of demo1
  - membrane args plumbed through fwd + bwd
  - "minimize overshoot" loss: (mean(liq_y_end) - target_y)²
  - trainable: spring_K (K) only by default

Loss intent: when membrane is too aggressive, liquid bounces into
space (mean y >> target). When membrane is missing or too weak,
liquid crashes to floor (mean y ≈ 2 << target). Target picks an
intermediate "pooling" altitude, so SGD pushes K toward the magnitude
that retains liquid near the sheet without overshoot.
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


def load_demo2_inputs(h=3.34):
    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers
    return load_to_metal_buffers('demo2', out_dir='/tmp/demo2_metal', h=h)


def run_fwd(info, rho, K, v, alpha, n_steps, dt=2e-5, sim_scale=7.4e-6,
            mass=2e-12, h=3.34, gravity_y=-9.8, floor_y=2.0, r0=1.67):
    """xpbd_full_fwd with membranes; return final pos [n_active, 3]."""
    state_path = "/tmp/sgd_demo2_state.bin"
    cmd = [BIN, "xpbd_full_fwd",
           str(info['n_active']), str(info['n_static']), str(n_steps),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity_y), f"{alpha:.15e}",
           info['paths']['pos_active'],
           info['paths']['vel_active'],
           info['paths']['pos_static'],
           state_path,
           f"{sim_scale:.15e}", f"{v:.15e}",
           f"{K:.15e}", info['paths']['bonds'],
           f"{floor_y:.6f}",
           "0.0",   # restitution
           str(info['n_membranes']), str(info['n_elastic']),
           f"{r0:.15e}",
           info['paths']['membranes'], info['paths']['pmem_index']]
    subprocess.run(cmd, check=True, capture_output=True)
    n_a = info['n_active']
    extra = (1 if v != 0 else 0) + 1 + (3 if info['n_membranes'] > 0 else 0)
    per_step = n_a * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(state_path, dtype=np.float32)
    traj_off = n_steps * per_step
    traj = raw[traj_off:traj_off + (n_steps + 1) * n_a * 3].reshape(
        n_steps + 1, n_a, 3)
    return traj[-1]


def run_bwd(info, rho, K, v, alpha, n_steps, grad_x_final, dt=2e-5,
            sim_scale=7.4e-6, mass=2e-12, h=3.34, gravity_y=-9.8,
            floor_y=2.0, r0=1.67):
    """xpbd_full_bwd with membranes."""
    np.asarray(grad_x_final, dtype=np.float32).tofile('/tmp/sgd_demo2_gxf.bin')
    cmd = [BIN, "xpbd_full_bwd",
           str(info['n_active']), str(info['n_static']), str(n_steps),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity_y), f"{alpha:.15e}",
           '/tmp/sgd_demo2_state.bin',
           info['paths']['pos_static'],
           '/tmp/sgd_demo2_gxf.bin',
           '/tmp/sgd_demo2_gx0.bin',
           '/tmp/sgd_demo2_gv0.bin',
           '/tmp/sgd_demo2_grho.bin',
           f"{sim_scale:.15e}", f"{v:.15e}",
           f"{K:.15e}", info['paths']['bonds'],
           '/tmp/sgd_demo2_gK.bin',
           '/tmp/sgd_demo2_gvK.bin',
           f"{floor_y:.6f}",
           '/tmp/sgd_demo2_galpha.bin',
           "0.0",   # restitution
           str(info['n_membranes']), str(info['n_elastic']),
           f"{r0:.15e}",
           info['paths']['membranes'], info['paths']['pmem_index']]
    subprocess.run(cmd, check=True, capture_output=True)
    g_K = float(np.fromfile('/tmp/sgd_demo2_gK.bin', dtype=np.float32)[0])
    g_alpha = float(np.fromfile('/tmp/sgd_demo2_galpha.bin', dtype=np.float32)[0])
    return g_K, g_alpha


def loss_and_grad_x(traj_final, n_elastic, target_y, w=1.0):
    """Loss = w · (mean(liq_y_end) - target_y)²
    where liquid is the active particles after the elastic block.

    ∂L/∂x[i, 1] = 2w(mean - target) · (1/n_liquid) for liquid particles
                  (other components and elastic particles get 0)
    """
    n_active = traj_final.shape[0]
    liq = traj_final[n_elastic:]
    n_liq = len(liq)
    mean_y = float(liq[:, 1].mean())
    err = mean_y - target_y
    L = w * err * err

    grad = np.zeros_like(traj_final)
    grad[n_elastic:, 1] = 2.0 * w * err / n_liq

    return L, grad, mean_y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-steps', type=int, default=10)
    ap.add_argument('--lr', type=float, default=0.05)
    ap.add_argument('--n-sim-steps', type=int, default=1500,
                    help="enough to capture the bounce window (1500 × 2e-5 = 30ms)")
    ap.add_argument('--K-init', type=float, default=2220.0)
    ap.add_argument('--alpha-init', type=float, default=1e-8)
    ap.add_argument('--target-y', type=float, default=15.0,
                    help="target final liquid mean y (sheet altitude)")
    ap.add_argument('--tbptt', type=int, default=3,
                    help="truncated BPTT — TBPTT>3 explodes at demo1+sib-scale")
    ap.add_argument('--clip-norm', type=float, default=10.0)
    ap.add_argument('--freeze-alpha', action='store_true')
    ap.add_argument('--rho', type=float, default=1000.0)
    ap.add_argument('--v', type=float, default=5e-5)
    args = ap.parse_args()

    os.environ['BWD_TBPTT'] = str(args.tbptt)
    info = load_demo2_inputs()
    n_elastic = info['n_elastic']

    K = args.K_init
    alpha = args.alpha_init
    rho = args.rho
    v = args.v

    print(f"# SGD demo2 + M10 membranes")
    print(f"# n_active={info['n_active']} n_elastic={n_elastic} n_mem={info['n_membranes']}")
    print(f"# Initial: K={K:.2f} alpha_dens={alpha:.3e} rho={rho:.0f} v={v:.3e}")
    print(f"# Target liq_mean_y={args.target_y}, TBPTT={args.tbptt}, lr={args.lr}, n_steps={args.n_sim_steps}")
    print()

    m = np.zeros(2); var = np.zeros(2)
    beta1, beta2, eps = 0.9, 0.999, 1e-8
    history = []
    best_L = float('inf'); best_K = K; best_alpha = alpha

    for it in range(args.max_steps):
        t0 = time.perf_counter()
        try:
            traj_final = run_fwd(info, rho, K, v, alpha, args.n_sim_steps)
        except subprocess.CalledProcessError as e:
            print(f"[step {it}] FWD FAIL: K={K:.2f}, skipping")
            break
        L, grad_x, mean_y = loss_and_grad_x(traj_final, n_elastic, args.target_y)
        elapsed_fwd = time.perf_counter() - t0
        if L < best_L: best_L = L; best_K = K; best_alpha = alpha; star = " *"
        else: star = ""

        print(f"[step {it}] L={L:.4e}  liq_y_end={mean_y:.2f}  "
              f"K={K:.2f} alpha={alpha:.3e}  fwd={elapsed_fwd:.1f}s{star}")

        t0 = time.perf_counter()
        g_K, g_alpha = run_bwd(info, rho, K, v, alpha, args.n_sim_steps, grad_x)
        elapsed_bwd = time.perf_counter() - t0

        # Log-space gradients
        g = np.array([g_K * K, g_alpha * alpha])
        n_bad = int(np.sum(~np.isfinite(g)))
        if n_bad: g = np.where(np.isfinite(g), g, 0.0)
        gn = float(np.linalg.norm(g))
        if gn > args.clip_norm:
            g = g * (args.clip_norm / gn)
            clip = f" (clipped {gn:.2e}→{args.clip_norm})"
        else: clip = ""
        print(f"   grad_log: K={g[0]:+.3e} alpha={g[1]:+.3e}  bwd={elapsed_bwd:.1f}s{clip}")
        if n_bad:
            print(f"   WARN: {n_bad} non-finite grad components — zeroed")

        history.append({'step': it, 'L': float(L), 'liq_y_end': float(mean_y),
                        'K': float(K), 'alpha': float(alpha)})

        m = beta1 * m + (1 - beta1) * g
        var = beta2 * var + (1 - beta2) * (g * g)
        m_hat = m / (1 - beta1 ** (it + 1))
        v_hat = var / (1 - beta2 ** (it + 1))
        delta = -args.lr * m_hat / (np.sqrt(v_hat) + eps)
        K = float(np.exp(np.log(K) + delta[0]))
        if not args.freeze_alpha:
            alpha = float(np.exp(np.log(alpha) + delta[1]))

    print(f"\n=== Best: L={best_L:.4e} ===")
    print(f"  K={best_K:.3f} alpha={best_alpha:.6e}")

    with open('/tmp/sgd_demo2_history.json', 'w') as f:
        json.dump(history, f, indent=2)
    print("History: /tmp/sgd_demo2_history.json")


if __name__ == "__main__":
    main()
