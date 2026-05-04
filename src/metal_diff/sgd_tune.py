"""SGD over Sibernetic-equivalent physics parameters (rho_rest,
spring_K, visc_pair_coef) to match OpenCL's demo1 cube-fall trajectory.

Loss = w1*(Δmean_y_metal - Δmean_y_target)^2
     + w2*(extent_y_metal - extent_y_target)^2

Targets (from OpenCL on Cloud Run L4 — see DEVELOPMENT_LOG):
- Δmean_y target: -36.75 (cube center fell this much in 1.0s)
- extent_y target: 1.085 (cube y-extent ratio at t=1.0s)

Gradients via central finite differences in log-space (params are
all positive scalars spanning many orders of magnitude). Adam-style
update with separate learning rates per parameter scale.

Per-eval cost: ~25s (50000 steps × 0.5 ms). 3 params × 2 perturbations
+ 1 baseline = 7 evals per gradient step ≈ 3 min/step.

Usage:
    python3 src/metal_diff/sgd_tune.py [--max-steps N] [--lr LR]
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
TARGET_EXTENT_Y = 1.085


def run_metal(rho_rest, spring_K, visc, n_steps=50000, dt=2e-5):
    """Run xpbd_step and return (Δmean_y, extent_y_ratio)."""
    cmd = [
        "python3", os.path.join(HERE, "run_demo1_via_metal.py"),
        "--steps", str(n_steps),
        "--dt", str(dt),
        "--rho-rest", f"{rho_rest:.15e}",
        "--spring-k", f"{spring_K:.15e}",
        "--visc-pair-coef", f"{visc:.15e}",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True,
                       cwd=os.path.dirname(HERE).rsplit('/', 1)[0])
    # Tolerate non-zero exit (script's "bounded < 1000" check trips with
    # certain param combos but the metrics may still be parseable). Only
    # fail hard if we can't even parse the output.
    delta_y, extent_y = None, None
    has_nan = False
    for line in r.stdout.split('\n'):
        if 'Δmean_y' in line:
            try:
                delta_y = float(line.split('=')[1].strip())
            except (ValueError, IndexError):
                pass
        elif 'final=(' in line:
            f_str = line.split('final=(')[1].split(')')[0]
            try:
                extent_y = float(f_str.split(',')[1].strip())
            except (ValueError, IndexError):
                pass
        elif 'NaN' in line and 'no NaN' not in line:
            has_nan = True
        elif 'mean_y=nan' in line:
            has_nan = True
    if has_nan or delta_y is None or extent_y is None or np.isnan(delta_y) or np.isnan(extent_y):
        # Penalty: return very-bad but finite metrics so SGD can recover
        return 200.0, 0.001  # Δmean_y=+200 (cube ejected up), tiny extent
    # extent ratio = final / initial (10.02)
    return delta_y, extent_y / 10.02


def loss(rho, K, v, n_steps=50000, dt=2e-5,
         w_delta=1.0, w_extent=1.0):
    """Loss in percentage-error space — both metrics equally weighted at 1%.

    Using *percentage* errors (not raw differences) means a 1% miss on
    Δmean_y contributes the same as a 1% miss on extent_y. Avoids the
    raw-magnitude trap where Δmean_y (~37 unit difference scale) dominates
    extent_y (~0.04 ratio difference scale).
    """
    dm, ext = run_metal(rho, K, v, n_steps=n_steps, dt=dt)
    err_dm = (dm - TARGET_DELTA_Y) / abs(TARGET_DELTA_Y)
    err_ext = (ext - TARGET_EXTENT_Y) / abs(TARGET_EXTENT_Y)
    L = w_delta * err_dm * err_dm + w_extent * err_ext * err_ext
    return L, dm, ext, abs(err_dm), abs(err_ext)


def fd_gradient(rho, K, v, h_rel=0.05, **loss_kwargs):
    """Central finite-diff in log-space; gradient is ∂L/∂(log param)."""
    log_rho = np.log(rho)
    log_K = np.log(K)
    log_v = np.log(v)
    h = h_rel  # ~5% perturbation in log space

    L_p_rho, *_ = loss(np.exp(log_rho + h), K, v, **loss_kwargs)
    L_m_rho, *_ = loss(np.exp(log_rho - h), K, v, **loss_kwargs)
    g_log_rho = (L_p_rho - L_m_rho) / (2 * h)

    L_p_K, *_ = loss(rho, np.exp(log_K + h), v, **loss_kwargs)
    L_m_K, *_ = loss(rho, np.exp(log_K - h), v, **loss_kwargs)
    g_log_K = (L_p_K - L_m_K) / (2 * h)

    L_p_v, *_ = loss(rho, K, np.exp(log_v + h), **loss_kwargs)
    L_m_v, *_ = loss(rho, K, np.exp(log_v - h), **loss_kwargs)
    g_log_v = (L_p_v - L_m_v) / (2 * h)

    return np.array([g_log_rho, g_log_K, g_log_v])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-steps', type=int, default=20)
    ap.add_argument('--lr', type=float, default=0.1,
                    help="learning rate in log-param space")
    ap.add_argument('--n-sim-steps', type=int, default=50000,
                    help="K — XPBD sim steps per evaluation")
    ap.add_argument('--w-delta', type=float, default=1.0)
    ap.add_argument('--w-extent', type=float, default=1.0)
    # Initial parameters (best so far from prior tuning).
    ap.add_argument('--rho-rest-init', type=float, default=4.05e-13)
    ap.add_argument('--spring-k-init', type=float, default=1000.0)
    ap.add_argument('--visc-init', type=float, default=1e-5)
    args = ap.parse_args()

    rho = args.rho_rest_init
    K = args.spring_k_init
    v = args.visc_init

    loss_kwargs = dict(n_steps=args.n_sim_steps,
                       w_delta=args.w_delta,
                       w_extent=args.w_extent)

    print(f"# SGD tuning Sibernetic-equivalent params to match OpenCL")
    print(f"# Targets: Δmean_y={TARGET_DELTA_Y}, extent_y={TARGET_EXTENT_Y}")
    print(f"# Initial: rho={rho:.3e}, K={K:.1f}, v={v:.3e}")
    print(f"# K_steps={args.n_sim_steps}, lr={args.lr}")
    print()

    # Adam state
    m = np.zeros(3)
    var = np.zeros(3)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    history = []
    for it in range(args.max_steps):
        t0 = time.perf_counter()
        L_cur, dm, ext, err_dm, err_ext = loss(rho, K, v, **loss_kwargs)
        elapsed_eval = time.perf_counter() - t0
        print(f"[step {it}] L={L_cur:.4e}  "
              f"Δm_y={dm:+.3f} (err {err_dm*100:.2f}%)  "
              f"ext_y={ext:.3f} (err {err_ext*100:.2f}%)  "
              f"rho={rho:.3e} K={K:.1f} v={v:.3e}  ({elapsed_eval:.1f}s)")
        history.append({
            'step': it, 'L': float(L_cur),
            'delta_y': float(dm), 'extent_y': float(ext),
            'err_delta_pct': float(err_dm * 100), 'err_extent_pct': float(err_ext * 100),
            'rho': float(rho), 'K': float(K), 'v': float(v),
        })

        # Convergence check: BOTH errors under 1%
        if err_dm < 0.01 and err_ext < 0.01:
            print(f"\n[CONVERGED] Both metrics under 1% at step {it}")
            break

        t0 = time.perf_counter()
        g = fd_gradient(rho, K, v, h_rel=0.05, **loss_kwargs)
        elapsed_grad = time.perf_counter() - t0
        print(f"   grad in log-space: rho={g[0]:+.3e} K={g[1]:+.3e} v={g[2]:+.3e}  ({elapsed_grad:.1f}s)")

        # Adam update
        m = beta1 * m + (1 - beta1) * g
        var = beta2 * var + (1 - beta2) * (g * g)
        m_hat = m / (1 - beta1 ** (it + 1))
        v_hat = var / (1 - beta2 ** (it + 1))
        delta = -args.lr * m_hat / (np.sqrt(v_hat) + eps)
        rho = float(np.exp(np.log(rho) + delta[0]))
        K = float(np.exp(np.log(K) + delta[1]))
        v = float(np.exp(np.log(v) + delta[2]))

    # Save history
    out_path = "/tmp/sgd_tune_history.json"
    with open(out_path, 'w') as f:
        json.dump(history, f, indent=2)
    print(f"\nHistory saved to {out_path}")
    print(f"Final: rho={rho:.15e} K={K:.3f} v={v:.15e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
