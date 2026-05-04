"""Random-search fallback for parameter tuning.

Useful when SGD gets stuck (e.g., the cube-fall loss landscape has
chaotic regions where finite-diff gradients are noisy). Samples
random parameter triples in log-space within bounds and keeps the
best seen so far.

Less efficient than Adam at smooth-loss problems, but more robust
when the loss has discontinuities or chaotic regions.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

# Reuse the loss function from sgd_tune.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sgd_tune import loss, TARGET_DELTA_Y, TARGET_EXTENT_Y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-trials', type=int, default=40)
    ap.add_argument('--n-sim-steps', type=int, default=50000)
    # Sample bounds in log space (decimal log).
    ap.add_argument('--rho-log10-min', type=float, default=-13.7)
    ap.add_argument('--rho-log10-max', type=float, default=-12.0)
    ap.add_argument('--K-log10-min', type=float, default=2.0)
    ap.add_argument('--K-log10-max', type=float, default=3.5)
    ap.add_argument('--v-log10-min', type=float, default=-6.0)
    ap.add_argument('--v-log10-max', type=float, default=-3.5)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    rng = np.random.RandomState(args.seed)

    print(f"# Random search ({args.n_trials} trials, K_steps={args.n_sim_steps})")
    print(f"# Targets: Δmean_y={TARGET_DELTA_Y}, extent_y={TARGET_EXTENT_Y}")
    print(f"# Bounds (log10): rho [{args.rho_log10_min},{args.rho_log10_max}] "
          f"K [{args.K_log10_min},{args.K_log10_max}] v [{args.v_log10_min},{args.v_log10_max}]")
    print()

    best_L = float('inf')
    best_params = None
    history = []
    for it in range(args.n_trials):
        rho = 10 ** rng.uniform(args.rho_log10_min, args.rho_log10_max)
        K = 10 ** rng.uniform(args.K_log10_min, args.K_log10_max)
        v = 10 ** rng.uniform(args.v_log10_min, args.v_log10_max)

        t0 = time.perf_counter()
        L, dm, ext, err_dm, err_ext = loss(rho, K, v, n_steps=args.n_sim_steps)
        elapsed = time.perf_counter() - t0
        is_best = L < best_L
        marker = " ✓" if is_best else ""
        print(f"[{it:3d}] L={L:.4e}  Δm_y err {err_dm*100:5.2f}%  "
              f"ext_y err {err_ext*100:5.2f}%  "
              f"rho={rho:.3e} K={K:.1f} v={v:.3e}  ({elapsed:.0f}s){marker}")
        history.append({
            'trial': it, 'L': float(L),
            'delta_y': float(dm), 'extent_y': float(ext),
            'err_delta_pct': float(err_dm * 100), 'err_extent_pct': float(err_ext * 100),
            'rho': float(rho), 'K': float(K), 'v': float(v),
        })
        if is_best:
            best_L = L
            best_params = (rho, K, v)
            best_metrics = (dm, ext, err_dm, err_ext)
        if err_dm < 0.01 and err_ext < 0.01:
            print(f"\n[CONVERGED] both errors under 1% at trial {it}")
            break

    print()
    print(f"Best L={best_L:.4e}")
    if best_params is not None:
        rho, K, v = best_params
        dm, ext, err_dm, err_ext = best_metrics
        print(f"  rho={rho:.6e} K={K:.3f} v={v:.6e}")
        print(f"  Δm_y={dm:+.3f} (err {err_dm*100:.2f}%)")
        print(f"  ext_y={ext:.3f} (err {err_ext*100:.2f}%)")

    out = "/tmp/random_search_history.json"
    with open(out, 'w') as f:
        json.dump(history, f, indent=2)
    print(f"\nHistory: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
