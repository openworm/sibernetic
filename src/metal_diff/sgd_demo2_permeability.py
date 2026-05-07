"""SGD on demo2 membrane permeability via finite-difference gradients.

Uses dump_metal_trajectory.py (which has anchors + box-clamp + initial
vel + all the recent fixes) and computes gradients by central FD over
trainable params. Each iter: 1 baseline + 2*n_params perturbed runs.

Target: match OpenCL's per-side liquid retention at t=16ms.
  - OpenCL: mem-side y_mean=56, por-side y_mean=6
  - Loss: w_m * (y_mem - 56)^2 + w_p * (y_por - 6)^2

Trainable: spring_K, anchor_k. Both in log-space (positive params).

Each iter ~5x40s (one baseline + 4 FD perturbations) = ~3min/iter.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def run_sim(K, anchor_k, n_steps=1500, scenario='demo2', tag='base'):
    """Run xpbd_step via dump_metal_trajectory; return liquid positions at end."""
    out = f"/tmp/sgd_perm_{tag}.txt"
    cmd = ["python3", f"{HERE}/dump_metal_trajectory.py",
           "--scenario", scenario,
           "--steps", str(n_steps),
           "--chunk", str(min(n_steps, 100)),
           "--rho-rest", "8e-13",
           "--visc-pair-coef", "5e-5",
           "--spring-k", f"{K:.6f}",
           "--anchor-k", f"{anchor_k:.6f}",
           "--alpha-dist", "1e-8",
           "--floor-y", "2.0",
           "--vel-clamp", "1.0",
           "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sim failed (K={K}, ak={anchor_k}): {r.stderr[-300:]}")
    # Read just the last frame
    return parse_last_frame(out)


def parse_last_frame(path):
    """Read final-frame liquid positions only."""
    with open(path) as f:
        for _ in range(6): f.readline()
        n_e = int(f.readline()); n_l = int(f.readline()); n_b = int(f.readline())
        f.readline(); f.readline()
        n_a = n_e + n_l
        # Skip to last frame: read all but track last frame's particles only
        last_frame = None
        while True:
            xyz, ts = [], []
            for _ in range(n_a):
                ln = f.readline()
                if not ln: break
                w = ln.split()
                if len(w) < 4: continue
                xyz.append([float(w[0]),float(w[1]),float(w[2])])
                ts.append(float(w[3]))
            if not xyz: break
            for _ in range(n_b): f.readline()
            last_frame = (np.array(xyz), np.array(ts))
    return last_frame


def loss(liq_xyz, init_x_above_mem, target_mem=56.0, target_por=6.0,
         w_m=1.0, w_p=1.0):
    """Loss matching OpenCL per-side y_mean at t=16ms."""
    y = liq_xyz[:, 1]
    y_mem = y[init_x_above_mem].mean()
    y_por = y[~init_x_above_mem].mean()
    err_m = y_mem - target_mem
    err_p = y_por - target_por
    L = w_m * err_m**2 + w_p * err_p**2
    return float(L), float(y_mem), float(y_por)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-steps', type=int, default=8)
    ap.add_argument('--lr', type=float, default=0.2,
                    help="log-space learning rate (since params span orders of magnitude)")
    ap.add_argument('--n-sim-steps', type=int, default=800,
                    help="800 = 16ms (OpenCL ref window)")
    ap.add_argument('--K-init', type=float, default=2200.0)
    ap.add_argument('--ak-init', type=float, default=200.0)
    ap.add_argument('--fd-eps', type=float, default=0.05,
                    help="log-space FD step (5 percent)")
    ap.add_argument('--target-mem', type=float, default=56.0)
    ap.add_argument('--target-por', type=float, default=6.0)
    args = ap.parse_args()

    # Load init liquid x to compute the membrane-side mask
    sys.path.insert(0, HERE)
    from load_config import parse_config
    cfg = parse_config('configuration/demo2')
    pos0 = cfg['pos']
    liq_init_x = pos0[(pos0[:,3]>=1.0)&(pos0[:,3]<2.0), 0]
    above_mem = liq_init_x < 60.0  # split between membraned vs porous halves
    print(f"# n_liquid={len(liq_init_x)}, above_mem={above_mem.sum()}, above_por={(~above_mem).sum()}")
    print(f"# Targets: y_mem={args.target_mem}, y_por={args.target_por}")
    print(f"# n_sim_steps={args.n_sim_steps} (= {args.n_sim_steps*2e-5*1000:.1f}ms)")
    print()

    K = args.K_init
    AK = args.ak_init
    history = []
    best_L = float('inf'); best_K = K; best_AK = AK

    for it in range(args.max_steps):
        t0 = time.perf_counter()

        # Baseline
        last = run_sim(K, AK, args.n_sim_steps, tag=f'i{it}_base')
        liq = last[0][(last[1]>=1.0)&(last[1]<2.0)]
        if len(liq) != len(liq_init_x):
            print(f"[step {it}] particle count mismatch ({len(liq)} vs {len(liq_init_x)}), skipping")
            break
        L_base, ym_base, yp_base = loss(liq, above_mem, args.target_mem, args.target_por)

        if L_base < best_L: best_L = L_base; best_K = K; best_AK = AK; star = ' *'
        else: star = ''
        print(f"[step {it}] L={L_base:.4e}  y_mem={ym_base:.2f} (target {args.target_mem})  "
              f"y_por={yp_base:.2f} (target {args.target_por})  K={K:.2f} AK={AK:.2f}  fwd1={time.perf_counter()-t0:.1f}s{star}")

        # FD gradient in log space
        gK = gAK = 0.0
        # K +/-
        K_p = K * np.exp(args.fd_eps); K_m = K * np.exp(-args.fd_eps)
        last_Kp = run_sim(K_p, AK, args.n_sim_steps, tag=f'i{it}_Kp')
        liq_Kp = last_Kp[0][(last_Kp[1]>=1.0)&(last_Kp[1]<2.0)]
        L_Kp, _, _ = loss(liq_Kp, above_mem, args.target_mem, args.target_por)
        last_Km = run_sim(K_m, AK, args.n_sim_steps, tag=f'i{it}_Km')
        liq_Km = last_Km[0][(last_Km[1]>=1.0)&(last_Km[1]<2.0)]
        L_Km, _, _ = loss(liq_Km, above_mem, args.target_mem, args.target_por)
        gK = (L_Kp - L_Km) / (2 * args.fd_eps)

        # AK +/-
        AK_p = AK * np.exp(args.fd_eps); AK_m = AK * np.exp(-args.fd_eps)
        last_AKp = run_sim(K, AK_p, args.n_sim_steps, tag=f'i{it}_AKp')
        liq_AKp = last_AKp[0][(last_AKp[1]>=1.0)&(last_AKp[1]<2.0)]
        L_AKp, _, _ = loss(liq_AKp, above_mem, args.target_mem, args.target_por)
        last_AKm = run_sim(K, AK_m, args.n_sim_steps, tag=f'i{it}_AKm')
        liq_AKm = last_AKm[0][(last_AKm[1]>=1.0)&(last_AKm[1]<2.0)]
        L_AKm, _, _ = loss(liq_AKm, above_mem, args.target_mem, args.target_por)
        gAK = (L_AKp - L_AKm) / (2 * args.fd_eps)

        elapsed = time.perf_counter() - t0
        # Clip gradients (log-space norm)
        g_norm = float(np.sqrt(gK**2 + gAK**2))
        max_norm = 5.0
        if g_norm > max_norm:
            gK *= max_norm / g_norm
            gAK *= max_norm / g_norm
            clip_msg = f' (clipped {g_norm:.2e})'
        else:
            clip_msg = ''
        print(f"   gradlog: gK={gK:+.3e} gAK={gAK:+.3e}  total={elapsed:.1f}s{clip_msg}")

        history.append({'step': it, 'L': L_base, 'y_mem': ym_base, 'y_por': yp_base,
                        'K': float(K), 'AK': float(AK)})

        # Update in log space
        K = float(K * np.exp(-args.lr * gK))
        AK = float(AK * np.exp(-args.lr * gAK))
        # Hard floor on parameters
        K = max(K, 100.0)
        AK = max(AK, 10.0)

    print(f"\n=== Best: L={best_L:.4e} ===")
    print(f"  K={best_K:.3f} AK={best_AK:.3f}")
    with open('/tmp/sgd_perm_history.json', 'w') as f:
        json.dump(history, f, indent=2)
    print("History: /tmp/sgd_perm_history.json")


if __name__ == "__main__":
    main()
