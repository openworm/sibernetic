"""SGD on one_sprig_test (single elastic on a single anchor spring) to match
OpenCL period + amplitude spec.

OpenCL gold standard:
  period = 0.725 ms (1379 Hz)
  half-amplitude = 1.694 sim units
  x, z stationary (no horizontal motion — built into the kernel structure)
  y oscillation around 18.35 sim units, no damping over 70+ cycles

Trainable: anchor_k (the anchor spring stiffness — the only spring in
the scenario). Loss: weighted squared error in (period, amplitude) vs
the OpenCL targets, both normalized.

This uses finite-diff gradients over the dump_metal_trajectory CLI
(rather than analytic xpbd_full_bwd) because the loss is a
MEASUREMENT on the trajectory (zero-crossings → period; range → amp),
not a closed-form function of the final state.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def run_sim(anchor_k, n_steps=2500, chunk=5, tag='base'):
    """Run xpbd_step on one_sprig_test; return elastic y trajectory + frame_dt."""
    out = f"/tmp/sgd_sprig_{tag}.txt"
    cmd = [
        "python3", f"{HERE}/dump_metal_trajectory.py",
        "--scenario", "one_sprig_test",
        "--steps", str(n_steps),
        "--chunk", str(chunk),
        "--rho-rest", "1.0",
        "--visc-pair-coef", "0",
        "--spring-k", "1",   # not used (no elastic-elastic bonds)
        "--anchor-k", f"{anchor_k:.6f}",
        "--alpha-dist", "0",
        "--floor-y", "0",
        "--no-membranes",
        "--vel-clamp", "0",
        "--no-box-clamp",
        "--out", out,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sim failed (anchor_k={anchor_k}): {r.stderr[-300:]}")
    return parse_elastic_y(out, chunk)


def parse_elastic_y(path, chunk):
    """Read a one_sprig_test position_buffer and return (ey, frame_dt_ms)."""
    dt = 2e-5
    frame_dt_ms = chunk * dt * 1000
    with open(path) as f:
        for _ in range(6): f.readline()
        n_e = int(f.readline()); n_l = int(f.readline()); n_b = int(f.readline())
        f.readline(); f.readline()
        n_a = n_e + n_l
        ey = []
        # First frame: n_a active + n_b boundary
        for _ in range(n_a):
            w = f.readline().split()
            if len(w) < 4: break
            ey.append(float(w[1]))
        for _ in range(n_b): f.readline()
        # Subsequent frames: each has n_a + n_b
        cur_a = 0; cur_b = 0; in_active = True
        for line in f:
            w = line.split()
            if len(w) < 4: continue
            if in_active:
                ey.append(float(w[1]))
                cur_a += 1
                if cur_a == n_a:
                    cur_a = 0; cur_b = 0; in_active = False
            else:
                cur_b += 1
                if cur_b == n_b:
                    cur_b = 0; in_active = True
    return np.array(ey), frame_dt_ms


def measure(ey, frame_dt_ms):
    """Extract (period_ms, half_amplitude) from elastic y trajectory."""
    if len(ey) < 4 or np.any(np.isnan(ey)):
        return None, None
    half_amp = (ey.max() - ey.min()) / 2.0
    y_centered = ey - ey.mean()
    zc = np.where(np.diff(np.sign(y_centered)) > 0)[0]
    if len(zc) < 2:
        return None, half_amp
    periods_ms = np.diff(zc) * frame_dt_ms
    return float(periods_ms.mean()), float(half_amp)


def loss(period_ms, half_amp, target_period=0.725, target_amp=1.694,
        w_p=1.0, w_a=0.25):
    if period_ms is None or half_amp is None:
        return 1e6
    err_p = (period_ms - target_period) / target_period
    err_a = (half_amp - target_amp) / target_amp
    return float(w_p * err_p**2 + w_a * err_a**2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-steps', type=int, default=10)
    ap.add_argument('--lr', type=float, default=0.5,
                    help="log-space lr; one_sprig_test loss is well-conditioned, can use big lr")
    ap.add_argument('--ak-init', type=float, default=2200.0,
                    help="initial anchor_k. start far from optimum (~555) to demo SGD convergence")
    ap.add_argument('--fd-eps', type=float, default=0.05)
    ap.add_argument('--n-steps', type=int, default=2500,
                    help="50 ms sim time, enough for ~70 oscillations")
    ap.add_argument('--target-period', type=float, default=0.725)
    ap.add_argument('--target-amp', type=float, default=1.694)
    args = ap.parse_args()

    AK = args.ak_init
    print(f"# SGD on one_sprig_test anchor_k")
    print(f"# Targets: period={args.target_period} ms, half_amp={args.target_amp}")
    print(f"# Initial AK={AK:.2f}, lr={args.lr} (log-space), fd_eps={args.fd_eps}")
    print()

    history = []
    best_L = float('inf'); best_AK = AK

    for it in range(args.max_steps):
        t0 = time.perf_counter()
        try:
            ey, frame_dt = run_sim(AK, args.n_steps, tag=f'i{it}_base')
        except RuntimeError as e:
            print(f"[step {it}] sim crashed: {e}")
            break
        period, amp = measure(ey, frame_dt)
        L = loss(period, amp, args.target_period, args.target_amp)

        if L < best_L: best_L = L; best_AK = AK; star = ' *'
        else: star = ''
        if period is not None:
            err_p_pct = (period - args.target_period) / args.target_period * 100
            err_a_pct = (amp - args.target_amp) / args.target_amp * 100
            print(f"[step {it}] AK={AK:>9.2f}  period={period:.4f}ms ({err_p_pct:+.1f}%)  half_amp={amp:.3f} ({err_a_pct:+.1f}%)  L={L:.4e}{star}")
        else:
            print(f"[step {it}] AK={AK:>9.2f}  diverged  L={L:.4e}{star}")

        # FD gradient in log space
        try:
            ey_p, _ = run_sim(AK * np.exp(args.fd_eps), args.n_steps, tag=f'i{it}_p')
            ey_m, _ = run_sim(AK * np.exp(-args.fd_eps), args.n_steps, tag=f'i{it}_m')
        except RuntimeError as e:
            print(f"  FD failed: {e}")
            break
        L_p = loss(*measure(ey_p, frame_dt), args.target_period, args.target_amp)
        L_m = loss(*measure(ey_m, frame_dt), args.target_period, args.target_amp)
        gAK = (L_p - L_m) / (2 * args.fd_eps)

        elapsed = time.perf_counter() - t0
        print(f"   gAK_log={gAK:+.3e}  ({elapsed:.1f}s)")

        history.append({'step': it, 'AK': float(AK), 'period_ms': period,
                        'half_amp': amp, 'L': L, 'gAK': float(gAK)})

        # Step in log space (with simple gradient clip)
        gAK_clipped = max(-10.0, min(10.0, gAK))
        AK = float(AK * np.exp(-args.lr * gAK_clipped))
        AK = max(AK, 10.0)

        # Convergence check
        if L < 1e-4:
            print(f"\n[CONVERGED at L={L:.2e}]")
            break

    print(f"\n=== Best ===  AK={best_AK:.3f}  L={best_L:.4e}")
    with open('/tmp/sgd_sprig_history.json', 'w') as f:
        json.dump(history, f, indent=2)
    print("History → /tmp/sgd_sprig_history.json")


if __name__ == "__main__":
    main()
