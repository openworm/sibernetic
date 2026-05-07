"""SGD harness for worm-config Metal port (e.g. worm_alone_half_resolution).

Loss = weighted L2 between Metal and OpenCL elastic-particle trajectories,
sampled at multiple time points and aggregated.

Trainable: spring_K, anchor_k, rho_rest. Use FD-SGD over the
dump_metal_trajectory CLI (same pattern as sgd_one_sprig and
sgd_demo2_permeability).

OpenCL reference: a position_buffer.txt at /tmp/worm_opencl_position.txt
with the same number of frames as the Metal output.

Targets (computed from OpenCL):
  - Per-particle mean position at frame f for f in {f0, f1, f2, ...}
  - Loss = mean_f mean_i ||x_metal[f, i] - x_opencl[f, i]||^2

We weight by elastic-only (skip liquid, skip boundary) since the worm
shape parity is the visual goal.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def parse_traj(path):
    """Return (active_traj [n_frames, n_active, 3], n_e, n_l, frame_dt_ms)."""
    with open(path) as f:
        n_lines = sum(1 for _ in f)
    with open(path) as f:
        for _ in range(6): f.readline()
        n_e = int(f.readline()); n_l = int(f.readline()); n_b = int(f.readline())
        dt = float(f.readline()); ls = float(f.readline())
        n_a = n_e + n_l
        frame_dt_ms = dt * ls * 1000

        remaining = n_lines - 11 - n_a - n_b
        per_frame = (n_a + n_b) if (remaining > 0 and (n_a + n_b) > 1 and
                                     remaining % (n_a + n_b) == 0) else n_a
        n_extra = remaining // per_frame

        # Frame 0
        traj = []
        f0 = []
        for _ in range(n_a):
            w = f.readline().split()
            f0.append([float(w[0]), float(w[1]), float(w[2])])
        traj.append(f0)
        for _ in range(n_b): f.readline()
        # Subsequent
        for _ in range(n_extra):
            fk = []
            for _ in range(n_a):
                line = f.readline()
                if not line: break
                w = line.split()
                if len(w) < 4: break
                fk.append([float(w[0]), float(w[1]), float(w[2])])
            if len(fk) != n_a: break
            traj.append(fk)
            if per_frame == n_a + n_b:
                for _ in range(n_b): f.readline()
    return np.array(traj), n_e, n_l, frame_dt_ms


def loss_vs_opencl(metal_traj, opencl_traj, n_e, n_frames=None):
    """Mean per-elastic-particle squared distance, averaged across frames."""
    n_f = min(metal_traj.shape[0], opencl_traj.shape[0])
    if n_frames is not None:
        n_f = min(n_f, n_frames)
    me = metal_traj[:n_f, :n_e, :]
    op = opencl_traj[:n_f, :n_e, :]
    diff = me - op
    return float((diff ** 2).sum(axis=-1).mean())


def run_metal(scenario, params, n_steps, chunk, tag):
    out = f"/tmp/sgd_worm_{tag}.txt"
    cmd = [
        "python3", str(HERE / "dump_metal_trajectory.py"),
        "--scenario", scenario,
        "--steps", str(n_steps),
        "--chunk", str(chunk),
        "--rho-rest", f"{params['rho_rest']:.6e}",
        "--visc-pair-coef", f"{params['visc']:.6e}",
        "--spring-k", f"{params['spring_K']:.6f}",
        "--anchor-k", f"{params['anchor_k']:.6f}",
        "--alpha-dist", f"{params['alpha_dist']:.6e}",
        "--floor-y", str(params.get('floor_y', 0)),
        "--vel-clamp", "0",
        "--no-box-clamp",
        "--out", out,
    ]
    if not params.get('membranes', True):
        cmd.append("--no-membranes")
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
    if r.returncode != 0:
        raise RuntimeError(f"sim failed: {r.stderr[-500:]}")
    return out


def fd_grad(scenario, base_params, target_traj, n_e, key, eps=0.05, **run_kwargs):
    """Centered FD on log(param[key])."""
    p_p = dict(base_params); p_p[key] = base_params[key] * np.exp(eps)
    p_m = dict(base_params); p_m[key] = base_params[key] * np.exp(-eps)
    out_p = run_metal(scenario, p_p, tag=f'fd_p_{key}', **run_kwargs)
    out_m = run_metal(scenario, p_m, tag=f'fd_m_{key}', **run_kwargs)
    L_p = loss_vs_opencl(parse_traj(out_p)[0], target_traj, n_e)
    L_m = loss_vs_opencl(parse_traj(out_m)[0], target_traj, n_e)
    return (L_p - L_m) / (2 * eps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default='worm_alone_half_resolution')
    ap.add_argument('--opencl-traj', required=True,
                    help='OpenCL position_buffer.txt')
    ap.add_argument('--n-steps', type=int, default=250,
                    help='5 ms at dt=2e-5')
    ap.add_argument('--chunk', type=int, default=25)
    ap.add_argument('--max-iters', type=int, default=15)
    ap.add_argument('--lr', type=float, default=0.3)
    ap.add_argument('--fd-eps', type=float, default=0.05)
    ap.add_argument('--init-spring-k', type=float, default=5500)
    ap.add_argument('--init-anchor-k', type=float, default=1.0,
                    help='No anchors in this config but kernel still reads it')
    ap.add_argument('--init-rho-rest', type=float, default=1000.0)
    ap.add_argument('--init-visc', type=float, default=5e-5)
    ap.add_argument('--init-alpha-dist', type=float, default=3.3e-9)
    ap.add_argument('--trainable', default='spring_K,rho_rest,visc',
                    help='comma-separated subset of {spring_K,anchor_k,rho_rest,visc,alpha_dist}')
    ap.add_argument('--history-out', default='/tmp/sgd_worm_history.json')
    args = ap.parse_args()

    target_traj, n_e_t, n_l_t, _ = parse_traj(args.opencl_traj)
    print(f'OpenCL target: shape={target_traj.shape}, n_e={n_e_t}')

    params = {
        'spring_K': args.init_spring_k,
        'anchor_k': args.init_anchor_k,
        'rho_rest': args.init_rho_rest,
        'visc': args.init_visc,
        'alpha_dist': args.init_alpha_dist,
        'floor_y': 0,
        'membranes': True,
    }
    trainable = args.trainable.split(',')
    print(f'Trainable: {trainable}')

    history = []
    best_L = float('inf'); best_params = dict(params)

    for it in range(args.max_iters):
        t0 = time.perf_counter()
        try:
            out = run_metal(args.scenario, params, args.n_steps, args.chunk,
                             tag=f'i{it}_base')
            metal_traj, n_e, n_l, _ = parse_traj(out)
        except RuntimeError as e:
            print(f'[step {it}] sim failed: {e}'); break

        n_f = min(metal_traj.shape[0], target_traj.shape[0])
        if n_f == 0:
            print(f'[step {it}] no overlapping frames'); break

        L = loss_vs_opencl(metal_traj, target_traj, n_e)
        if not np.isfinite(L):
            print(f'[step {it}] L diverged'); break

        if L < best_L:
            best_L = L; best_params = dict(params); star = ' *'
        else:
            star = ''

        line = (f'[step {it}] L={L:.4e}{star}  ' +
                ' '.join(f'{k}={params[k]:.4g}' for k in trainable))
        print(line)

        # FD over each trainable
        grads = {}
        for k in trainable:
            try:
                g = fd_grad(args.scenario, params, target_traj, n_e, k,
                             eps=args.fd_eps, n_steps=args.n_steps,
                             chunk=args.chunk)
                grads[k] = g
            except RuntimeError as e:
                print(f'  FD {k} failed: {e}'); grads[k] = 0
        elapsed = time.perf_counter() - t0
        print(f'  grads={grads}  ({elapsed:.1f}s)')

        history.append({'step': it, 'L': L, 'params': dict(params), 'grads': grads})
        # Step in log space (clipped)
        for k in trainable:
            g = max(-10, min(10, grads.get(k, 0)))
            params[k] = float(params[k] * np.exp(-args.lr * g))
            params[k] = max(params[k], 1e-12)

        if L < 1e-3:
            print(f'\n[CONVERGED at L={L:.2e}]')
            break

    print(f'\n=== Best === L={best_L:.4e}')
    for k, v in best_params.items():
        print(f'  {k} = {v}')
    with open(args.history_out, 'w') as f:
        json.dump(history, f, indent=2, default=str)
    print(f'History → {args.history_out}')


if __name__ == '__main__':
    main()
