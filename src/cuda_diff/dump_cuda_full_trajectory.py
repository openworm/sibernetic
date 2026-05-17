"""Run xpbd_full_fwd once for K steps, dump Sibernetic position_buffer.txt.

Parallel to dump_cuda_trajectory.py but uses the differentiable
xpbd_full_fwd op (single CLI call, K steps in one shot, full trajectory
tape on disk) rather than the imperative xpbd_step path (chunked).

This is needed for cross-backend regression vs OpenCL Sibernetic: the
OpenCL `position_buffer.txt` is the input to scripts/measure_cube_stability.py,
so to compare CUDA against OpenCL we need our K-step XPBD trajectory in
the same line-based file format.

state_out.bin layout (from src/cuda_diff/ops_xpbd_full.cu around l.397):
    K  * per_step_floats float32  (per-step intermediate tape — unused here)
    (K+1) * n_active * 3 float32  (trajectory positions, frame 0 = pos0)
    n_active * 3 float32          (final velocity)

  per_step_floats = n_active * (3+3+1+3+1 + extra)
    extra = (use_pair ? 1 : 0) + (use_floor ? 1 : 0)

Format (matches scripts/measure_cube_stability.py header parsing):
  line 0..5   box bounds
  line 6      numOfElasticP
  line 7      numOfLiquidP
  line 8      numOfBoundaryP
  line 9      time_step
  line 10     log_step
  line 11+    particles concatenated (elastic -> liquid -> boundary),
              4 cols: x y z type (1.x liquid, 2.x elastic, 3.x boundary)

Usage:
    python src/cuda_diff/dump_cuda_full_trajectory.py \
        --scenario demo1 --K 50000 --logstep 500 --dt 2e-5 \
        --out /tmp/cuda_full_position_buffer.txt
"""
import argparse
import os
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE,
                   "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default='demo1')
    ap.add_argument('--K', type=int, default=50000,
                    help="total XPBD steps (default 50000 ~ 1.0s at dt=2e-5)")
    ap.add_argument('--logstep', type=int, default=500,
                    help="write every N-th trajectory frame to output "
                         "(default 500). Matches Sibernetic's logstep.")
    ap.add_argument('--dt', type=float, default=2e-5)
    ap.add_argument('--gravity-y', type=float, default=-9.8)
    ap.add_argument('--h', type=float, default=3.34)
    ap.add_argument('--mass', type=float, default=2e-12)
    ap.add_argument('--sim-scale', type=float, default=7.4e-6)
    ap.add_argument('--rho-rest', type=float, default=2.5e-13)
    ap.add_argument('--alpha-dens', type=float, default=1e-3)
    ap.add_argument('--visc-pair-coef', type=float, default=5e-5)
    ap.add_argument('--spring-k', type=float, default=1000.0)
    ap.add_argument('--floor-y', type=float, default=0.0)
    ap.add_argument('--restitution', type=float, default=0.0)
    ap.add_argument('--out',
                    default=os.path.join(tempfile.gettempdir(),
                                         "cuda_full_position_buffer.txt"))
    ap.add_argument('--keep-state', action='store_true',
                    help="don't delete state_out.bin after parsing")
    args = ap.parse_args()

    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers, parse_config

    tmp_dir = os.path.join(tempfile.gettempdir(), "dcft_cuda")
    os.makedirs(tmp_dir, exist_ok=True)
    info = load_to_metal_buffers(args.scenario, out_dir=tmp_dir, h=args.h)
    cfg = parse_config(os.path.join("configuration", args.scenario))
    box = cfg['box']
    pos_all = cfg['pos']
    types = pos_all[:, 3]
    n_elastic = int(((types >= 2.0) & (types < 3.0)).sum())
    n_liquid = int(((types >= 1.0) & (types < 2.0)).sum())
    n_static = int((types >= 3.0).sum())
    n_active = info['n_active']

    active_types = types[types < 3.0]
    static_types = types[types >= 3.0]
    pos_static = pos_all[types >= 3.0, :3].astype(np.float32)

    state_out = os.path.join(tmp_dir, "state_out.bin")

    cmd = [BIN, "xpbd_full_fwd",
           str(n_active), str(n_static), str(args.K),
           str(args.h), str(args.mass), f"{args.rho_rest:.15e}",
           str(args.dt), str(args.gravity_y), str(args.alpha_dens),
           info['paths']['pos_active'], info['paths']['vel_active'],
           info['paths']['pos_static'], state_out,
           f"{args.sim_scale:.15e}",
           f"{args.visc_pair_coef:.15e}",
           f"{args.spring_k:.15e}",
           info['paths']['bonds'],
           str(args.floor_y),
           f"{args.restitution:.6f}"]

    print(f"Running xpbd_full_fwd: K={args.K} steps "
          f"({args.K * args.dt:.3f}s sim time)")
    print(f"  n_active={n_active} n_static={n_static} "
          f"n_bonds={info['n_bonds']}")
    print(f"  rho={args.rho_rest:.3e} K={args.spring_k:.1f} "
          f"visc={args.visc_pair_coef:.3e}")
    print(f"  state_out: {state_out}")

    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if r.returncode != 0:
        print(f"sib_cuda failed (exit {r.returncode}):", file=sys.stderr)
        print(r.stderr, file=sys.stderr)
        return 1
    print(f"  done in {elapsed:.1f}s")
    if r.stderr.strip():
        print(f"  stderr: {r.stderr.strip()[:200]}")

    # Parse state_out.bin: skip per-step tape, read (K+1) trajectory frames.
    # extra_per_step = (use_pair?1:0) + (use_floor?1:0).
    # We always pass visc_pair_coef>0 and floor_y so extra=2 here.
    use_pair = args.visc_pair_coef != 0.0
    use_floor = True   # we always pass floor_y
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step_floats = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    state_offset_bytes = args.K * per_step_floats * 4
    traj_floats = (args.K + 1) * n_active * 3

    with open(state_out, 'rb') as f:
        f.seek(state_offset_bytes)
        traj = np.fromfile(f, dtype=np.float32, count=traj_floats)
    if traj.size != traj_floats:
        print(f"unexpected trajectory size: got {traj.size}, "
              f"expected {traj_floats}", file=sys.stderr)
        return 1
    traj = traj.reshape(args.K + 1, n_active, 3)

    # Subsample frames at logstep granularity (frame 0 + every logstep-th).
    frame_indices = list(range(0, args.K + 1, args.logstep))
    if frame_indices[-1] != args.K:
        frame_indices.append(args.K)
    n_frames = len(frame_indices)
    print(f"  trajectory: {args.K + 1} raw frames -> "
          f"{n_frames} logged frames (logstep={args.logstep})")

    # Write Sibernetic format. Active buffer ordering is elastic -> liquid
    # (load_config splits by type < 3.0 with elastic first because they
    # appear first in the config).
    #
    # Match OpenCL Sibernetic's exact format: boundary particles are written
    # ONLY in the first logged frame; subsequent frames contain only
    # elastic+liquid (n_active rows). See sibernetic/src/owHelper.cpp
    # loadConfigurationToFile() — the `firstIteration` branch writes
    # boundary; later frames skip BOUNDARY_PARTICLE rows.
    with open(args.out, 'w') as f:
        for v in box:
            f.write(f"{v}\n")
        f.write(f"{n_elastic}\n")
        f.write(f"{n_liquid}\n")
        f.write(f"{n_static}\n")
        f.write(f"{args.dt}\n")
        f.write(f"{args.logstep}\n")
        for fi_idx, fi in enumerate(frame_indices):
            snap = traj[fi]
            for i in range(n_active):
                t = active_types[i]
                f.write(f"{snap[i, 0]} {snap[i, 1]} {snap[i, 2]} {t}\n")
            if fi_idx == 0:
                for i in range(n_static):
                    t = static_types[i]
                    f.write(f"{pos_static[i, 0]} {pos_static[i, 1]} "
                            f"{pos_static[i, 2]} {t}\n")
    sz_mb = os.path.getsize(args.out) / 1e6
    print(f"Wrote {args.out} ({sz_mb:.1f} MB, {n_frames} frames)")

    if not args.keep_state:
        try:
            os.remove(state_out)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
