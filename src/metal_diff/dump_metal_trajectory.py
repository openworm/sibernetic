"""Run xpbd_step in chunks and write Sibernetic-format position_buffer.txt.

Each chunk runs `chunk_steps` XPBD steps; we capture the final position
snapshot after each chunk. Resulting frames are concatenated into a
Sibernetic-format file that `render_movie.py` consumes.

Format (matches scripts/measure_cube_stability.py header parsing):
  line 0..5   box bounds (xmin, xmax, ymin, ymax, zmin, zmax)
  line 6      numOfElasticP
  line 7      numOfLiquidP
  line 8      numOfBoundaryP
  line 9      time_step
  line 10     log_step
  line 11+    particles concatenated (elastic → liquid → boundary)
              4 cols: x y z type (1.x liquid, 2.x elastic, 3.x boundary)

Usage:
    python3 src/metal_diff/dump_metal_trajectory.py \\
        --steps 150000 --chunk 150 --dt 2e-5 \\
        --rho-rest 2.5e-13 --spring-k 1000 --visc-pair-coef 5e-5 \\
        --out /tmp/metal_position_buffer.txt
"""
import argparse
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def run_chunk(args, info, pos_in, vel_in, n_steps):
    """Run xpbd_step for n_steps, return (final_pos, final_vel)."""
    pos_in.astype(np.float32).tofile('/tmp/dmt_pos_in.bin')
    vel_in.astype(np.float32).tofile('/tmp/dmt_vel_in.bin')
    cmd = [BIN, "xpbd_step",
           str(info['n_active']), str(info['n_static']),
           str(args.h), str(args.mass), f"{args.rho_rest:.15e}",
           str(args.dt), str(args.gravity_y), str(args.floor_y),
           str(args.alpha_dens), str(args.n_iters),
           '/tmp/dmt_pos_in.bin', '/tmp/dmt_vel_in.bin',
           info['paths']['pos_static'],
           str(info['n_bonds']), info['paths']['bonds'],
           str(args.alpha_dist), str(n_steps),
           f"{args.sim_scale:.15e}",
           f"{args.visc_pair_coef:.15e}",
           f"{args.spring_k:.15e}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"xpbd_step failed: {r.stderr[:300]}")
    pos_out = np.fromfile('/tmp/xpbd_pos_out.bin',
                          dtype=np.float32).reshape(-1, 3)
    vel_out = np.fromfile('/tmp/xpbd_vel_out.bin',
                          dtype=np.float32).reshape(-1, 3)
    return pos_out, vel_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default='demo1')
    ap.add_argument('--steps', type=int, default=150000,
                    help="total XPBD steps (3s sim at dt=2e-5)")
    ap.add_argument('--chunk', type=int, default=150,
                    help="steps per snapshot — matches Sibernetic logstep")
    ap.add_argument('--dense-steps', type=int, default=0,
                    help="if > 0, run first --dense-steps with --dense-chunk "
                         "(e.g. --dense-steps 5000 --dense-chunk 5 for descent capture)")
    ap.add_argument('--dense-chunk', type=int, default=5,
                    help="chunk size for the dense descent phase (default 5)")
    ap.add_argument('--dt', type=float, default=2e-5)
    ap.add_argument('--gravity-y', type=float, default=-9.8)
    ap.add_argument('--n-iters', type=int, default=3)
    ap.add_argument('--h', type=float, default=3.34)
    ap.add_argument('--mass', type=float, default=2e-12)
    ap.add_argument('--sim-scale', type=float, default=7.4e-6)
    ap.add_argument('--rho-rest', type=float, default=1000.0,
                    help="rest density (kg/m³); Sibernetic uses 1000 (water).")
    ap.add_argument('--alpha-dens', type=float, default=1e-3)
    ap.add_argument('--alpha-dist', type=float, default=3.3e-9)
    ap.add_argument('--visc-pair-coef', type=float, default=5e-5)
    ap.add_argument('--spring-k', type=float, default=1000.0)
    ap.add_argument('--floor-y', type=float, default=0.0)
    ap.add_argument('--out', default='/tmp/metal_position_buffer.txt')
    args = ap.parse_args()

    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers, parse_config

    info = load_to_metal_buffers(args.scenario,
                                  out_dir='/tmp/dmt_metal',
                                  h=args.h)
    cfg = parse_config(os.path.join("configuration", args.scenario))
    box = cfg['box']
    pos_all = cfg['pos']                       # [N, 4] x,y,z,type — full original ordering
    n_total = len(pos_all)

    types = pos_all[:, 3]
    n_elastic = int(((types >= 2.0) & (types < 3.0)).sum())
    n_liquid  = int(((types >= 1.0) & (types < 2.0)).sum())
    n_static  = int((types >= 3.0).sum())

    # Active buffer ordering: elastic (0..n_elastic-1) + liquid (n_elastic..n_active-1).
    # Static buffer: boundary in original order.
    pos_active = pos_all[types < 3.0, :3].astype(np.float32)
    vel_active = np.zeros_like(pos_active)
    pos_static = pos_all[types >= 3.0, :3].astype(np.float32)
    static_types = types[types >= 3.0]
    active_types = types[types < 3.0]

    # Build a schedule of chunks. If --dense-steps > 0, the first
    # --dense-steps run at --dense-chunk, then the rest at --chunk.
    schedule = []
    if args.dense_steps > 0:
        n_dense = args.dense_steps // args.dense_chunk
        schedule.extend([args.dense_chunk] * n_dense)
        rem = args.steps - args.dense_steps
        if rem > 0:
            n_sparse = rem // args.chunk
            schedule.extend([args.chunk] * n_sparse)
    else:
        n_chunks = args.steps // args.chunk
        schedule.extend([args.chunk] * n_chunks)

    print(f"Dumping Metal trajectory: {args.steps} steps in {len(schedule)} chunks")
    if args.dense_steps > 0:
        print(f"  dense: first {args.dense_steps} steps in chunks of {args.dense_chunk}")
        print(f"  sparse: rest in chunks of {args.chunk}")
    else:
        print(f"  uniform: chunks of {args.chunk}")
    print(f"  rho={args.rho_rest:.3e} K={args.spring_k:.1f} v={args.visc_pair_coef:.3e}")
    print(f"  Output: {args.out}")

    snapshots = [pos_active.copy()]   # frame 0
    chunk_sizes_used = []             # per-snapshot chunk size, for header
    t0 = time.perf_counter()
    cur_pos, cur_vel = pos_active.copy(), vel_active.copy()
    for k, n_steps in enumerate(schedule):
        cur_pos, cur_vel = run_chunk(args, info, cur_pos, cur_vel, n_steps)
        snapshots.append(cur_pos.copy())
        chunk_sizes_used.append(n_steps)
        if (k + 1) % 100 == 0:
            elapsed = time.perf_counter() - t0
            mean_y = float(cur_pos[:n_elastic, 1].mean())
            print(f"  chunk {k+1}/{len(schedule)}: t={elapsed:.0f}s  cube mean_y={mean_y:.2f}")
    elapsed_total = time.perf_counter() - t0
    print(f"  done: {elapsed_total:.1f}s for {len(schedule)} chunks")

    # Write Sibernetic-format file.
    # If chunk sizes are variable (dense + sparse), we still write a single
    # log_step in the header (the smallest, dense_chunk) so legacy parsers
    # don't crash, but the real per-frame times go into a sidecar file at
    # <out>.times.txt — one float per line, the cumulative sim time of each
    # snapshot in seconds. Readers that want correct timing should prefer
    # the sidecar when it exists.
    header_log_step = args.dense_chunk if args.dense_steps > 0 else args.chunk
    with open(args.out, 'w') as f:
        for v in box: f.write(f"{v}\n")               # 6 lines: box bounds
        f.write(f"{n_elastic}\n")
        f.write(f"{n_liquid}\n")
        f.write(f"{n_static}\n")
        f.write(f"{args.dt}\n")
        f.write(f"{header_log_step}\n")               # log_step (best-effort)
        for snap in snapshots:
            for i in range(len(snap)):
                t = active_types[i]
                f.write(f"{snap[i, 0]} {snap[i, 1]} {snap[i, 2]} {t}\n")
            for i in range(len(pos_static)):
                t = static_types[i]
                f.write(f"{pos_static[i, 0]} {pos_static[i, 1]} {pos_static[i, 2]} {t}\n")

    # Sidecar with explicit per-snapshot sim times (handles variable chunks).
    sidecar = args.out + ".times.txt"
    cumulative_steps = 0
    with open(sidecar, 'w') as f:
        f.write(f"{0.0}\n")  # frame 0
        for n_steps in chunk_sizes_used:
            cumulative_steps += n_steps
            f.write(f"{cumulative_steps * args.dt}\n")
    print(f"Wrote {args.out} ({os.path.getsize(args.out)/1e6:.1f} MB, {len(snapshots)} frames)")
    print(f"Wrote {sidecar} ({len(snapshots)} frame times)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
