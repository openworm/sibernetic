"""Run xpbd_step in chunks and write a Sibernetic-format position_buffer.txt.

Each chunk runs `chunk_steps` XPBD steps via the CUDA `xpbd_step`
orchestrator; we capture the final position snapshot after each chunk.
Resulting frames are concatenated into a Sibernetic-format file that
`render_movie.py` consumes.

CUDA `xpbd_step` covers the cube-drop scope (density + distance
constraints + floor). Pair forces, spring forces, membranes, and
anchors are NOT exposed by this op — they live behind `xpbd_full_fwd`,
which `dump_cuda_full_trajectory.py` drives. Use that script when you
need a worm-physics or fluid+visc trajectory.

Format (matches scripts/measure_cube_stability.py header parsing):
  line 0..5   box bounds (xmin, xmax, ymin, ymax, zmin, zmax)
  line 6      numOfElasticP
  line 7      numOfLiquidP
  line 8      numOfBoundaryP
  line 9      time_step
  line 10     log_step
  line 11+    particles concatenated (elastic -> liquid -> boundary)
              4 cols: x y z type (1.x liquid, 2.x elastic, 3.x boundary)

Usage:
    python src/cuda_diff/dump_cuda_trajectory.py \\
        --steps 12500 --chunk 5 --dt 2e-5 \\
        --rho-rest 2.5e-13 \\
        --out cuda_position_buffer.txt
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
TMP = tempfile.gettempdir()


def _tmp(name):
    """Build a tempdir path; isolated per-script so concurrent dumps don't clash."""
    return os.path.join(TMP, "dump_cuda_xpbd", name)


def run_chunk(args, info, pos_in, vel_in, n_steps):
    """Run xpbd_step for n_steps, return (final_pos, final_vel).

    argv layout for CUDA xpbd_step (see ops_xpbd_step.cu:8):
      argv[2..10]   n_active, n_static, h, mass, rho_rest, dt,
                    gravity_y, floor_y, alpha_density
      argv[11..12]  n_iters, n_steps
      argv[13..15]  pos_in.bin, vel_in.bin, pos_static.bin
      argv[16..18]  n_bonds, bonds.bin, alpha_dist
      argv[19..20]  pos_out.bin, vel_out.bin
      argv[21..26]  [traj_out] [traj_every] [sim_scale] [vel_clamp]
                    [restitution] [skip_density]
    """
    p_pos_in = _tmp("pos_in.bin")
    p_vel_in = _tmp("vel_in.bin")
    p_pos_out = _tmp("pos_out.bin")
    p_vel_out = _tmp("vel_out.bin")
    pos_in.astype(np.float32).tofile(p_pos_in)
    vel_in.astype(np.float32).tofile(p_vel_in)

    cmd = [BIN, "xpbd_step",
           str(info['n_active']), str(info['n_static']),
           str(args.h), str(args.mass), f"{args.rho_rest:.15e}",
           str(args.dt), str(args.gravity_y), str(args.floor_y),
           str(args.alpha_dens),
           str(args.n_iters), str(n_steps),
           p_pos_in, p_vel_in, info['paths']['pos_static'],
           str(info['n_bonds']), info['paths']['bonds'],
           str(args.alpha_dist),
           p_pos_out, p_vel_out]

    # Optional positional args. The CUDA xpbd_step CLI consumes them in
    # order: traj_out, traj_every, sim_scale, vel_clamp, restitution,
    # skip_density. Each is independently optional, but you can't pass
    # a later one without filling the earlier slots, so fill with
    # safe defaults whenever a later non-default value is needed.
    needs_optional = (args.sim_scale != 1.0 or args.vel_clamp > 0
                      or args.restitution != 0.0 or args.skip_density)
    if needs_optional:
        # traj_out + traj_every: this script captures frames itself between
        # chunks, so disable the C++ traj writer (empty path, 0 frequency).
        cmd.extend(["", "0",
                    f"{args.sim_scale:.15e}",
                    f"{args.vel_clamp:.15e}",
                    f"{args.restitution:.6f}",
                    "1" if args.skip_density else "0"])

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            f"sib_cuda xpbd_step failed (exit {r.returncode}):\n"
            f"  stderr: {r.stderr[:400]}")

    pos_out = np.fromfile(p_pos_out, dtype=np.float32).reshape(-1, 3)
    vel_out = np.fromfile(p_vel_out, dtype=np.float32).reshape(-1, 3)
    return pos_out, vel_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default='demo1')
    ap.add_argument('--steps', type=int, default=12500,
                    help="total XPBD steps (default 12500 = 0.25s at dt=2e-5)")
    ap.add_argument('--chunk', type=int, default=5,
                    help="steps per snapshot (matches Sibernetic logstep)")
    ap.add_argument('--dt', type=float, default=2e-5)
    ap.add_argument('--gravity-y', type=float, default=-9.8)
    ap.add_argument('--n-iters', type=int, default=3,
                    help="inner XPBD solver iterations per step")
    ap.add_argument('--h', type=float, default=3.34)
    ap.add_argument('--mass', type=float, default=2e-12)
    ap.add_argument('--sim-scale', type=float, default=1.0)
    ap.add_argument('--rho-rest', type=float, default=2.5e-13,
                    help="rest density (kg/sim_unit^3, NOT kg/m^3)")
    ap.add_argument('--alpha-dens', type=float, default=1e-3,
                    help="XPBD density compliance (smaller = stiffer)")
    ap.add_argument('--alpha-dist', type=float, default=3.3e-9,
                    help="XPBD distance-constraint compliance")
    ap.add_argument('--floor-y', type=float, default=0.0)
    ap.add_argument('--restitution', type=float, default=0.0,
                    help="floor restitution coefficient e in [0,1]")
    ap.add_argument('--vel-clamp', type=float, default=0.0,
                    help="max velocity magnitude (0 = disabled)")
    ap.add_argument('--skip-density', action='store_true',
                    help="bypass the SPH density chain (use for pure-bond cube tests)")
    ap.add_argument('--out',
                    default=os.path.join(TMP, "cuda_position_buffer.txt"))
    args = ap.parse_args()

    if not os.path.exists(BIN):
        print(f"ERROR: sib_cuda binary missing at {BIN}", file=sys.stderr)
        print(f"Build first: src/cuda_diff/build."
              f"{'bat' if platform.system() == 'Windows' else 'sh'}",
              file=sys.stderr)
        return 1

    os.makedirs(os.path.dirname(_tmp("dummy")), exist_ok=True)

    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers, parse_config

    info = load_to_metal_buffers(args.scenario,
                                  out_dir=os.path.join(TMP, "dump_cuda_buffers"),
                                  h=args.h)
    cfg = parse_config(os.path.join("configuration", args.scenario))
    box = cfg['box']
    pos_all = cfg['pos']
    types = pos_all[:, 3]
    n_elastic = int(((types >= 2.0) & (types < 3.0)).sum())
    n_liquid = int(((types >= 1.0) & (types < 2.0)).sum())
    n_static = int((types >= 3.0).sum())

    # Active buffer ordering matches load_to_metal_buffers: elastic first,
    # then liquid, in original config order.
    pos_active = pos_all[types < 3.0, :3].astype(np.float32)
    vel_all = cfg['vel'][:, :3].astype(np.float32)
    vel_active = vel_all[types < 3.0].astype(np.float32)
    pos_static = pos_all[types >= 3.0, :3].astype(np.float32)
    static_types = types[types >= 3.0]
    active_types = types[types < 3.0]

    n_chunks = args.steps // args.chunk
    print(f"Dumping CUDA xpbd_step trajectory: {args.steps} steps "
          f"in {n_chunks} chunks of {args.chunk}")
    print(f"  rho={args.rho_rest:.3e}  alpha_dens={args.alpha_dens:.3e}  "
          f"alpha_dist={args.alpha_dist:.3e}")
    print(f"  Output: {args.out}")

    snapshots = [pos_active.copy()]
    t0 = time.perf_counter()
    cur_pos, cur_vel = pos_active.copy(), vel_active.copy()
    for k in range(n_chunks):
        cur_pos, cur_vel = run_chunk(args, info, cur_pos, cur_vel, args.chunk)
        snapshots.append(cur_pos.copy())
        if (k + 1) % 100 == 0:
            elapsed = time.perf_counter() - t0
            mean_y = float(cur_pos[:n_elastic, 1].mean()) if n_elastic else float('nan')
            print(f"  chunk {k+1}/{n_chunks}: t={elapsed:.0f}s  "
                  f"cube mean_y={mean_y:.2f}")
    elapsed_total = time.perf_counter() - t0
    print(f"  done: {elapsed_total:.1f}s for {n_chunks} chunks")

    with open(args.out, 'w') as f:
        for v in box:
            f.write(f"{v}\n")
        f.write(f"{n_elastic}\n")
        f.write(f"{n_liquid}\n")
        f.write(f"{n_static}\n")
        f.write(f"{args.dt}\n")
        f.write(f"{args.chunk}\n")
        for snap in snapshots:
            for i in range(len(snap)):
                t = active_types[i]
                f.write(f"{snap[i, 0]} {snap[i, 1]} {snap[i, 2]} {t}\n")
            for i in range(len(pos_static)):
                t = static_types[i]
                f.write(f"{pos_static[i, 0]} {pos_static[i, 1]} "
                        f"{pos_static[i, 2]} {t}\n")
    print(f"Wrote {args.out} "
          f"({os.path.getsize(args.out)/1e6:.1f} MB, {len(snapshots)} frames)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
