"""Run demo1 (the same scenario OpenCL runs) on the hand-written Metal substrate.

This proves we can use Sibernetic's existing configuration files
instead of hardcoded toy cubes — the user's stated requirement after
M8 landed.

Pipeline:
  1. Parse configuration/demo1 via load_config.py
  2. Hand the binary buffers to sib_metal xpbd_step
  3. Run K simulation steps with full XPBD physics
     (predict + density constraint + elastic bonds + floor + update_vel)
  4. Report cube extent, mean Y, etc. as a sanity check

Initial physics tuning is approximate — we don't yet auto-derive
rho_rest, h, etc. from the config file. These are passed as CLI flags.

Run:
    .venv/bin/python src/metal_diff/run_demo1_via_metal.py [--steps N] [--dt DT]
"""
import argparse
import os
import subprocess
import time

import numpy as np

from load_config import load_to_metal_buffers

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default='demo1')
    ap.add_argument('--steps', type=int, default=100,
                    help="K — number of XPBD timesteps")
    ap.add_argument('--dt', type=float, default=2e-5,
                    help="dt — same as Sibernetic OpenCL default")
    ap.add_argument('--gravity-y', type=float, default=-9.8)
    ap.add_argument('--n-iters', type=int, default=3,
                    help="XPBD inner iteration count")
    # Defaults below reflect Sibernetic's `owPhysicsConstant.h`:
    #   h = 3.34 (smoothing radius in particle units)
    #   mass = 20.00e-13 kg (per-particle mass)
    #   sim_scale = 7.4e-6 (m / particle-unit; depends on mass)
    #   rho0 = 1000 kg/m^3   →   rho_rest_units = rho0 * sim_scale^3 = 4.05e-13
    # The user can override any of these to reproduce a different
    # Sibernetic configuration.
    ap.add_argument('--h', type=float, default=3.34,
                    help="SPH smoothing radius (Sibernetic h, particle units)")
    ap.add_argument('--mass', type=float, default=2.0e-12,
                    help="per-particle mass in kg (Sibernetic 20.00e-13)")
    ap.add_argument('--sim-scale', type=float, default=7.4e-6,
                    help="m / particle-unit. Sibernetic ≈ 7.4e-6 for "
                         "20.00e-13 kg particle mass")
    ap.add_argument('--rho-rest', type=float, default=4.05e-13,
                    help="rest density (mass/unit^3). For Sibernetic-style "
                         "units this is rho0_kg_per_m3 * sim_scale^3.")
    ap.add_argument('--alpha-dens', type=float, default=1e-3,
                    help="density compliance (smaller = stiffer)")
    ap.add_argument('--alpha-dist', type=float, default=3.3e-9,
                    help="bond compliance (smaller = stiffer; "
                         "Sibernetic elasticityCoefficient = 1/3.3e-9)")
    ap.add_argument('--floor-y', type=float, default=0.0)
    args = ap.parse_args()

    # Step 1: parse config + write binary buffers
    info = load_to_metal_buffers(args.scenario, out_dir="/tmp/demo1_metal",
                                  h=args.h)
    print(f"Loaded configuration/{args.scenario}:")
    print(f"  active = {info['n_active']} (elastic={info['n_elastic']}, "
          f"liquid={info['n_liquid']})")
    print(f"  static (boundary) = {info['n_static']}")
    print(f"  bonds = {info['n_bonds']}")
    print()

    # Step 2: load initial state (for measuring deltas)
    pos_init = np.fromfile(info['paths']['pos_active'], dtype=np.float32) \
                 .reshape(info['n_active'], 3)
    init_mean_y = float(pos_init[:, 1].mean())
    init_min_y = float(pos_init[:, 1].min())
    init_max_y = float(pos_init[:, 1].max())
    print(f"Initial state:")
    print(f"  mean_y={init_mean_y:.3f}  min_y={init_min_y:.3f}  max_y={init_max_y:.3f}")
    print()

    # Step 3: invoke sib_metal xpbd_step with the loaded buffers
    cmd = [
        BINARY, "xpbd_step",
        str(info['n_active']), str(info['n_static']),
        str(args.h), str(args.mass), str(args.rho_rest),
        str(args.dt), str(args.gravity_y), str(args.floor_y),
        str(args.alpha_dens), str(args.n_iters),
        info['paths']['pos_active'],
        info['paths']['vel_active'],
        info['paths']['pos_static'],
        str(info['n_bonds']), info['paths']['bonds'],
        str(args.alpha_dist),
        str(args.steps),
        str(args.sim_scale),
    ]
    print(f"Running xpbd_step for {args.steps} steps (dt={args.dt}, "
          f"sim_time={args.steps * args.dt:.4f}s)...")
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    print(f"  Wall: {wall:.2f}s")
    if result.returncode != 0:
        print(f"  FAILED: {result.stderr}")
        return 1
    # The xpbd_step driver prints per-step ms on stderr.
    for line in result.stderr.strip().split('\n'):
        if line:
            print(f"  {line}")
    print()

    # Step 4: read output positions
    pos_out = np.fromfile("/tmp/xpbd_pos_out.bin", dtype=np.float32) \
                .reshape(info['n_active'], 3)
    final_mean_y = float(pos_out[:, 1].mean())
    final_min_y = float(pos_out[:, 1].min())
    final_max_y = float(pos_out[:, 1].max())

    # Cube extent (look at elastic particles only, indices 0..n_elastic-1)
    pos_elastic_init = pos_init[:info['n_elastic']]
    pos_elastic_final = pos_out[:info['n_elastic']]
    init_extent = pos_elastic_init.max(axis=0) - pos_elastic_init.min(axis=0)
    final_extent = pos_elastic_final.max(axis=0) - pos_elastic_final.min(axis=0)

    print(f"Final state:")
    print(f"  mean_y={final_mean_y:.3f}  min_y={final_min_y:.3f}  max_y={final_max_y:.3f}")
    print(f"  Δmean_y = {final_mean_y - init_mean_y:+.3f}")
    print(f"  Cube (elastic only) extent: "
          f"init=({init_extent[0]:.2f}, {init_extent[1]:.2f}, {init_extent[2]:.2f})  "
          f"final=({final_extent[0]:.2f}, {final_extent[1]:.2f}, {final_extent[2]:.2f})")
    integrity = (final_extent / init_extent)
    print(f"  Cube integrity ratio: ({integrity[0]:.2f}, {integrity[1]:.2f}, {integrity[2]:.2f})")

    # Sanity checks
    has_nan = bool(np.isnan(pos_out).any())
    bounded = bool(np.abs(pos_out).max() < 1000.0)
    print()
    print(f"  no NaN:        {'[OK]' if not has_nan else '[FAIL]'}")
    print(f"  positions bounded (<1000): {'[OK]' if bounded else '[FAIL]'}")
    return 0 if (not has_nan and bounded) else 1


if __name__ == "__main__":
    raise SystemExit(main())
