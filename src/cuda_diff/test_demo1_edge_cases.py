"""V5 edge-case sweep on demo1 through the CUDA xpbd_full_fwd driver.

Runs the differentiable forward driver directly (no SGD wrapper) for K=1000
steps on the demo1 scenario, varying one physics parameter per variant. Each
run must:
  - exit 0,
  - produce a state_out.bin of the expected size,
  - contain no NaN/Inf anywhere in the trajectory frames,
  - have a finite mean Y in the final frame.

The 7 variants exercise the kernels that the standard demo1 SGD run doesn't
stress at their parameter extremes: stiff springs, no springs, high viscous
pair forces, no pair forces, no floor, and near-elastic floor restitution.

Pass criterion: all 7 variants meet the four checks above.
This script does NOT do gradient validation (FD parity is V4's job); the
goal here is shake-out for crashes, allocator failures, and NaN propagation.

Exits 0 iff all variants pass.
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
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP_DIR = os.path.join(tempfile.gettempdir(), "v5_edge_cases")


# Demo1 baseline params (match V3 SGD run + V2 forward-parity check).
SCENARIO = "demo1"
K = 1000
H = 3.34
MASS = 2e-12
RHO_REST = 2.5e-13
DT = 2e-5
G_Y = -9.8
ALPHA_DENS = 1e-3
SIM_SCALE = 7.4e-6
VISC_PAIR_COEF = 5e-5
SPRING_K = 1000.0
FLOOR_Y = 0.0
RESTITUTION = 0.0


def load_demo1_buffers():
    """Run load_config to prepare the binary buffers demo1 needs."""
    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers
    os.makedirs(TMP_DIR, exist_ok=True)
    info = load_to_metal_buffers(SCENARIO, out_dir=TMP_DIR, h=H)
    return info


def build_cmd(info, state_out_path, *, visc_pair_coef, spring_k, use_floor,
              floor_y, restitution):
    """Assemble the xpbd_full_fwd CLI exactly as Metal's positional contract.

    Arg slots beyond the required 14 are positional, so optional args before
    the slot we want must be filled (even if zero/dummy). Layout:
      15: sim_scale
      16: visc_pair_coef
      17: spring_K
      18: bonds.bin (only loaded when spring_K > 0; empty file is fine when
          spring_K = 0 but the slot is still positionally required to reach
          floor_y at slot 19).
      19: floor_y (presence ALONE enables the floor; passing 0.0 is "floor
          at y=0", not "no floor")
      20: restitution

    To express "no floor" we simply omit args 19 and 20 (and any slots after
    them). To express "no bonds" we omit args 17 onward -- but if we want
    floor support without springs, we still pass an empty bonds.bin so slot
    18 is filled.
    """
    cmd = [BINARY, "xpbd_full_fwd",
           str(info['n_active']), str(info['n_static']), str(K),
           str(H), str(MASS), f"{RHO_REST:.15e}",
           str(DT), str(G_Y), str(ALPHA_DENS),
           info['paths']['pos_active'], info['paths']['vel_active'],
           info['paths']['pos_static'], state_out_path,
           f"{SIM_SCALE:.15e}",
           f"{visc_pair_coef:.15e}",
           f"{spring_k:.15e}"]

    # Bonds slot is required positionally if we want to pass floor_y after.
    # When spring_K == 0 the file is not opened, but we still need a path.
    # When neither spring_K nor floor is needed, we can stop here.
    if spring_k != 0.0:
        cmd.append(info['paths']['bonds'])
    elif use_floor:
        # spring_K = 0 but we need floor_y: positionally fill the bonds slot
        # with a path the driver won't open (since spring_K == 0).
        dummy = os.path.join(TMP_DIR, "empty_bonds.bin")
        if not os.path.exists(dummy):
            open(dummy, "wb").close()
        cmd.append(dummy)

    if use_floor:
        cmd.append(str(floor_y))
        cmd.append(f"{restitution:.6f}")

    return cmd


def expected_state_size_bytes(n_active, K_, use_pair, use_floor):
    """Per ops_xpbd_full.cu line 132-138:
       per_step = n_active * (3+3+1+3+1 + (use_pair?1:0) + (use_floor?1:0))
       total = K * per_step + (K+1)*n_active*3 + n_active*3
    """
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    floats = K_ * per_step + (K_ + 1) * n_active * 3 + n_active * 3
    return 4 * floats


def parse_trajectory(state_path, n_active, K_, use_pair, use_floor):
    """Read just the trajectory portion (skip per-step tape)."""
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    state_offset_bytes = K_ * per_step * 4
    traj_floats = (K_ + 1) * n_active * 3
    with open(state_path, "rb") as f:
        f.seek(state_offset_bytes)
        traj = np.fromfile(f, dtype=np.float32, count=traj_floats)
    if traj.size != traj_floats:
        return None
    return traj.reshape(K_ + 1, n_active, 3)


def run_variant(name, info, *, visc_pair_coef, spring_k, use_floor,
                floor_y=0.0, restitution=0.0):
    """Run xpbd_full_fwd for one variant, return result dict."""
    state_out = os.path.join(TMP_DIR, f"state_{name}.bin")
    if os.path.exists(state_out):
        os.remove(state_out)
    cmd = build_cmd(info, state_out,
                    visc_pair_coef=visc_pair_coef, spring_k=spring_k,
                    use_floor=use_floor, floor_y=floor_y, restitution=restitution)
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    rc = r.returncode

    n_nans = -1
    final_mean_y = float('nan')
    size_ok = False
    actual_size = 0
    expected_size = 0

    use_pair = (visc_pair_coef != 0.0)
    if rc == 0 and os.path.exists(state_out):
        actual_size = os.path.getsize(state_out)
        expected_size = expected_state_size_bytes(
            info['n_active'], K, use_pair, use_floor)
        size_ok = (actual_size == expected_size)
        if size_ok:
            traj = parse_trajectory(state_out, info['n_active'], K,
                                    use_pair, use_floor)
            if traj is not None:
                n_nans = int(np.isnan(traj).sum() + np.isinf(traj).sum())
                if n_nans == 0:
                    final_mean_y = float(traj[-1, :, 1].mean())
                else:
                    final_mean_y = float('nan')

    # Cleanup state file once parsed -- these are ~700MB each for n_active=343.
    if os.path.exists(state_out):
        os.remove(state_out)

    return {
        'name': name,
        'rc': rc,
        'elapsed_s': elapsed,
        'size_ok': size_ok,
        'size_actual': actual_size,
        'size_expected': expected_size,
        'n_nans': n_nans,
        'final_mean_y': final_mean_y,
        'stderr_tail': r.stderr[-200:] if rc != 0 else '',
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true',
                    help="reduce K for smoke testing")
    args = ap.parse_args()

    if args.quick:
        global K
        K = 100

    if not os.path.exists(BINARY):
        print(f"ERROR: sib_cuda binary missing at {BINARY}", file=sys.stderr)
        return 1

    print(f"=== V5 edge-case sweep on demo1 (K={K}, ~{K * DT:.4f}s sim) ===")
    print(f"Binary: {BINARY}")
    print()

    info = load_demo1_buffers()
    print(f"demo1 loaded: n_active={info['n_active']} n_static={info['n_static']} "
          f"n_bonds={info['n_bonds']}")
    print()

    variants = [
        # (name, visc_pair_coef, spring_k, use_floor, floor_y, restitution)
        dict(name="1_baseline",       visc_pair_coef=VISC_PAIR_COEF, spring_k=SPRING_K, use_floor=True,  floor_y=0.0,  restitution=0.0),
        dict(name="2_stiff_spring",   visc_pair_coef=VISC_PAIR_COEF, spring_k=1e5,      use_floor=True,  floor_y=0.0,  restitution=0.0),
        dict(name="3_no_spring",      visc_pair_coef=VISC_PAIR_COEF, spring_k=0.0,      use_floor=True,  floor_y=0.0,  restitution=0.0),
        dict(name="4_high_visc",      visc_pair_coef=1e-2,           spring_k=SPRING_K, use_floor=True,  floor_y=0.0,  restitution=0.0),
        dict(name="5_no_pair_forces", visc_pair_coef=0.0,            spring_k=SPRING_K, use_floor=True,  floor_y=0.0,  restitution=0.0),
        dict(name="6_no_floor",       visc_pair_coef=VISC_PAIR_COEF, spring_k=SPRING_K, use_floor=False),
        dict(name="7_high_restitution", visc_pair_coef=VISC_PAIR_COEF, spring_k=SPRING_K, use_floor=True, floor_y=0.0, restitution=0.95),
    ]

    results = []
    for v in variants:
        print(f"  [{v['name']}] running...")
        res = run_variant(info=info, **v)
        results.append(res)

    # Table.
    print()
    print(f"{'variant':<22s} {'rc':>3s}  {'size_ok':>8s}  {'nans':>5s}  "
          f"{'final_y':>12s}  {'time_s':>7s}")
    print("-" * 72)
    for r in results:
        size_str = "ok" if r['size_ok'] else f"{r['size_actual']}!={r['size_expected']}"
        nan_str = str(r['n_nans']) if r['n_nans'] >= 0 else "-"
        y_str = f"{r['final_mean_y']:.4g}" if np.isfinite(r['final_mean_y']) else "NaN"
        print(f"{r['name']:<22s} {r['rc']:>3d}  {size_str:>8s}  {nan_str:>5s}  "
              f"{y_str:>12s}  {r['elapsed_s']:7.2f}")

    # Pass criterion per variant: rc=0, size_ok, n_nans=0, final_mean_y finite.
    fails = []
    for r in results:
        ok = (r['rc'] == 0 and r['size_ok'] and r['n_nans'] == 0
              and np.isfinite(r['final_mean_y']))
        if not ok:
            fails.append(r)

    print()
    if not fails:
        print(f"[OVERALL PASS] all {len(results)}/{len(results)} edge-case variants finite + finished")
        return 0
    print(f"[OVERALL FAIL] {len(fails)}/{len(results)} variants failed:")
    for r in fails:
        print(f"  - {r['name']}: rc={r['rc']} size_ok={r['size_ok']} "
              f"nans={r['n_nans']} final_y={r['final_mean_y']}")
        if r['stderr_tail']:
            print(f"    stderr: {r['stderr_tail']!r}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
